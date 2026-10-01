"""Protein and dataset enrichments for the browser: similar proteins, structures, findings.

``register(app, ctx)`` adds the routes and a Jinja global ``insight`` that page
templates query through small includes:

- **Similar proteins** come from the dataset's persistent FAISS index
  (``sharur build-vector-index``), opened read-only on first use; without an
  index the panel says how to build one.
- **Structures** are indexed only from explicit records: a canonical
  ``structures/<protein_id>.pdb`` (unsafe characters as ``_``), or JSON result
  files under ``structures/`` whose records name both ``protein_id`` and
  ``pdb_path``. Foldseek hits come from the same record or from result files
  keyed by the record's label. Models shorter than the protein are labeled as
  partial.
- **Findings** are read from ``findings.jsonl`` in the dataset directory and its
  immediate subdirectories, linked to the proteins, contigs and genomes they
  reference, and to system types named in their title. Verification queries
  rerun on a separate read-only cursor: one SELECT statement, file-reading
  functions rejected, with a timeout.
"""

from __future__ import annotations

import json
import logging
import os
import re
import threading
import time
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.responses import HTMLResponse, JSONResponse, PlainTextResponse

logger = logging.getLogger(__name__)

VERIFY_TIMEOUT_S = 10.0
MAX_RESULT_ROWS = 5
# table functions that read files, environment or secrets; statement types are checked separately
_FORBIDDEN_SQL = re.compile(
    r"\b(read_\w+|glob|sniff_csv|\w+_scan|query_table|query|getenv|duckdb_secrets|duckdb_settings)\s*\(",
    re.IGNORECASE,
)


def safe_name(protein_id: str) -> str:
    """Canonical structure file stem for a protein ID."""
    return re.sub(r"[^A-Za-z0-9_.-]", "_", protein_id)


# --------------------------------------------------------------------------- #
# Structures
# --------------------------------------------------------------------------- #


@dataclass
class StructureRecord:
    key: int
    protein_id: str
    path: Path
    label: str
    source: str
    plddt: float | None = None
    ptm: float | None = None
    model_length: int | None = None
    hits: list[dict[str, Any]] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {"key": self.key, "protein_id": self.protein_id, "label": self.label, "source": self.source,
                "plddt": self.plddt, "ptm": self.ptm, "model_length": self.model_length, "hits": self.hits}


def _normalize_hit(hit: dict[str, Any], database: str | None = None) -> dict[str, Any] | None:
    target = hit.get("target") or hit.get("target_id")
    if not target:
        return None
    identity = hit.get("identity", hit.get("seq_identity", hit.get("prob")))
    return {"database": hit.get("database") or hit.get("db") or database or "",
            "target": str(target), "description": hit.get("description") or hit.get("target_description") or "",
            "evalue": hit.get("evalue"), "identity": identity, "taxon": hit.get("taxon") or hit.get("taxName") or ""}


def _hits_from(value: Any) -> list[dict[str, Any]]:
    """Foldseek hits from a list of hit dicts or a {database: [hits]} mapping."""
    out: list[dict[str, Any]] = []
    if isinstance(value, list):
        out = [h for h in (_normalize_hit(x) for x in value if isinstance(x, dict)) if h]
    elif isinstance(value, dict):
        for database, hits in value.items():
            if isinstance(hits, list):
                out.extend(h for h in (_normalize_hit(x, database) for x in hits if isinstance(x, dict)) if h)
    return sorted(out, key=lambda h: (h["evalue"] if isinstance(h["evalue"], (int, float)) else 1e9))


def _records_in(data: Any, label: str | None = None) -> list[tuple[str | None, dict[str, Any]]]:
    """(label, record) pairs for dicts that carry a protein_id, at any depth."""
    found: list[tuple[str | None, dict[str, Any]]] = []
    if isinstance(data, dict):
        if "protein_id" in data:
            found.append((label or data.get("label"), data))
        for key, value in data.items():
            if isinstance(value, (dict, list)):
                found.extend(_records_in(value, key if isinstance(value, dict) else label))
    elif isinstance(data, list):
        for item in data:
            found.extend(_records_in(item, label))
    return found


def index_structures(dataset_dir: Path, known_proteins) -> dict[str, list[StructureRecord]]:
    """protein_id -> structure records, from explicit records only."""
    root = dataset_dir / "structures"
    by_protein: dict[str, list[StructureRecord]] = defaultdict(list)
    if not root.is_dir():
        return {}
    seen_paths: set[tuple[str, Path]] = set()
    json_docs: dict[Path, Any] = {}
    for path in sorted(root.glob("*.json")):
        try:
            json_docs[path] = json.loads(path.read_text())
        except (OSError, ValueError):
            continue
    # hits keyed by label in separate result files: {label: [hits]} or {label: {db: [hits]}}
    hits_by_label: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for doc in json_docs.values():
        if isinstance(doc, dict):
            for label, value in doc.items():
                hits = _hits_from(value)
                if hits:
                    hits_by_label[label].extend(hits)
    key = 0

    def add(protein_id: str, pdb: Path, label: str, source: str, record: dict[str, Any] | None = None):
        nonlocal key
        if (protein_id, pdb) in seen_paths:
            return
        seen_paths.add((protein_id, pdb))
        record = record or {}
        hits = _hits_from(record.get("foldseek_hits") or record.get("hits") or []) or hits_by_label.get(label, [])
        plddt = record.get("plddt", record.get("plddt_mean"))
        if isinstance(plddt, (int, float)) and 0 <= plddt <= 1:
            plddt = plddt * 100  # ESM3 reports pLDDT as a fraction; AlphaFold as 0-100
        length = record.get("length") or record.get("tip_length")
        by_protein[protein_id].append(StructureRecord(
            key, protein_id, pdb, label, source,
            float(plddt) if isinstance(plddt, (int, float)) else None,
            float(record["ptm"]) if isinstance(record.get("ptm"), (int, float)) else None,
            int(length) if isinstance(length, int) else None, hits[:15]))
        key += 1

    for path, doc in json_docs.items():
        for label, record in _records_in(doc):
            pdb_value = record.get("pdb_path")
            protein_id = record.get("protein_id")
            if not pdb_value or not isinstance(protein_id, str):
                continue
            candidates = [Path(pdb_value), dataset_dir.parent.parent / pdb_value, root / Path(pdb_value).name]
            pdb = next((c for c in candidates if c.is_file() and c.suffix.lower() in (".pdb", ".cif")), None)
            if pdb is not None:
                add(protein_id, pdb.resolve(), str(label or pdb.stem), path.name, record)
    # canonical file names: the protein ID, with characters outside [A-Za-z0-9_.-] as "_"
    stems = {p.stem: p for p in root.glob("*.pdb")}
    for stem, protein_ids in known_proteins(list(stems)).items():
        for protein_id in protein_ids:
            add(protein_id, stems[stem].resolve(), stem, "structures/" + stems[stem].name)
    return dict(by_protein)


# --------------------------------------------------------------------------- #
# Findings
# --------------------------------------------------------------------------- #


def load_findings(dataset_dir: Path) -> list[dict[str, Any]]:
    paths = [dataset_dir / "findings.jsonl"] + sorted(dataset_dir.glob("*/findings.jsonl"))
    findings, seen = [], set()
    for path in paths:
        if not path.is_file():
            continue
        try:
            lines = path.read_text().splitlines()
        except OSError:
            continue
        for n, line in enumerate(lines, 1):
            if not line.strip():
                continue
            try:
                row = json.loads(line)
            except ValueError:
                continue
            if not isinstance(row, dict):
                continue
            fid = str(row.get("id") or f"{path.parent.name}:{n}")
            if fid in seen:
                fid = f"{fid}@{path.parent.name}"
            seen.add(fid)
            row["_id"], row["_source"] = fid, str(path.relative_to(dataset_dir))
            findings.append(row)
    return findings


def _strings(value: Any, limit: int = 4000) -> list[str]:
    out: list[str] = []

    def walk(v: Any) -> None:
        if len(out) >= limit:
            return
        if isinstance(v, str):
            if 3 <= len(v) <= 300 and not any(ch.isspace() for ch in v):
                out.append(v)
        elif isinstance(v, dict):
            for item in v.values():
                walk(item)
        elif isinstance(v, (list, tuple)):
            for item in v:
                walk(item)

    walk(value)
    return out


def is_safe_select(sql: str) -> tuple[bool, str]:
    import duckdb  # noqa: PLC0415

    text = sql.strip().rstrip(";").strip()
    if not text:
        return False, "empty query"
    try:
        statements = duckdb.extract_statements(text)
    except Exception as exc:  # noqa: BLE001 - parser errors are reported to the user
        return False, f"not valid SQL ({type(exc).__name__})"
    if len(statements) != 1:
        return False, "only a single statement may run"
    if statements[0].type != duckdb.StatementType.SELECT:
        return False, f"only SELECT queries run here ({statements[0].type.name})"
    if not re.match(r"(?is)^\s*(--[^\n]*\n\s*)*(select|with|from|values|\()", text):
        return False, "only SELECT queries run here"
    if _FORBIDDEN_SQL.search(text):
        return False, "file, extension and settings functions are not allowed"
    return True, ""


def _compare(expected: Any, rows: list[tuple]) -> tuple[str, Any]:
    """('pass'|'fail'|'unchecked', observed)."""
    if len(rows) == 1 and len(rows[0]) == 1:
        observed: Any = rows[0][0]
    elif rows and all(len(r) == 1 for r in rows):
        observed = [r[0] for r in rows]
    else:
        observed = [list(r) for r in rows]
    if expected is None:
        return "unchecked", observed

    def same(a: Any, b: Any) -> bool:
        if isinstance(a, bool) or isinstance(b, bool):
            return a == b
        if isinstance(a, (int, float)) and isinstance(b, (int, float)):
            return abs(float(a) - float(b)) <= 1e-6 * max(1.0, abs(float(a)))
        if isinstance(a, list) and isinstance(b, list):
            return len(a) == len(b) and all(same(x, y) for x, y in zip(a, b, strict=True))
        return str(a) == str(b)

    if isinstance(expected, dict):
        return "unchecked", observed
    return ("pass" if same(expected, observed) else "fail"), observed


def _clip(value: Any) -> Any:
    if isinstance(value, str):
        return value if len(value) <= 80 else value[:77] + "…"
    if isinstance(value, list):
        return [_clip(v) for v in value[:MAX_RESULT_ROWS]] + (["…"] if len(value) > MAX_RESULT_ROWS else [])
    if isinstance(value, float):
        return round(value, 6)
    return value if isinstance(value, (int, bool)) or value is None else str(value)


class Insight:
    def __init__(self, ctx) -> None:
        self.ctx = ctx
        self.store, self.lock, self.catalog = ctx.store, ctx.lock, ctx.catalog
        self.dataset_dir = Path(ctx.db_path).resolve().parent
        self._similar = None
        self._similar_state: str | None = None
        self._similar_lock = threading.Lock()
        override = os.environ.get("SHARUR_BROWSE_EMBEDDINGS") or getattr(ctx, "embeddings_h5", None)
        self.embeddings_h5 = Path(override) if override else self.dataset_dir / "embeddings" / "protein_embeddings.h5"
        t0 = time.time()
        self.structures = index_structures(self.dataset_dir, self._proteins_for_stem)
        self.structure_files = {r.key: r for records in self.structures.values() for r in records}
        self.findings = load_findings(self.dataset_dir)
        self.finding_by_id = {f["_id"]: f for f in self.findings}
        self._index_findings()
        logger.info("insight loaded in %.2fs: %d structures, %d findings", time.time() - t0,
                    len(self.structure_files), len(self.findings))

    # -- structures ----------------------------------------------------- #

    def _proteins_for_stem(self, stems: list[str]) -> dict[str, list[str]]:
        """Canonical structure file stems -> protein IDs (one scan for all files)."""
        if not stems:
            return {}
        with self.lock:
            rows = self.store.execute(
                """SELECT regexp_replace(protein_id, '[^A-Za-z0-9_.-]', '_', 'g') AS stem, protein_id FROM proteins
                   WHERE regexp_replace(protein_id, '[^A-Za-z0-9_.-]', '_', 'g') IN (SELECT UNNEST(?::VARCHAR[]))""",
                [stems])
        out: dict[str, list[str]] = defaultdict(list)
        for stem, pid in rows:
            out[stem].append(pid)
        return out

    def structures_for(self, protein_id: str) -> list[dict[str, Any]]:
        return [r.to_dict() for r in self.structures.get(protein_id, [])]

    # -- similarity ----------------------------------------------------- #

    def similarity_state(self) -> str:
        if self._similar_state is None:
            try:
                from sharur.storage.vector_store import inspect_vector_index  # noqa: PLC0415

                self._similar_state = (inspect_vector_index(self.embeddings_h5).state
                                       if self.embeddings_h5.is_file() else "no_embeddings")
            except Exception:  # noqa: BLE001 - optional dependency
                self._similar_state = "unavailable"
        return self._similar_state

    def similar(self, protein_id: str, k: int = 12) -> dict[str, Any]:
        state = self.similarity_state()
        if state != "available":
            return {"state": state, "neighbors": []}
        with self._similar_lock:
            if self._similar is None:
                from sharur.storage.vector_store import FAISSStore  # noqa: PLC0415

                self._similar = FAISSStore(self.embeddings_h5, build_if_missing=False)
        t0 = time.perf_counter()
        pairs = self._similar.query(protein_id, k=k, include_distances=True)
        elapsed = time.perf_counter() - t0
        if not pairs:
            return {"state": "available", "neighbors": [], "ms": round(elapsed * 1000, 1)}
        ids = [p for p, _ in pairs]
        from sharur.architecture import _hits, compact, resolve  # noqa: PLC0415

        with self.lock:
            meta = {pid: (b, n) for pid, b, n in self.store.execute(
                "SELECT protein_id, bin_id, sequence_length FROM proteins WHERE protein_id IN "
                "(SELECT UNNEST(?::VARCHAR[]))", [ids])}
            best = dict(self.store.execute(
                """SELECT protein_id, ARG_MIN(COALESCE(NULLIF(name, ''), accession), COALESCE(evalue, 1))
                   FROM annotations WHERE protein_id IN (SELECT UNNEST(?::VARCHAR[])) GROUP BY 1""", [ids]))
            domains = _hits(self.store, ["pfam"], "AND a.protein_id IN (SELECT UNNEST(?::VARCHAR[]))", [ids])
        neighbors = []
        for pid, score in pairs:
            bin_id, length = meta.get(pid, (None, None))
            genome = self.catalog.by_bin.get(bin_id)
            neighbors.append({
                "protein_id": pid, "similarity": round(score, 4), "bin_id": bin_id,
                "lineage": genome.label if genome else "", "length": length,
                "architecture": compact([d.name for d in resolve(domains.get(pid, []))]),
                "best_hit": self.ctx.describe_hit(best.get(pid)) or "",
                "url": self.ctx.url("protein", pid), "genome_url": self.ctx.url("genome", bin_id) if bin_id else ""})
        return {"state": "available", "neighbors": neighbors, "ms": round(elapsed * 1000, 1)}

    # -- findings ------------------------------------------------------- #

    def _index_findings(self) -> None:
        self.by_protein: dict[str, list[str]] = defaultdict(list)
        self.by_genome: dict[str, list[str]] = defaultdict(list)
        self.by_system: dict[str, list[str]] = defaultdict(list)
        candidates: dict[str, set[str]] = defaultdict(set)
        for f in self.findings:
            fields = [f.get("protein_ids"), f.get("genes"), f.get("contigs"), f.get("genomes"),
                      f.get("bin_ids"), f.get("evidence")]
            for s in _strings(fields):
                candidates[s].add(f["_id"])
        names = list(candidates)
        proteins: dict[str, str] = {}
        contigs: dict[str, str] = {}
        if names:
            with self.lock:
                for start in range(0, len(names), 5000):
                    chunk = names[start:start + 5000]
                    proteins.update(dict(self.store.execute(
                        "SELECT protein_id, bin_id FROM proteins WHERE protein_id IN (SELECT UNNEST(?::VARCHAR[]))",
                        [chunk])))
                    contigs.update(dict(self.store.execute(
                        "SELECT contig_id, bin_id FROM contigs WHERE contig_id IN (SELECT UNNEST(?::VARCHAR[]))",
                        [chunk])))
        for name, ids in candidates.items():
            for fid in ids:
                if name in proteins:
                    self.by_protein[name].append(fid)
                    self.by_genome[proteins[name]].append(fid)
                if name in contigs:
                    self.by_genome[contigs[name]].append(fid)
                if name in self.catalog.by_bin:
                    self.by_genome[name].append(fid)
        types = {s["type"] for s in self.catalog.systems if s.get("type") and len(str(s["type"])) >= 4}
        member_type = {p: s["type"] for s in self.catalog.systems for p in s["proteins"]}
        for f in self.findings:
            text = f"{f.get('title', '')} {f.get('description', '')}"
            for t in types:
                if re.search(rf"(?<![\w-]){re.escape(str(t))}(?![\w-])", text):
                    self.by_system[str(t)].append(f["_id"])
        for pid, fids in self.by_protein.items():
            if pid in member_type:
                self.by_system[member_type[pid]].extend(fids)
        for table in (self.by_protein, self.by_genome, self.by_system):
            for k in list(table):
                table[k] = list(dict.fromkeys(table[k]))

    def _summaries(self, ids: list[str], limit: int = 12) -> list[dict[str, Any]]:
        return [self.summary(self.finding_by_id[i]) for i in ids[:limit] if i in self.finding_by_id]

    def summary(self, f: dict[str, Any]) -> dict[str, Any]:
        verification = f.get("verification") if isinstance(f.get("verification"), list) else []
        return {"id": f["_id"], "title": f.get("title") or f["_id"], "category": f.get("category") or "",
                "phase": f.get("phase") or "", "n_genomes": f.get("n_genomes"),
                "verifications": len(verification), "source": f["_source"],
                "url": self.ctx.url("finding", f["_id"])}

    def findings_for_protein(self, protein_id: str) -> list[dict[str, Any]]:
        return self._summaries(self.by_protein.get(protein_id, []))

    def findings_for_genome(self, bin_id: str) -> list[dict[str, Any]]:
        return self._summaries(self.by_genome.get(bin_id, []))

    def findings_for_system(self, system_type: str) -> list[dict[str, Any]]:
        return self._summaries(self.by_system.get(system_type, []))

    def findings_for_clade(self, rank: str | None, name: str | None) -> list[dict[str, Any]]:
        if not rank:
            return []
        members = {g.bin_id for g in self.catalog.clade(rank, name)}
        scores = defaultdict(int)
        for bin_id in members:
            for fid in self.by_genome.get(bin_id, []):
                scores[fid] += 1
        return self._summaries(sorted(scores, key=lambda fid: -scores[fid]), limit=10)

    def verify(self, finding_id: str) -> list[dict[str, Any]]:
        f = self.finding_by_id.get(finding_id)
        if f is None:
            raise KeyError(finding_id)
        records = f.get("verification") if isinstance(f.get("verification"), list) else []
        results = []
        for record in records:
            if not isinstance(record, dict):
                continue
            query = str(record.get("query") or "")
            entry = {"claim": record.get("claim", ""), "expected": _clip(record.get("expected"))}
            ok, why = is_safe_select(query)
            if not ok:
                entry.update(status="skipped", detail=why)
                results.append(entry)
                continue
            cursor = self.store.conn.cursor()
            timer = threading.Timer(VERIFY_TIMEOUT_S, cursor.interrupt)
            t0 = time.perf_counter()
            try:
                timer.start()
                rows = cursor.execute(query).fetchmany(200)
                status, observed = _compare(record.get("expected"), rows)
                entry.update(status=status, observed=_clip(observed))
            except Exception as exc:  # noqa: BLE001 - surfaced to the page
                timed_out = time.perf_counter() - t0 >= VERIFY_TIMEOUT_S - 0.05
                entry.update(status="error", detail=f"timed out after {VERIFY_TIMEOUT_S:.0f}s" if timed_out
                             else f"{type(exc).__name__}: {str(exc)[:200]}")
            finally:
                timer.cancel()
                cursor.close()
            entry["ms"] = round((time.perf_counter() - t0) * 1000, 1)
            results.append(entry)
        return results


# --------------------------------------------------------------------------- #
# Routes
# --------------------------------------------------------------------------- #


def register(app: FastAPI, ctx) -> Insight:
    insight = Insight(ctx)
    ctx.templates.env.globals["insight"] = insight
    app.state.insight = insight

    @app.get("/api/similar/{protein_id:path}")
    def similar(protein_id: str, k: int = Query(12, ge=1, le=50)):
        try:
            return JSONResponse(insight.similar(protein_id, k=k))
        except Exception as exc:  # noqa: BLE001 - panel shows the reason
            logger.exception("similarity query failed")
            return JSONResponse({"state": "error", "detail": type(exc).__name__, "neighbors": []})

    @app.get("/structure-file/{key}")
    def structure_file(key: int):
        record = insight.structure_files.get(key)
        if record is None or not record.path.is_file():
            raise HTTPException(404, "Structure not found")
        return PlainTextResponse(record.path.read_text(), media_type="chemical/x-pdb")

    @app.get("/findings", response_class=HTMLResponse)
    def findings_page(request: Request, q: str = Query("", max_length=200), category: str = Query("")):
        rows = [insight.summary(f) for f in insight.findings]
        needle = q.strip().lower()
        if needle:
            rows = [r for r in rows if needle in r["title"].lower() or needle in r["id"].lower()]
        if category:
            rows = [r for r in rows if r["category"] == category]
        categories = sorted({insight.summary(f)["category"] for f in insight.findings if f.get("category")})
        return ctx.render(request, "findings.html", "findings", rows=rows, q=q, category=category,
                          categories=categories, total=len(insight.findings))

    @app.get("/finding/{finding_id:path}", response_class=HTMLResponse)
    def finding_page(request: Request, finding_id: str):
        f = insight.finding_by_id.get(finding_id)
        if f is None:
            raise HTTPException(404, "Finding not found")
        proteins = [p for p, ids in insight.by_protein.items() if finding_id in ids]
        genomes = [b for b, ids in insight.by_genome.items() if finding_id in ids]
        systems = [t for t, ids in insight.by_system.items() if finding_id in ids]
        evidence = json.dumps(f.get("evidence"), indent=2, default=str) if f.get("evidence") else ""
        verification = [r for r in (f.get("verification") or []) if isinstance(r, dict)]
        return ctx.render(request, "finding.html", "findings", f=f, s=insight.summary(f), proteins=proteins[:60],
                          genomes=genomes[:60], systems=systems, evidence=evidence[:20000],
                          verification=verification)

    @app.post("/api/finding/{finding_id:path}/verify")
    def verify(finding_id: str):
        try:
            return JSONResponse({"results": insight.verify(finding_id)})
        except KeyError as exc:
            raise HTTPException(404, "Finding not found") from exc

    return insight
