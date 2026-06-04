"""Cross-dataset finding mirror.

The per-dataset OpsStore (`data/<dataset>/sharur_ops.db`) is siloed by design:
each analysis run has its own coordination history. The cost is that
**patterns across datasets do not surface**. If `DUF6088` shows up adjacent
to AbiEii toxin in CoronaMine, *and* the same DUF surfaces in srvp_bacteria_pb,
no agent in either dataset will notice — the local OpsStores never talk.

`GlobalOpsStore` is a per-user SQLite database at `~/.sharur/global_ops.db`
that mirrors high-novelty findings from every local OpsStore and exposes
cross-dataset queries:

  - `find_by_accession(acc)` — every finding citing PF13437 across datasets
  - `find_by_category(cat)` — defense-system findings across all datasets
  - `co_occurring_accessions(acc)` — accessions that show up alongside `acc`
    in findings from multiple datasets

The mirror is lossy on purpose: we copy summary, novelty, accession(s), and
a link back to the source OpsStore. The full evidence stays in the local
store. The global view answers "where else have we seen this?" — not "what
exactly happened that day."

Usage from a coordinator:

    from sharur.ops.global_store import GlobalOpsStore, mirror_from_local
    g = GlobalOpsStore()
    mirror_from_local(g, "data/coronamine_boiler_100nm/sharur_ops.db", min_novelty=2)
    print(g.find_by_accession("DUF6088"))

Or from the CLI:

    python scripts/global_ops_sync.py --dataset data/coronamine_boiler_100nm
"""

from __future__ import annotations

import json
import os
import re
import sqlite3
import time
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Iterator, Optional

DEFAULT_GLOBAL_DB = Path(
    os.environ.get("SHARUR_GLOBAL_OPS", str(Path.home() / ".sharur" / "global_ops.db"))
)

# Accession patterns we recognize when scanning evidence/summary JSON blobs.
# Order matters: PFAM short form before generic uppercase, KO before generic.
ACCESSION_PATTERNS = (
    re.compile(r"\bPF\d{5}\b"),                # PFAM
    re.compile(r"\bK\d{5}\b"),                  # KEGG KO
    re.compile(r"\bDUF\d{2,5}\b"),              # DUF
    re.compile(r"\bCOG\d{4}\b"),                # COG
    re.compile(r"\bGH\d{1,3}[a-z]?\b"),         # CAZy GH families
    re.compile(r"\bGT\d{1,3}[a-z]?\b"),         # CAZy GT
    re.compile(r"\bCBM\d{1,3}\b"),              # CAZy CBM
    re.compile(r"\bVOG\d+\b"),                  # VOGdb
)


# ---------------------------------------------------------------------------
# Schema
# ---------------------------------------------------------------------------

SCHEMA = """
CREATE TABLE IF NOT EXISTS datasets (
    dataset_path TEXT PRIMARY KEY,
    dataset_name TEXT NOT NULL,
    last_synced_ts REAL,
    n_findings INTEGER DEFAULT 0
);

CREATE TABLE IF NOT EXISTS findings_mirror (
    -- composite key: source dataset + local finding id
    source_dataset TEXT NOT NULL,
    finding_id TEXT NOT NULL,
    summary TEXT NOT NULL,
    domain TEXT,
    finding_type TEXT,
    novelty INTEGER NOT NULL,
    confidence REAL,
    agent_id TEXT,
    ts REAL,
    evidence_json TEXT,
    PRIMARY KEY (source_dataset, finding_id),
    FOREIGN KEY (source_dataset) REFERENCES datasets(dataset_path)
);
CREATE INDEX IF NOT EXISTS idx_fm_novelty ON findings_mirror(novelty);
CREATE INDEX IF NOT EXISTS idx_fm_domain ON findings_mirror(domain);
CREATE INDEX IF NOT EXISTS idx_fm_ts ON findings_mirror(ts);

CREATE TABLE IF NOT EXISTS accession_mentions (
    accession TEXT NOT NULL,
    source_dataset TEXT NOT NULL,
    finding_id TEXT NOT NULL,
    novelty INTEGER NOT NULL,
    PRIMARY KEY (accession, source_dataset, finding_id),
    FOREIGN KEY (source_dataset, finding_id)
        REFERENCES findings_mirror(source_dataset, finding_id) ON DELETE CASCADE
);
CREATE INDEX IF NOT EXISTS idx_am_accession ON accession_mentions(accession);
CREATE INDEX IF NOT EXISTS idx_am_dataset ON accession_mentions(source_dataset);

CREATE TABLE IF NOT EXISTS sync_log (
    ts REAL NOT NULL,
    source_dataset TEXT NOT NULL,
    n_new INTEGER NOT NULL,
    note TEXT
);
"""


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class MirroredFinding:
    source_dataset: str
    finding_id: str
    summary: str
    domain: Optional[str]
    finding_type: Optional[str]
    novelty: int
    confidence: Optional[float]
    agent_id: Optional[str]
    ts: Optional[float]
    accessions: list[str]


# ---------------------------------------------------------------------------
# Accession extraction
# ---------------------------------------------------------------------------

def extract_accessions(*texts: str) -> set[str]:
    """Pull PF/K/DUF/COG/CAZy/VOG accessions from any number of text blobs."""
    found: set[str] = set()
    for t in texts:
        if not t:
            continue
        for pat in ACCESSION_PATTERNS:
            for m in pat.findall(t):
                found.add(m)
    return found


# ---------------------------------------------------------------------------
# The store
# ---------------------------------------------------------------------------

class GlobalOpsStore:
    """Per-user mirror of high-novelty findings across all Sharur datasets."""

    def __init__(self, path: Path | str | None = None):
        self.path = Path(path) if path else DEFAULT_GLOBAL_DB
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._init_schema()

    def _init_schema(self) -> None:
        with self._conn() as conn:
            conn.executescript(SCHEMA)

    @contextmanager
    def _conn(self) -> Iterator[sqlite3.Connection]:
        conn = sqlite3.connect(self.path)
        conn.row_factory = sqlite3.Row
        try:
            yield conn
            conn.commit()
        finally:
            conn.close()

    # -- writes ----------------------------------------------------------

    def register_dataset(self, dataset_path: Path | str) -> None:
        p = Path(dataset_path).resolve()
        with self._conn() as conn:
            conn.execute(
                "INSERT OR IGNORE INTO datasets (dataset_path, dataset_name) VALUES (?, ?)",
                (str(p), p.name),
            )

    def upsert_finding(self, f: MirroredFinding) -> None:
        with self._conn() as conn:
            conn.execute(
                """
                INSERT INTO findings_mirror
                    (source_dataset, finding_id, summary, domain, finding_type,
                     novelty, confidence, agent_id, ts, evidence_json)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                ON CONFLICT (source_dataset, finding_id) DO UPDATE SET
                    summary = excluded.summary,
                    domain = excluded.domain,
                    finding_type = excluded.finding_type,
                    novelty = excluded.novelty,
                    confidence = excluded.confidence,
                    agent_id = excluded.agent_id,
                    ts = excluded.ts
                """,
                (f.source_dataset, f.finding_id, f.summary, f.domain,
                 f.finding_type, f.novelty, f.confidence, f.agent_id, f.ts,
                 json.dumps([])),
            )
            # Refresh accession mentions
            conn.execute(
                "DELETE FROM accession_mentions WHERE source_dataset = ? AND finding_id = ?",
                (f.source_dataset, f.finding_id),
            )
            for acc in f.accessions:
                conn.execute(
                    """
                    INSERT OR IGNORE INTO accession_mentions
                        (accession, source_dataset, finding_id, novelty)
                    VALUES (?, ?, ?, ?)
                    """,
                    (acc, f.source_dataset, f.finding_id, f.novelty),
                )

    def log_sync(self, source_dataset: str, n_new: int, note: str = "") -> None:
        with self._conn() as conn:
            conn.execute(
                "INSERT INTO sync_log (ts, source_dataset, n_new, note) VALUES (?, ?, ?, ?)",
                (time.time(), source_dataset, n_new, note),
            )
            conn.execute(
                """UPDATE datasets
                   SET last_synced_ts = ?, n_findings = (
                       SELECT COUNT(*) FROM findings_mirror WHERE source_dataset = ?
                   )
                   WHERE dataset_path = ?""",
                (time.time(), source_dataset, source_dataset),
            )

    # -- queries ---------------------------------------------------------

    def find_by_accession(self, accession: str) -> list[dict]:
        """Every finding that references this accession across datasets."""
        with self._conn() as conn:
            rows = conn.execute(
                """
                SELECT fm.source_dataset, fm.finding_id, fm.summary, fm.domain,
                       fm.novelty, fm.agent_id, am.accession
                FROM accession_mentions am
                JOIN findings_mirror fm
                  ON fm.source_dataset = am.source_dataset
                 AND fm.finding_id = am.finding_id
                WHERE am.accession = ?
                ORDER BY fm.novelty DESC, fm.ts DESC
                """,
                (accession,),
            ).fetchall()
        return [dict(r) for r in rows]

    def find_by_domain(self, domain: str, min_novelty: int = 0) -> list[dict]:
        with self._conn() as conn:
            rows = conn.execute(
                """
                SELECT source_dataset, finding_id, summary, domain, novelty, agent_id, ts
                FROM findings_mirror
                WHERE domain = ? AND novelty >= ?
                ORDER BY novelty DESC, ts DESC
                """,
                (domain, min_novelty),
            ).fetchall()
        return [dict(r) for r in rows]

    def co_occurring_accessions(
        self, accession: str, min_datasets: int = 2,
    ) -> list[dict]:
        """Accessions that co-occur with `accession` in findings across N+ datasets.

        Useful for spotting recurring guilt-by-association patterns:
        "every time we see DUF6088, we also see AbiEii — and it's in 3 datasets."
        """
        with self._conn() as conn:
            rows = conn.execute(
                """
                WITH target_findings AS (
                    SELECT DISTINCT source_dataset, finding_id
                    FROM accession_mentions
                    WHERE accession = ?
                )
                SELECT am.accession,
                       COUNT(DISTINCT am.source_dataset) AS n_datasets,
                       COUNT(*) AS n_findings
                FROM accession_mentions am
                JOIN target_findings tf
                  ON tf.source_dataset = am.source_dataset
                 AND tf.finding_id = am.finding_id
                WHERE am.accession != ?
                GROUP BY am.accession
                HAVING n_datasets >= ?
                ORDER BY n_datasets DESC, n_findings DESC
                """,
                (accession, accession, min_datasets),
            ).fetchall()
        return [dict(r) for r in rows]

    def recent(self, min_novelty: int = 2, limit: int = 50) -> list[dict]:
        with self._conn() as conn:
            rows = conn.execute(
                """
                SELECT source_dataset, finding_id, summary, domain, novelty, agent_id, ts
                FROM findings_mirror
                WHERE novelty >= ?
                ORDER BY ts DESC NULLS LAST
                LIMIT ?
                """,
                (min_novelty, limit),
            ).fetchall()
        return [dict(r) for r in rows]

    def datasets(self) -> list[dict]:
        with self._conn() as conn:
            rows = conn.execute(
                "SELECT dataset_path, dataset_name, last_synced_ts, n_findings FROM datasets ORDER BY dataset_name"
            ).fetchall()
        return [dict(r) for r in rows]


# ---------------------------------------------------------------------------
# Sync helpers
# ---------------------------------------------------------------------------

def _local_findings(local_ops_path: Path, min_novelty: int) -> Iterable[sqlite3.Row]:
    conn = sqlite3.connect(local_ops_path)
    conn.row_factory = sqlite3.Row
    try:
        yield from conn.execute(
            """
            SELECT id, agent_id, ts, finding_type, domain, summary,
                   evidence, confidence, novelty
            FROM findings
            WHERE novelty >= ?
            """,
            (min_novelty,),
        ).fetchall()
    finally:
        conn.close()


def mirror_from_local(
    g: GlobalOpsStore,
    local_ops_path: Path | str,
    *,
    min_novelty: int = 2,
    dataset_root: Optional[Path | str] = None,
) -> int:
    """Pull high-novelty findings from a per-dataset OpsStore into the global mirror.

    Args:
        g: Global store to write into.
        local_ops_path: Path to `data/<name>/sharur_ops.db`.
        min_novelty: Minimum novelty score to mirror. Default 2.
        dataset_root: Override what to register as the dataset path. Defaults
            to the parent of `local_ops_path`.

    Returns the number of findings mirrored (new or updated).
    """
    local = Path(local_ops_path).resolve()
    if not local.exists():
        return 0
    ds = Path(dataset_root).resolve() if dataset_root else local.parent
    g.register_dataset(ds)

    n = 0
    for row in _local_findings(local, min_novelty):
        evidence_text = row["evidence"] or ""
        acc = extract_accessions(row["summary"] or "", evidence_text)
        g.upsert_finding(MirroredFinding(
            source_dataset=str(ds),
            finding_id=row["id"],
            summary=row["summary"],
            domain=row["domain"],
            finding_type=row["finding_type"],
            novelty=row["novelty"],
            confidence=row["confidence"],
            agent_id=row["agent_id"],
            ts=row["ts"],
            accessions=sorted(acc),
        ))
        n += 1
    g.log_sync(str(ds), n, note=f"min_novelty={min_novelty}")
    return n


def mirror_from_findings_jsonl(
    g: GlobalOpsStore,
    findings_path: Path | str,
    *,
    dataset_root: Optional[Path | str] = None,
    min_novelty: int = 2,
) -> int:
    """Backfill from a findings.jsonl when no per-dataset OpsStore exists.

    Useful for older datasets that wrote findings.jsonl but never ran the
    coordinator (no sharur_ops.db).
    """
    fp = Path(findings_path).resolve()
    if not fp.exists():
        return 0
    ds = Path(dataset_root).resolve() if dataset_root else fp.parent.parent
    g.register_dataset(ds)

    n = 0
    for line in fp.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            f = json.loads(line)
        except json.JSONDecodeError:
            continue
        novelty = int(f.get("novelty") or 0)
        if novelty < min_novelty:
            continue
        evidence_text = json.dumps(f.get("evidence") or {})
        title = f.get("title") or f.get("summary") or ""
        description = f.get("description") or ""
        acc = extract_accessions(title, description, evidence_text)
        g.upsert_finding(MirroredFinding(
            source_dataset=str(ds),
            finding_id=str(f.get("id") or f"backfill_{n}"),
            summary=title,
            domain=f.get("category"),
            finding_type=f.get("phase"),
            novelty=novelty,
            confidence=f.get("confidence"),
            agent_id=(f.get("provenance") or {}).get("agent_id"),
            ts=None,
            accessions=sorted(acc),
        ))
        n += 1
    g.log_sync(str(ds), n, note=f"backfill from findings.jsonl, min_novelty={min_novelty}")
    return n


__all__ = [
    "GlobalOpsStore",
    "MirroredFinding",
    "mirror_from_local",
    "mirror_from_findings_jsonl",
    "extract_accessions",
    "DEFAULT_GLOBAL_DB",
]
