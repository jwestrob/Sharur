"""Where genes sit on their contigs: distance to contig ends and truncation.

Metagenome-assembled genomes are fragmented, so a gene near a contig end may
have neighbors (or its own missing half) on another contig or outside the bin.
This module reports, per protein:

- ``genes_to_start`` / ``genes_to_end``: genes between this one and each contig
  end, ranked by coordinate (robust to gaps in ``gene_index``);
- ``bp_to_start`` / ``bp_to_end``: base pairs to each end; ``bp_to_end`` needs
  the assembly contig length (``contigs.length_source = 'assembly'``);
- ``truncated_start`` / ``truncated_end``: Prodigal ``partial=`` flags, i.e.
  the gene runs off that contig end;
- ``edge_status``: ``truncated``, ``contig_edge``, ``interior``, ``circular``
  or ``no_coordinates``.

It also reads true contig lengths from assemblies and Prodigal flags from
protein FASTA headers, for ingest and for backfilling existing databases.
"""

from __future__ import annotations

import gzip
import json
import re
from collections.abc import Iterable, Iterator
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

EDGE_GENES = 3
EDGE_BP = 5000
_PARTIAL = re.compile(r"(?:^|;)partial=([01]{2})(?:;|$)")


@dataclass(frozen=True)
class EdgeContext:
    protein_id: str
    contig_id: str
    contig_genes: int
    genes_to_start: int | None
    genes_to_end: int | None
    bp_to_start: int | None
    bp_to_end: int | None
    truncated_start: bool
    truncated_end: bool
    contig_length: int | None
    length_source: str | None
    edge_status: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @property
    def near_edge(self) -> bool:
        return self.edge_status in ("truncated", "contig_edge")


# --------------------------------------------------------------------------- #
# Inputs: assemblies and Prodigal headers
# --------------------------------------------------------------------------- #


def _open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path)


def fasta_lengths(path: Path) -> dict[str, int]:
    """Record ID (first header word) -> sequence length, streamed."""
    lengths: dict[str, int] = {}
    current, length = None, 0
    with _open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if current is not None:
                    lengths[current] = length
                fields = line[1:].split(None, 1)
                current, length = (fields[0] if fields else ""), 0
            else:
                length += len(line.strip())
    if current is not None:
        lengths[current] = length
    return lengths


def parse_partial(header_fields: Iterable[str] | None) -> str | None:
    """Prodigal ``partial=XY`` from the `` # ``-separated header fields."""
    for field in header_fields or ():
        match = _PARTIAL.search(field.strip())
        if match:
            return match.group(1)
    return None


def assemblies_from_stage00(stage00_dir: Path) -> dict[str, Path]:
    """genome_id -> prepared assembly FASTA, from the stage 00 manifest."""
    manifest = stage00_dir / "processing_manifest.json"
    if not manifest.is_file():
        return {}
    data = json.loads(manifest.read_text())
    paths = {}
    for entry in data.get("genomes", []):
        path = entry.get("output_path") or (stage00_dir / "genomes" / entry.get("filename", ""))
        if entry.get("genome_id") and Path(path).is_file():
            paths[entry["genome_id"]] = Path(path)
    return paths


def assemblies_from_dir(directory: Path) -> dict[str, Path]:
    """genome_id -> assembly FASTA for ``*.fna|fa|fasta[.gz]`` files in a directory."""
    suffixes = (".fna", ".fa", ".fasta", ".fas")
    paths = {}
    for path in sorted(directory.iterdir()):
        name = path.name[:-3] if path.name.endswith(".gz") else path.name
        for suffix in suffixes:
            if name.endswith(suffix):
                paths[name[: -len(suffix)]] = path
    return paths


def prodigal_partials(faa: Path) -> Iterator[tuple[str, str]]:
    """(protein_id, partial flag) for each Prodigal protein in a FASTA."""
    with _open_text(faa) as handle:
        for line in handle:
            if line.startswith(">"):
                parts = line[1:].rstrip("\n").split(" # ")
                flag = parse_partial(parts[4:]) if len(parts) > 4 else None
                if flag is not None:
                    yield parts[0].split()[0], flag


# --------------------------------------------------------------------------- #
# Queries
# --------------------------------------------------------------------------- #


def _columns(store, table: str) -> set[str]:
    return {r[0] for r in store.execute(
        "SELECT column_name FROM information_schema.columns WHERE table_name = ?", [table])}


def _status(rank: int, n: int, start: int, bp_end: int | None, truncated: tuple[bool, bool],
            circular: bool, has_coordinates: bool, edge_genes: int, edge_bp: int) -> str:
    if not has_coordinates:
        return "no_coordinates"
    if circular:
        return "circular"
    if any(truncated):
        return "truncated"
    near_start = rank < edge_genes and start - 1 < edge_bp
    near_end = (n - 1 - rank) < edge_genes and (bp_end is None or bp_end < edge_bp)
    return "contig_edge" if near_start or near_end else "interior"


def edge_context(store, protein_ids: Iterable[str], *, edge_genes: int = EDGE_GENES,
                 edge_bp: int = EDGE_BP) -> dict[str, EdgeContext]:
    """Contig-edge context for each protein (missing IDs are left out)."""
    ids = sorted(set(protein_ids))
    if not ids:
        return {}
    has_partial = "partial" in _columns(store, "proteins")
    has_source = "length_source" in _columns(store, "contigs")
    partial = "p.partial" if has_partial else "NULL"
    source = "c.length_source" if has_source else "NULL"
    rows = store.execute(
        f"""
        WITH _edge_targets AS (SELECT UNNEST(?::VARCHAR[]) AS protein_id),
        target_contigs AS (
            SELECT DISTINCT p.contig_id FROM proteins p JOIN _edge_targets t USING (protein_id)
        ),
        ranked AS (
            SELECT p.protein_id, p.contig_id, p.start, p.end_coord, {partial} AS partial,
                   ROW_NUMBER() OVER w - 1 AS rank, COUNT(*) OVER (PARTITION BY p.contig_id) AS n
            FROM proteins p JOIN target_contigs USING (contig_id)
            WINDOW w AS (PARTITION BY p.contig_id ORDER BY p.start, p.end_coord, p.protein_id)
        )
        SELECT r.protein_id, r.contig_id, r.rank, r.n, r.start, r.end_coord, r.partial,
               c.length, {source} AS length_source, COALESCE(c.is_circular, FALSE)
        FROM ranked r JOIN _edge_targets t USING (protein_id)
        LEFT JOIN contigs c ON c.contig_id = r.contig_id
        """, [ids])
    result = {}
    for pid, contig, rank, n, start, end, flag, length, length_source, circular in rows:
        has_coordinates = bool(start and start > 0)
        truncated = (bool(flag) and flag[0] == "1", bool(flag) and flag[1] == "1")
        bp_end = length - end if (length_source == "assembly" and length is not None and has_coordinates) else None
        result[pid] = EdgeContext(
            protein_id=pid, contig_id=contig, contig_genes=int(n),
            genes_to_start=int(rank) if has_coordinates else None,
            genes_to_end=int(n - 1 - rank) if has_coordinates else None,
            bp_to_start=start - 1 if has_coordinates else None,
            bp_to_end=bp_end, truncated_start=truncated[0], truncated_end=truncated[1],
            contig_length=length if length_source == "assembly" else None,
            length_source=length_source,
            edge_status=_status(int(rank), int(n), start or 0, bp_end, truncated, bool(circular),
                                has_coordinates, edge_genes, edge_bp),
        )
    return result


def describe_edge(context: EdgeContext) -> str:
    """One-line human summary."""
    if context.edge_status == "no_coordinates":
        return "no contig coordinates"
    if context.edge_status == "circular":
        return "circular contig"
    parts = []
    if context.truncated_start or context.truncated_end:
        sides = [s for s, t in (("start", context.truncated_start), ("end", context.truncated_end)) if t]
        parts.append(f"runs off contig {' and '.join(sides)}")
    end_bp = f", {context.bp_to_end:,} bp" if context.bp_to_end is not None else ""
    parts.append(f"{context.genes_to_start} genes ({context.bp_to_start:,} bp) from contig start, "
                 f"{context.genes_to_end} genes{end_bp} from contig end")
    if context.length_source != "assembly":
        parts.append("contig length from gene span")
    return f"{context.edge_status}: " + "; ".join(parts)


# --------------------------------------------------------------------------- #
# Backfill
# --------------------------------------------------------------------------- #


def backfill_contig_context(conn, *, assemblies: dict[str, Path] | None = None,
                            protein_faas: Iterable[Path] = (), threads: int | None = None) -> dict[str, int]:
    """Write assembly contig lengths and Prodigal partial flags into a database.

    ``assemblies`` maps bin_id -> assembly FASTA; only contigs already in that
    bin are updated. A contig whose assembly record is shorter than its gene
    span means the wrong assembly was supplied, and nothing is written.
    """
    import os  # noqa: PLC0415
    from concurrent.futures import ThreadPoolExecutor  # noqa: PLC0415

    import pandas as pd  # noqa: PLC0415

    from sharur.storage.migrations import run_migrations  # noqa: PLC0415

    run_migrations(conn)
    workers = threads or os.cpu_count() or 1
    stats = {"contigs_updated": 0, "contigs_unmatched": 0, "proteins_flagged": 0}
    with ThreadPoolExecutor(max_workers=workers) as executor:
        bins = sorted(assemblies or {})
        lengths = list(executor.map(lambda b: fasta_lengths(assemblies[b]), bins))
        flags = list(executor.map(lambda f: list(prodigal_partials(f)), list(protein_faas)))
    conn.execute("BEGIN TRANSACTION")
    try:
        if bins:
            asm = pd.DataFrame(
                [(b, cid, n) for b, table in zip(bins, lengths, strict=True) for cid, n in table.items()],
                columns=["bin_id", "contig_id", "length"])
            conn.register("_asm", asm)
            bad = conn.execute(
                """SELECT c.contig_id, a.length, MAX(p.end_coord) FROM contigs c
                   JOIN _asm a ON a.contig_id = c.contig_id AND a.bin_id = c.bin_id
                   JOIN proteins p ON p.contig_id = c.contig_id
                   GROUP BY ALL HAVING MAX(p.end_coord) > a.length LIMIT 5""").fetchall()
            if bad:
                raise ValueError(f"Assembly contigs shorter than their genes (wrong assembly?): {bad}")
            stats["contigs_updated"] = conn.execute(
                """UPDATE contigs SET length = a.length, length_source = 'assembly'
                   FROM _asm a WHERE a.contig_id = contigs.contig_id AND a.bin_id = contigs.bin_id"""
            ).fetchone()[0]
            stats["contigs_unmatched"] = conn.execute(
                """SELECT COUNT(*) FROM contigs WHERE bin_id IN (SELECT DISTINCT bin_id FROM _asm)
                   AND COALESCE(length_source, '') <> 'assembly'""").fetchone()[0]
            conn.unregister("_asm")
        conn.execute("UPDATE contigs SET length_source = 'gene_span' WHERE length_source IS NULL")
        pairs = [pair for chunk in flags for pair in chunk]
        if pairs:
            conn.register("_partial", pd.DataFrame(pairs, columns=["protein_id", "partial"]))
            stats["proteins_flagged"] = conn.execute(
                """UPDATE proteins SET partial = f.partial FROM _partial f
                   WHERE f.protein_id = proteins.protein_id""").fetchone()[0]
            conn.unregister("_partial")
        conn.execute("COMMIT")
    except Exception:
        conn.execute("ROLLBACK")
        raise
    return stats
