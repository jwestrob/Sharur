"""Build compact semantic-term artifacts from a dataset's ``semantic_terms`` table.

Ported from the certified isolated prototype builders with the same SQL, row order,
dictionary order and encodings, so a rebuild of the same source reproduces the
certified payload bytes. The source opens read-only; its size, mtime, inode and
full SHA-256 are checked before and after each build and recorded in the manifest.
Every output directory must be fresh.
"""

from __future__ import annotations

import array
import contextlib
import importlib.metadata
import json
import resource
import struct
import sys
import tempfile
import time
from pathlib import Path

import duckdb
import numpy as np
from pyroaring import BitMap

from sharur.semantic_index.formats import (
    ACTIVE_EXPRESSION,
    FIELDS,
    FORWARD_FORMAT,
    MEMBERSHIP_FORMAT,
    SCOPE_FORMAT,
    sha256_file,
)


FORWARD_ORDER = ("term_kind", "term_id", "source_db", "source_accession", "facet", "relation")


class BuildError(RuntimeError):
    """A source fails a build precondition or changed during a build."""


def source_guard(path: str | Path) -> dict[str, int]:
    """Size, mtime and inode of a checkpointed source (a pending WAL is refused)."""
    p = Path(path)
    if Path(str(p) + ".wal").exists():
        raise BuildError(f"Source has a pending WAL; checkpoint it first: {p}")
    s = p.stat()
    return {"bytes": s.st_size, "mtime_ns": s.st_mtime_ns, "inode": s.st_ino}


def _versions() -> dict:
    import pyroaring  # noqa: PLC0415

    return {"python": sys.version,
            "packages": {n: importlib.metadata.version(n) for n in ("duckdb", "numpy", "pyroaring")},
            "croaring": pyroaring.__croaring_version__}


def _peak_rss_bytes() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform == "darwin" else value * 1024)


def _fresh(output: str | Path) -> Path:
    out = Path(output).resolve()
    if out.exists():
        raise FileExistsError(f"Build output exists; choose a fresh directory: {out}")
    out.mkdir(parents=True)
    return out


def _connect(source: Path, threads: int, memory_limit: str):
    return duckdb.connect(str(source), read_only=True, config={"threads": threads, "memory_limit": memory_limit})


@contextlib.contextmanager
def _session(source: Path, threads: int, memory_limit: str, parent: Path, prefix: str):
    """Read-only source connection spilling under ``parent``; closed before its spill directory is removed."""
    with tempfile.TemporaryDirectory(prefix=prefix, dir=parent) as spill:
        con = _connect(source, threads, memory_limit)
        try:
            con.execute("SET temp_directory = ?", [spill])
            yield con
        finally:
            con.close()


def _payloads(out: Path) -> dict[str, dict]:
    return {p.name: {"bytes": p.stat().st_size, "sha256": sha256_file(p)} for p in sorted(out.iterdir()) if p.is_file()}


def _finish_source(source: Path, before: dict, source_hash: str) -> dict:
    after = source_guard(source)
    if after != before or sha256_file(source) != source_hash:
        raise BuildError("Source changed during a read-only build")
    return after


def _write_dictionary(con, table: str, column: str, out: Path, prefix: str) -> int:
    con.execute(f"SELECT {column} FROM {table} ORDER BY id")
    count, position = 0, 0
    with (out / f"{prefix}.utf8").open("wb") as blob, (out / f"{prefix}.offsets.u64").open("wb") as offsets:
        offsets.write(struct.pack("<Q", 0))
        while rows := con.fetchmany(65536):
            for (text,) in rows:
                encoded = text.encode("utf-8")
                blob.write(encoded)
                position += len(encoded)
                offsets.write(struct.pack("<Q", position))
                count += 1
    return count


def build_membership(source: str | Path, output: str | Path, *, seal: str | Path | None = None,
                     threads: int = 2, memory_limit: str = "4GB") -> dict:
    """Active protein/term membership postings plus the shared protein and term dictionaries."""
    source = Path(source).resolve()
    before = source_guard(source)
    source_hash = sha256_file(source)
    out = _fresh(output)
    started = time.perf_counter()
    with _session(source, threads, memory_limit, out.parent, "spill-") as con:
        n_proteins = con.execute("SELECT count(*) FROM proteins").fetchone()[0]
        if n_proteins > 2**32:
            raise BuildError("Protein universe exceeds uint32")
        invalid = con.execute(f"""SELECT count(*) FROM semantic_terms t
            LEFT JOIN proteins p ON p.protein_id=t.protein_id WHERE {ACTIVE_EXPRESSION}
            AND (t.term_id IS NULL OR t.protein_id IS NULL OR p.protein_id IS NULL)""").fetchone()[0]
        if invalid:
            raise BuildError("Invalid or orphaned active membership keys")
        if con.execute("SELECT count(*)-count(DISTINCT protein_id) FROM proteins").fetchone()[0]:
            raise BuildError("Protein identifiers must be unique and non-NULL")
        con.execute("""CREATE TEMP TABLE protein_dictionary AS
            SELECT protein_id, (row_number() OVER (ORDER BY encode(protein_id))-1)::UINTEGER id
            FROM proteins""")
        con.execute(f"""CREATE TEMP TABLE term_dictionary AS
            SELECT term_id, (row_number() OVER (ORDER BY encode(term_id))-1)::UINTEGER id
            FROM (SELECT DISTINCT term_id FROM semantic_terms WHERE {ACTIVE_EXPRESSION})""")
        n_terms = con.execute("SELECT count(*) FROM term_dictionary").fetchone()[0]
        _write_dictionary(con, "protein_dictionary", "protein_id", out, "proteins")
        _write_dictionary(con, "term_dictionary", "term_id", out, "terms")
        membership_sql = f"""SELECT DISTINCT td.id tid, pd.id pid
            FROM semantic_terms st JOIN protein_dictionary pd USING(protein_id)
            JOIN term_dictionary td USING(term_id)
            WHERE {ACTIVE_EXPRESSION} ORDER BY tid,pid"""
        con.execute(membership_sql)
        universe = BitMap()
        counts = {"memberships": 0, "list_terms": 0, "roaring_terms": 0, "list_bytes": 0, "roaring_bytes": 0}
        with (out / "postings.bin").open("wb") as postings, (out / "catalog.bin").open("wb") as catalog:
            def flush(tid, values):
                if tid != counts["list_terms"] + counts["roaring_terms"]:
                    raise BuildError("Non-dense or out-of-order term stream")
                raw = np.asarray(values, dtype="<u4").tobytes()
                bm = BitMap(values)
                bm.run_optimize()
                packed = bm.serialize()
                encoding, payload = (1, packed) if len(packed) < len(raw) else (0, raw)
                catalog.write(struct.pack("<QQIB3x", postings.tell(), len(payload), len(values), encoding))
                postings.write(payload)
                universe.update(bm)
                counts["roaring_terms" if encoding else "list_terms"] += 1
                counts["roaring_bytes" if encoding else "list_bytes"] += len(payload)
                counts["memberships"] += len(values)

            term, ids = None, array.array("I")
            while rows := con.fetchmany(65536):
                for tid, pid in rows:
                    if tid != term:
                        if term is not None:
                            flush(term, ids)
                        term, ids = tid, array.array("I")
                    ids.append(pid)
            if term is not None:
                flush(term, ids)
        (out / "active_universe.roaring").write_bytes(universe.serialize())
        if counts["list_terms"] + counts["roaring_terms"] != n_terms:
            raise BuildError("Term catalog count mismatch")
    after = _finish_source(source, before, source_hash)
    files = _payloads(out)
    seal_info = None
    if seal:
        seal_info = {"sha256": sha256_file(seal), "dataset_id": json.loads(Path(seal).read_text()).get("dataset_id")}
    manifest = {"format": MEMBERSHIP_FORMAT, "endianness": "little", "id_order": "UTF-8 byte lexicographic",
                "active_expression": ACTIVE_EXPRESSION, "universe": "proteins with at least one active membership",
                "source": {"name": source.name, "sha256": source_hash, "before": before, "after": after,
                           "seal": seal_info},
                "versions": _versions(), "protein_count": n_proteins, "term_count": n_terms,
                "active_protein_count": len(universe), **counts, "files": files,
                "builder_sha256": sha256_file(__file__), "membership_sql": membership_sql,
                "wall_seconds": time.perf_counter() - started, "peak_rss_bytes": _peak_rss_bytes(),
                "payload_bytes": sum(v["bytes"] for v in files.values()),
                "deployment_payload_bytes": sum(v["bytes"] for v in files.values()),
                "optional_numeric_pairs_included": False}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def _code_type(column: str, count: int) -> str:
    if count > 2**32 - 1:
        raise BuildError("Dictionary exceeds uint32 code space")
    if column in ("term_id", "source_accession"):
        return "<u4"
    return "u1" if count <= 255 else "<u2" if count <= 65535 else "<u4"


def _write_field_dictionary(con, column: str, root: Path) -> dict:
    table = "dict_" + column
    con.execute(f"""CREATE TEMP TABLE {table} AS SELECT text_value,
        (row_number() OVER (ORDER BY encode(text_value)))::UINTEGER code
        FROM (SELECT DISTINCT {column} AS text_value FROM semantic_terms WHERE {column} IS NOT NULL)""")
    con.execute(f"SELECT text_value FROM {table} ORDER BY code")
    offsets, count = [0], 0
    with (root / f"{column}.utf8").open("wb") as blob:
        while rows := con.fetchmany(65536):
            for (text,) in rows:
                encoded = text.encode("utf-8")
                blob.write(encoded)
                offsets.append(offsets[-1] + len(encoded))
                count += 1
    np.asarray(offsets, dtype="<u8").tofile(root / f"{column}.offsets.u64")
    return {"non_null_values": count, "null_code": 0, "code_type": _code_type(column, count),
            "id_order": "nonzero codes in UTF-8 byte lexicographic order"}


def _matching_binding(source: Path, matching: Path) -> tuple[dict, dict, dict]:
    mm = json.loads((matching / "manifest.json").read_text())
    before = source_guard(source)
    if before != mm["source"]["before"]:
        raise BuildError("Source differs from the membership snapshot")
    if sha256_file(source) != mm["source"]["sha256"]:
        raise BuildError("Source checksum differs from the membership snapshot")
    shared = {name: mm["files"][name] for name in ("proteins.utf8", "proteins.offsets.u64")}
    for name, expected in shared.items():
        path = matching / name
        if path.stat().st_size != expected["bytes"] or sha256_file(path) != expected["sha256"]:
            raise BuildError("Shared protein dictionary checksum mismatch")
    return mm, before, shared


def build_forward(source: str | Path, matching: str | Path, output: str | Path, *,
                  threads: int = 2, memory_limit: str = "4GB") -> dict:
    """Every original ``semantic_terms`` row, CSR by the membership protein IDs."""
    source, matching = Path(source).resolve(), Path(matching).resolve()
    mm, before, shared = _matching_binding(source, matching)
    source_hash = mm["source"]["sha256"]
    out = _fresh(output)
    started = time.perf_counter()
    with _session(source, threads, memory_limit, out.parent, "forward-spill-") as con:
        n = con.execute("SELECT count(*) FROM proteins").fetchone()[0]
        if n != mm["protein_count"]:
            raise BuildError("Protein dictionary count changed")
        if con.execute("""SELECT count(*) FROM semantic_terms st LEFT JOIN proteins p USING(protein_id)
                          WHERE p.protein_id IS NULL""").fetchone()[0]:
            raise BuildError("Source contains orphan/NULL protein term rows")
        source_count = con.execute("SELECT count(*) FROM semantic_terms").fetchone()[0]
        con.execute("""CREATE TEMP TABLE protein_dictionary AS SELECT protein_id,
            (row_number() OVER (ORDER BY encode(protein_id))-1)::UINTEGER pid FROM proteins""")
        dictionaries = {field: _write_field_dictionary(con, field, out) for field in FIELDS}
        dtype = np.dtype([(f, dictionaries[f]["code_type"]) for f in FIELDS])
        selection = ", ".join(f"coalesce(d_{f}.code,0)::UINTEGER AS {f}" for f in FIELDS)
        joins = " ".join(f"LEFT JOIN dict_{f} d_{f} ON st.{f}=d_{f}.text_value" for f in FIELDS)
        order = ", ".join(f"encode(st.{f}) ASC NULLS LAST" for f in FORWARD_ORDER)
        sql = f"""SELECT pd.pid, {selection} FROM semantic_terms st
                JOIN protein_dictionary pd USING(protein_id) {joins}
                ORDER BY pd.pid, {order}"""
        con.execute(sql)
        counts = np.zeros(n + 1, dtype="<u8")
        total = 0
        with (out / "rows.bin").open("wb") as rows_file:
            while rows := con.fetchmany(65536):
                values = np.asarray(rows, dtype="<u4")
                pids, frequency = np.unique(values[:, 0], return_counts=True)
                counts[pids.astype(np.uint64) + 1] += frequency.astype(np.uint64)
                packed = np.empty(len(rows), dtype=dtype)
                for i, f in enumerate(FIELDS):
                    packed[f] = values[:, i + 1]
                rows_file.write(packed.tobytes())
                total += len(rows)
        np.cumsum(counts, out=counts)
        counts.tofile(out / "protein_offsets.u64")
        if total != source_count or int(counts[-1]) != total:
            raise BuildError("Rich row total differs from the source count")
    after = _finish_source(source, before, source_hash)
    files = _payloads(out)
    manifest = {"format": FORWARD_FORMAT, "scope": "all original seven-column semantic_terms, exact tuple multiplicity",
                "source": {"name": source.name, "sha256": source_hash, "before": before, "after": after},
                "matching_manifest_sha256": sha256_file(matching / "manifest.json"),
                "shared_protein_dictionary": {"protein_count": n, "files": shared, "codes": "matching uint32 pid"},
                "fields": list(FIELDS), "dictionaries": dictionaries, "row_dtype": dtype.descr,
                "row_bytes": dtype.itemsize, "rows": total, "offsets": "little-endian uint64, n+1 sentinel",
                "null_semantics": "code0=NULL; empty string is its own nonzero dictionary entry",
                "api_order": list(FORWARD_ORDER[:4]), "tie_only_order_extension": list(FORWARD_ORDER[4:]),
                "null_sort": "NULLS LAST", "sql": sql, "files": files,
                "payload_bytes": sum(f["bytes"] for f in files.values()), "versions": _versions(),
                "builder_sha256": sha256_file(__file__), "wall_seconds": time.perf_counter() - started,
                "peak_rss_bytes": _peak_rss_bytes()}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def build_genome_scope(source: str | Path, matching: str | Path, output: str | Path, *,
                       threads: int = 2, memory_limit: str = "4GB") -> dict:
    """Per-genome protein postings over the membership protein IDs; NULL owners kept separately."""
    source, matching = Path(source).resolve(), Path(matching).resolve()
    mm, before, shared = _matching_binding(source, matching)
    source_hash = mm["source"]["sha256"]
    out = _fresh(output)
    started = time.perf_counter()
    with _session(source, threads, memory_limit, out.parent, "scope-spill-") as con:
        n = con.execute("SELECT count(*) FROM proteins").fetchone()[0]
        if n != mm["protein_count"]:
            raise BuildError("Protein universe differs from the membership snapshot")
        if con.execute("SELECT count(*)-count(DISTINCT protein_id) FROM proteins").fetchone()[0]:
            raise BuildError("Protein IDs must be unique and non-NULL")
        con.execute("""CREATE TEMP TABLE protein_dictionary AS SELECT protein_id,
            (row_number() OVER (ORDER BY encode(protein_id))-1)::UINTEGER pid FROM proteins""")
        con.execute("""CREATE TEMP TABLE genome_dictionary AS SELECT bin_id,
            (row_number() OVER (ORDER BY encode(bin_id))-1)::UINTEGER gid
            FROM (SELECT bin_id FROM bins WHERE bin_id IS NOT NULL
                  UNION SELECT bin_id FROM proteins WHERE bin_id IS NOT NULL)""")
        n_genomes = con.execute("SELECT count(*) FROM genome_dictionary").fetchone()[0]
        n_declared = con.execute("SELECT count(DISTINCT bin_id) FROM bins WHERE bin_id IS NOT NULL").fetchone()[0]
        outside = con.execute("""SELECT count(DISTINCT p.bin_id) FROM proteins p
            LEFT JOIN bins b USING(bin_id) WHERE p.bin_id IS NOT NULL AND b.bin_id IS NULL""").fetchone()[0]
        positions = [0]
        con.execute("SELECT bin_id FROM genome_dictionary ORDER BY gid")
        with (out / "genomes.utf8").open("wb") as blob:
            while rows := con.fetchmany(65536):
                for (name,) in rows:
                    raw = name.encode("utf-8")
                    blob.write(raw)
                    positions.append(positions[-1] + len(raw))
        np.asarray(positions, dtype="<u8").tofile(out / "genomes.offsets.u64")
        counts = {"associated_proteins": 0, "list_bins": 0, "roaring_bins": 0, "empty_bins": 0}
        sql = """SELECT gd.gid,pd.pid FROM proteins p JOIN protein_dictionary pd USING(protein_id)
               JOIN genome_dictionary gd ON p.bin_id=gd.bin_id ORDER BY gd.gid,pd.pid"""
        con.execute(sql)
        with (out / "postings.bin").open("wb") as postings, (out / "catalog.bin").open("wb") as catalog:
            next_gid = 0

            def emit(gid, ids):
                nonlocal next_gid
                if gid != next_gid:
                    raise BuildError("Scope stream genome IDs must be dense and ordered")
                raw = np.asarray(ids, dtype="<u4").tobytes()
                bitmap = BitMap(ids)
                bitmap.run_optimize()
                packed = bitmap.serialize()
                encoding, payload = (1, packed) if len(packed) < len(raw) else (0, raw)
                catalog.write(struct.pack("<QQIB3x", postings.tell(), len(payload), len(ids), encoding))
                postings.write(payload)
                counts["associated_proteins"] += len(ids)
                counts["roaring_bins" if encoding else "list_bins"] += 1
                counts["empty_bins"] += int(not ids)
                next_gid += 1

            group, values = None, array.array("I")
            while rows := con.fetchmany(65536):
                for gid, pid in rows:
                    if gid != group:
                        if group is not None:
                            emit(group, values)
                        while next_gid < gid:
                            emit(next_gid, array.array("I"))
                        group, values = gid, array.array("I")
                    values.append(pid)
            if group is not None:
                emit(group, values)
            while next_gid < n_genomes:
                emit(next_gid, array.array("I"))
        null_owners = BitMap()
        con.execute("""SELECT pd.pid FROM proteins p JOIN protein_dictionary pd USING(protein_id)
                       WHERE p.bin_id IS NULL ORDER BY pd.pid""")
        while rows := con.fetchmany(65536):
            null_owners.update(pid for (pid,) in rows)
        null_owners.run_optimize()
        (out / "null_ownership.roaring").write_bytes(null_owners.serialize())
        if counts["associated_proteins"] + len(null_owners) != n:
            raise BuildError("Genome scope associations fail the exhaustive ownership count")
    after = _finish_source(source, before, source_hash)
    files = _payloads(out)
    manifest = {"format": SCOPE_FORMAT,
                "source": {"name": source.name, "sha256": source_hash, "before": before, "after": after},
                "matching_manifest_sha256": sha256_file(matching / "manifest.json"),
                "shared_protein_dictionary": {"protein_count": n, "files": shared},
                "genome_count": n_genomes, "declared_source_bin_count": n_declared,
                "owners_outside_declared_bin_catalog": outside, "null_ownership_count": len(null_owners),
                "ownership_semantics": "NULL owners in separate bitmap; empty string remains a non-NULL bin ID",
                "genome_dictionary": "UTF-8 byte ordered non-NULL union of bins catalog and observed owners",
                "catalog": "<QQIB3x offset,bytecount,cardinality,encoding(0=LEu32 list,1=portableRoaring)",
                "scope_semantics": "union selected known non-NULL bins; unknown IDs and empty request select no proteins",
                "sql": sql, **counts, "files": files, "payload_bytes": sum(x["bytes"] for x in files.values()),
                "versions": _versions(), "builder_sha256": sha256_file(__file__),
                "wall_seconds": time.perf_counter() - started, "peak_rss_bytes": _peak_rss_bytes()}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest
