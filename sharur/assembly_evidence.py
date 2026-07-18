"""Optional contig-level assembly and host-assignment evidence.

Assembly evidence is intentionally stored beside, not inside, the canonical
Sharur DuckDB.  Existing datasets therefore remain valid and read-only
exploration does not acquire a mandatory schema migration.  The sidecar is
small for ordinary coverage/taxonomy metrics and never stores raw reads,
assemblies, or k-mer vectors.

Composition analysis is explicit opt-in because reading every assembly can be
expensive on shared filesystems.  When requested, vectors are held in memory
for one assembly and only scalar GC/tetranucleotide summaries are persisted.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import statistics
import uuid
from array import array
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING, Any

import duckdb

from sharur.core.case_models import AssemblyEvidenceRecord


if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator, Mapping


ASSEMBLY_EVIDENCE_SCHEMA_VERSION = 1
DEFAULT_ASSEMBLY_EVIDENCE_FILENAME = "assembly_evidence.duckdb"

_EVIDENCE_FIELDS: dict[str, str] = {
    "bin_id": "VARCHAR",
    "contig_id": "VARCHAR",
    "coverage_mean": "DOUBLE",
    "coverage_sd": "DOUBLE",
    "coverage_cv": "DOUBLE",
    "coverage_ratio_to_bin_median": "DOUBLE",
    "mapped_reads": "BIGINT",
    "proper_pair_fraction": "DOUBLE",
    "insert_size_median": "DOUBLE",
    "insert_size_mad": "DOUBLE",
    "snv_count": "BIGINT",
    "snv_density_per_kb": "DOUBLE",
    "assembly_graph_degree": "INTEGER",
    "assembly_graph_component": "VARCHAR",
    "taxonomy": "VARCHAR",
    "taxonomy_method": "VARCHAR",
    "taxonomy_confidence": "DOUBLE",
    "taxonomy_congruent": "BOOLEAN",
    "gc_zscore": "DOUBLE",
    "tetranucleotide_distance": "DOUBLE",
    "tetranucleotide_percentile": "DOUBLE",
    "source": "VARCHAR",
    "metadata": "JSON",
}

_FLOAT_FIELDS = {name for name, sql_type in _EVIDENCE_FIELDS.items() if sql_type == "DOUBLE"}
_INTEGER_FIELDS = {
    name for name, sql_type in _EVIDENCE_FIELDS.items() if sql_type in {"INTEGER", "BIGINT"}
}
_BOOLEAN_FIELDS = {name for name, sql_type in _EVIDENCE_FIELDS.items() if sql_type == "BOOLEAN"}


def default_assembly_evidence_path(db_path: str | Path) -> Path:
    """Return the conventional sidecar path for a Sharur database."""
    return Path(db_path).expanduser().resolve().parent / DEFAULT_ASSEMBLY_EVIDENCE_FILENAME


def discover_assembly_evidence(
    db_path: str | Path,
    explicit_path: str | Path | None = None,
) -> Path | None:
    """Return an explicitly configured or conventionally located sidecar."""
    candidate = (
        Path(explicit_path).expanduser().resolve()
        if explicit_path is not None
        else default_assembly_evidence_path(db_path)
    )
    return candidate if candidate.is_file() else None


class AssemblyEvidenceStore:
    """Small DuckDB sidecar for optional contig evidence."""

    def __init__(self, path: str | Path, *, read_only: bool = False):
        self.path = Path(path).expanduser().resolve()
        self.read_only = read_only
        self._conn: duckdb.DuckDBPyConnection | None = None

    @property
    def conn(self) -> duckdb.DuckDBPyConnection:
        if self._conn is None:
            if not self.read_only:
                self.path.parent.mkdir(parents=True, exist_ok=True)
            self._conn = duckdb.connect(str(self.path), read_only=self.read_only)
            if not self.read_only:
                self._initialize()
        return self._conn

    def _initialize(self) -> None:
        columns = ",\n".join(
            f"            {name} {sql_type}" for name, sql_type in _EVIDENCE_FIELDS.items()
        )
        self._conn.execute(
            f"""
            CREATE TABLE IF NOT EXISTS contig_evidence (
{columns},
                updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                PRIMARY KEY (bin_id, contig_id)
            );
            CREATE INDEX IF NOT EXISTS idx_contig_evidence_contig
                ON contig_evidence(contig_id);
            CREATE INDEX IF NOT EXISTS idx_contig_evidence_bin
                ON contig_evidence(bin_id);

            CREATE TABLE IF NOT EXISTS evidence_imports (
                import_id VARCHAR PRIMARY KEY,
                imported_at TIMESTAMP,
                input_path VARCHAR,
                input_sha256 VARCHAR,
                row_count BIGINT,
                options JSON
            );

            CREATE TABLE IF NOT EXISTS evidence_schema_version (
                version INTEGER PRIMARY KEY,
                applied_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            );
            INSERT OR IGNORE INTO evidence_schema_version(version)
                VALUES ({ASSEMBLY_EVIDENCE_SCHEMA_VERSION});
            """
        )

    def close(self) -> None:
        if self._conn is not None:
            self._conn.close()
            self._conn = None

    def __enter__(self) -> AssemblyEvidenceStore:
        _ = self.conn
        return self

    def __exit__(self, *_exc_info: object) -> None:
        self.close()

    def count(self) -> int:
        return int(self.conn.execute("SELECT COUNT(*) FROM contig_evidence").fetchone()[0])

    def schema_version(self) -> int:
        row = self.conn.execute("SELECT MAX(version) FROM evidence_schema_version").fetchone()
        return int(row[0]) if row and row[0] is not None else 0

    def get(self, bin_id: str, contig_id: str) -> AssemblyEvidenceRecord | None:
        """Get one record and derive its coverage ratio when possible."""
        columns = list(_EVIDENCE_FIELDS)
        row = self.conn.execute(
            f"""
            SELECT {", ".join(columns)}
            FROM contig_evidence
            WHERE bin_id = ? AND contig_id = ?
            """,
            [bin_id, contig_id],
        ).fetchone()
        if row is None:
            return None
        payload = dict(zip(columns, row, strict=True))
        payload["metadata"] = _decode_json(payload.get("metadata"))

        if (
            payload.get("coverage_mean") not in {None, 0}
            and payload.get("coverage_sd") is not None
            and payload.get("coverage_cv") is None
        ):
            payload["coverage_cv"] = float(payload["coverage_sd"]) / float(payload["coverage_mean"])
        if (
            payload.get("coverage_mean") is not None
            and payload.get("coverage_ratio_to_bin_median") is None
        ):
            median_row = self.conn.execute(
                """
                SELECT MEDIAN(coverage_mean)
                FROM contig_evidence
                WHERE bin_id = ? AND coverage_mean IS NOT NULL
                """,
                [bin_id],
            ).fetchone()
            median = median_row[0] if median_row else None
            if median not in {None, 0}:
                payload["coverage_ratio_to_bin_median"] = float(payload["coverage_mean"]) / float(
                    median
                )
        return AssemblyEvidenceRecord(**payload)

    def upsert(
        self,
        records: Iterable[AssemblyEvidenceRecord | Mapping[str, Any]],
        *,
        source_path: str | Path | None = None,
        source_sha256: str | None = None,
        options: Mapping[str, Any] | None = None,
    ) -> int:
        """Upsert normalized records in one transaction."""
        normalized = [
            record
            if isinstance(record, AssemblyEvidenceRecord)
            else AssemblyEvidenceRecord(**dict(record))
            for record in records
        ]
        if not normalized:
            return 0

        columns = list(_EVIDENCE_FIELDS)
        placeholders = ", ".join(["?"] * len(columns))
        keys = [(record.bin_id, record.contig_id) for record in normalized]
        if len(set(keys)) != len(keys):
            raise ValueError("One upsert batch contains duplicate (bin_id, contig_id) rows")
        update_columns = [column for column in columns if column not in {"bin_id", "contig_id"}]
        update_sql = ", ".join(f"{column} = excluded.{column}" for column in update_columns)

        self.conn.execute("BEGIN TRANSACTION")
        try:
            existing: dict[tuple[str, str], dict[str, Any]] = {}
            for start in range(0, len(keys), 5_000):
                batch = keys[start : start + 5_000]
                key_placeholders = ", ".join(["(?, ?)"] * len(batch))
                values = [value for key in batch for value in key]
                current_rows = self.conn.execute(
                    f"""
                    SELECT {", ".join(columns)}
                    FROM contig_evidence
                    WHERE (bin_id, contig_id) IN ({key_placeholders})
                    """,
                    values,
                ).fetchall()
                for current_row in current_rows:
                    current = dict(zip(columns, current_row, strict=True))
                    current["metadata"] = _decode_json(current.get("metadata"))
                    existing[(str(current["bin_id"]), str(current["contig_id"]))] = current

            rows = []
            for record in normalized:
                payload = record.model_dump()
                current = existing.get((record.bin_id, record.contig_id))
                if current is not None:
                    for column in update_columns:
                        if column == "metadata":
                            payload[column] = {
                                **current.get("metadata", {}),
                                **(payload.get("metadata") or {}),
                            }
                        elif payload.get(column) is None:
                            payload[column] = current.get(column)
                payload["metadata"] = (
                    json.dumps(payload["metadata"], sort_keys=True)
                    if payload.get("metadata")
                    else None
                )
                rows.append([payload.get(column) for column in columns])

            self.conn.executemany(
                f"""
                INSERT INTO contig_evidence ({", ".join(columns)})
                VALUES ({placeholders})
                ON CONFLICT (bin_id, contig_id) DO UPDATE SET
                    {update_sql},
                    updated_at = now()
                """,
                rows,
            )
            self.conn.execute(
                """
                INSERT INTO evidence_imports (
                    import_id, imported_at, input_path, input_sha256,
                    row_count, options
                ) VALUES (?, ?, ?, ?, ?, ?)
                """,
                [
                    uuid.uuid4().hex,
                    datetime.now(timezone.utc),
                    str(source_path) if source_path is not None else None,
                    source_sha256,
                    len(rows),
                    json.dumps(dict(options or {}), sort_keys=True),
                ],
            )
            self.conn.commit()
        except Exception:
            self.conn.rollback()
            raise
        return len(rows)


def _decode_json(value: Any) -> dict[str, Any]:
    if value is None:
        return {}
    if isinstance(value, dict):
        return value
    try:
        decoded = json.loads(str(value))
    except (TypeError, json.JSONDecodeError):
        return {"raw": str(value)}
    return decoded if isinstance(decoded, dict) else {"value": decoded}


def _coerce_value(field: str, value: Any) -> Any:
    if value is None:
        return None
    if isinstance(value, str):
        value = value.strip()
        if not value or value.lower() in {"na", "nan", "none", "null", "."}:
            return None
    if field in _FLOAT_FIELDS:
        return float(value)
    if field in _INTEGER_FIELDS:
        return int(float(value))
    if field in _BOOLEAN_FIELDS:
        if isinstance(value, bool):
            return value
        normalized = str(value).strip().lower()
        if normalized in {"1", "true", "yes", "y", "congruent", "pass"}:
            return True
        if normalized in {"0", "false", "no", "n", "discordant", "fail"}:
            return False
        raise ValueError(f"Cannot parse boolean value {value!r} for {field}")
    if field == "metadata":
        return _decode_json(value)
    return str(value)


def _read_tabular_records(path: Path) -> Iterator[dict[str, Any]]:
    suffix = path.suffix.lower()
    if suffix in {".jsonl", ".ndjson"}:
        with path.open(encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.strip():
                    continue
                value = json.loads(line)
                if not isinstance(value, dict):
                    raise ValueError(f"{path}:{line_number} is not a JSON object")
                yield value
        return

    delimiter = "\t" if suffix in {".tsv", ".tab"} else ","
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError(f"No header found in {path}")
        yield from reader


def _sha256_file(path: Path, *, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def import_contig_evidence(
    input_path: str | Path,
    sidecar_path: str | Path,
    *,
    source: str | None = None,
    validate_db_path: str | Path | None = None,
    hash_input: bool = True,
) -> dict[str, Any]:
    """Import TSV/CSV/JSONL contig evidence into a compact sidecar.

    Required columns are ``bin_id`` and ``contig_id``.  Recognized typed
    columns are listed in :data:`_EVIDENCE_FIELDS`; unrecognized columns are
    retained under ``metadata``.
    """
    input_path = Path(input_path).expanduser().resolve()
    if not input_path.is_file():
        raise FileNotFoundError(input_path)
    resolved_sidecar = Path(sidecar_path).expanduser().resolve()
    if (
        validate_db_path is not None
        and resolved_sidecar == Path(validate_db_path).expanduser().resolve()
    ):
        raise ValueError("Assembly evidence must use a separate sidecar, not the core DuckDB")

    records: list[AssemblyEvidenceRecord] = []
    seen_keys: set[tuple[str, str]] = set()
    for row_number, raw in enumerate(_read_tabular_records(input_path), start=2):
        if not raw.get("bin_id") or not raw.get("contig_id"):
            raise ValueError(f"{input_path}:{row_number} requires non-empty bin_id and contig_id")
        key = (str(raw["bin_id"]), str(raw["contig_id"]))
        if key in seen_keys:
            raise ValueError(f"{input_path}:{row_number} repeats (bin_id, contig_id) {key}")
        seen_keys.add(key)
        payload: dict[str, Any] = {}
        metadata = _decode_json(raw.get("metadata"))
        for key, value in raw.items():
            if key in _EVIDENCE_FIELDS:
                payload[key] = _coerce_value(key, value)
            elif value not in {None, ""}:
                metadata[key] = value
        payload["metadata"] = metadata
        if source is not None:
            payload["source"] = source
            payload["metadata"]["import_source"] = source
        records.append(AssemblyEvidenceRecord(**payload))
    if not records:
        raise ValueError(f"No evidence records found in {input_path}")

    validation = _validate_core_contigs(records, validate_db_path)
    digest = _sha256_file(input_path) if hash_input else None
    with AssemblyEvidenceStore(sidecar_path) as evidence_store:
        written = evidence_store.upsert(
            records,
            source_path=input_path,
            source_sha256=digest,
            options={
                "source": source,
                "validated_against": (
                    str(Path(validate_db_path).expanduser().resolve())
                    if validate_db_path is not None
                    else None
                ),
                "hash_input": hash_input,
            },
        )
        total = evidence_store.count()

    return {
        "input_path": str(input_path),
        "sidecar_path": str(Path(sidecar_path).expanduser().resolve()),
        "input_sha256": digest,
        "rows_written": written,
        "total_rows": total,
        "validation": validation,
    }


def _validate_core_contigs(
    records: Iterable[AssemblyEvidenceRecord],
    db_path: str | Path | None,
) -> dict[str, Any]:
    if db_path is None:
        return {"state": "not_requested", "missing": []}

    db_path = Path(db_path).expanduser().resolve()
    if not db_path.is_file():
        raise FileNotFoundError(db_path)
    keys = sorted({(record.bin_id, record.contig_id) for record in records})
    connection = duckdb.connect(str(db_path), read_only=True)
    try:
        missing: list[dict[str, str]] = []
        for start in range(0, len(keys), 5_000):
            batch = keys[start : start + 5_000]
            placeholders = ", ".join(["(?, ?)"] * len(batch))
            values = [value for key in batch for value in key]
            rows = connection.execute(
                f"""
                WITH requested(bin_id, contig_id) AS (
                    VALUES {placeholders}
                )
                SELECT requested.bin_id, requested.contig_id
                FROM requested
                LEFT JOIN proteins
                  ON proteins.bin_id = requested.bin_id
                 AND proteins.contig_id = requested.contig_id
                GROUP BY requested.bin_id, requested.contig_id
                HAVING COUNT(proteins.protein_id) = 0
                ORDER BY requested.bin_id, requested.contig_id
                """,
                values,
            ).fetchall()
            missing.extend(
                {"bin_id": str(bin_id), "contig_id": str(contig_id)} for bin_id, contig_id in rows
            )
    finally:
        connection.close()
    if missing:
        preview = ", ".join(f"({item['bin_id']}, {item['contig_id']})" for item in missing[:5])
        raise ValueError(
            f"{len(missing)} evidence rows do not match a core (bin_id, contig_id), "
            f"including {preview}"
        )
    return {"state": "validated", "missing": []}


@contextmanager
def _fasta_records(path: Path) -> Iterator[Iterator[tuple[str, str]]]:
    """Yield a streaming FASTA iterator while keeping the handle scoped."""
    with (
        gzip.open(path, mode="rt", encoding="utf-8")
        if path.suffix.lower() == ".gz"
        else path.open(encoding="utf-8")
    ) as handle:

        def iterator() -> Iterator[tuple[str, str]]:
            identifier: str | None = None
            chunks: list[str] = []
            for line_number, line in enumerate(handle, start=1):
                if line.startswith(">"):
                    if identifier is not None:
                        yield identifier, "".join(chunks).upper()
                    header = line[1:].strip()
                    if not header:
                        raise ValueError(f"{path}:{line_number} has an empty FASTA header")
                    identifier = header.split()[0]
                    chunks = []
                elif line.strip():
                    if identifier is None:
                        raise ValueError(
                            f"{path}:{line_number} contains sequence before a FASTA header"
                        )
                    chunks.append(line.strip())
            if identifier is not None:
                yield identifier, "".join(chunks).upper()

        yield iterator()


_BASE_INDEX = {"A": 0, "C": 1, "G": 2, "T": 3}


def _reverse_complement_code(code: int) -> int:
    reverse_complement = 0
    for _ in range(4):
        reverse_complement = (reverse_complement << 2) | (3 - (code & 0b11))
        code >>= 2
    return reverse_complement


_CANONICAL_KMER_INDEX = tuple(min(code, _reverse_complement_code(code)) for code in range(256))


def _tetranucleotide_counts(sequence: str) -> array:
    """Count strand-invariant canonical 4-mers in a fixed 256-slot vector."""
    counts = array("I", [0]) * 256
    rolling = 0
    valid = 0
    for base in sequence:
        code = _BASE_INDEX.get(base)
        if code is None:
            rolling = 0
            valid = 0
            continue
        rolling = ((rolling << 2) | code) & 0xFF
        valid += 1
        if valid >= 4:
            counts[_CANONICAL_KMER_INDEX[rolling]] += 1
    return counts


def _cosine_distance(left: array, right: array) -> float | None:
    dot = sum(int(a) * int(b) for a, b in zip(left, right, strict=True))
    left_norm = math.sqrt(sum(int(value) ** 2 for value in left))
    right_norm = math.sqrt(sum(int(value) ** 2 for value in right))
    if left_norm == 0 or right_norm == 0:
        return None
    return 1.0 - dot / (left_norm * right_norm)


def compute_composition_evidence(
    assemblies: Mapping[str, str | Path],
    sidecar_path: str | Path,
    *,
    validate_db_path: str | Path | None = None,
    source: str = "sharur_composition_opt_in",
) -> dict[str, Any]:
    """Explicitly scan assemblies and persist scalar composition summaries.

    ``assemblies`` maps ``bin_id`` to FASTA path.  Full tetranucleotide vectors
    are kept only in memory for the current assembly and are discarded after
    scalar cosine distances and within-assembly percentiles are computed.
    """
    if not assemblies:
        raise ValueError("At least one bin_id-to-FASTA mapping is required")
    resolved_sidecar = Path(sidecar_path).expanduser().resolve()
    if (
        validate_db_path is not None
        and resolved_sidecar == Path(validate_db_path).expanduser().resolve()
    ):
        raise ValueError("Assembly evidence must use a separate sidecar, not the core DuckDB")
    all_records: list[AssemblyEvidenceRecord] = []
    assembly_summaries: list[dict[str, Any]] = []

    for bin_id, raw_path in assemblies.items():
        fasta_path = Path(raw_path).expanduser().resolve()
        if not fasta_path.is_file():
            raise FileNotFoundError(fasta_path)

        vectors: dict[str, array] = {}
        gc_values: dict[str, float | None] = {}
        assembly_vector = array("Q", [0]) * 256
        content_digest = hashlib.sha256()
        with _fasta_records(fasta_path) as records:
            for contig_id, sequence in records:
                if contig_id in vectors:
                    raise ValueError(f"Duplicate FASTA identifier {contig_id!r} in {fasta_path}")
                if not sequence:
                    raise ValueError(f"Empty FASTA record {contig_id!r} in {fasta_path}")
                content_digest.update(f">{contig_id}\n".encode())
                content_digest.update(sequence.encode())
                content_digest.update(b"\n")
                vector = _tetranucleotide_counts(sequence)
                vectors[contig_id] = vector
                for index, count in enumerate(vector):
                    assembly_vector[index] += count
                canonical = sum(sequence.count(base) for base in "ACGT")
                gc_values[contig_id] = (
                    (sequence.count("G") + sequence.count("C")) / canonical if canonical else None
                )

        if not vectors:
            raise ValueError(f"No FASTA records found in {fasta_path}")

        populated_gc_values = [value for value in gc_values.values() if value is not None]
        gc_mean = statistics.fmean(populated_gc_values) if populated_gc_values else None
        gc_sd = statistics.pstdev(populated_gc_values) if populated_gc_values else None
        distances: dict[str, float | None] = {}
        for contig_id, vector in vectors.items():
            if len(vectors) == 1:
                distances[contig_id] = None
                continue
            leave_one_out = array(
                "Q",
                (int(assembly_vector[index]) - int(vector[index]) for index in range(256)),
            )
            distances[contig_id] = _cosine_distance(vector, leave_one_out)
        sorted_distances = sorted(
            distance for distance in distances.values() if distance is not None
        )

        for contig_id, distance in distances.items():
            if distance is None:
                percentile = None
            else:
                lower = sum(value < distance for value in sorted_distances)
                equal = sum(value == distance for value in sorted_distances)
                percentile = 100.0 * (lower + 0.5 * equal) / len(sorted_distances)
            gc_value = gc_values[contig_id]
            gc_zscore = (
                None
                if gc_value is None or gc_mean is None
                else (gc_value - gc_mean) / gc_sd
                if gc_sd
                else 0.0
            )
            all_records.append(
                AssemblyEvidenceRecord(
                    bin_id=str(bin_id),
                    contig_id=contig_id,
                    gc_zscore=gc_zscore,
                    tetranucleotide_distance=distance,
                    tetranucleotide_percentile=percentile,
                    source=source,
                    metadata={
                        "composition_fasta": str(fasta_path),
                        "composition_content_sha256": content_digest.hexdigest(),
                        "composition_hash_basis": ("normalized_fasta_ids_and_uppercase_sequences"),
                        "composition_reference": "leave_one_contig_out_bin_aggregate",
                        "composition_reference_weighting": "tetranucleotide_count_weighted",
                        "kmer_canonicalization": "reverse_complement_collapsed",
                        "percentile_method": "within_bin_midrank",
                        "composition_source": source,
                        "vectors_persisted": False,
                    },
                )
            )
        assembly_summaries.append(
            {
                "bin_id": str(bin_id),
                "fasta": str(fasta_path),
                "content_sha256": content_digest.hexdigest(),
                "hash_basis": "normalized_fasta_ids_and_uppercase_sequences",
                "contigs": len(vectors),
                "vectors_persisted": False,
            }
        )

    validation = _validate_core_contigs(all_records, validate_db_path)
    with AssemblyEvidenceStore(sidecar_path) as evidence_store:
        written = evidence_store.upsert(
            all_records,
            options={
                "operation": "compute_composition",
                "assemblies": assembly_summaries,
                "vectors_persisted": False,
            },
        )
        total = evidence_store.count()
    return {
        "sidecar_path": str(Path(sidecar_path).expanduser().resolve()),
        "rows_written": written,
        "total_rows": total,
        "assemblies": assembly_summaries,
        "validation": validation,
    }


__all__ = [
    "ASSEMBLY_EVIDENCE_SCHEMA_VERSION",
    "DEFAULT_ASSEMBLY_EVIDENCE_FILENAME",
    "AssemblyEvidenceStore",
    "compute_composition_evidence",
    "default_assembly_evidence_path",
    "discover_assembly_evidence",
    "import_contig_evidence",
]
