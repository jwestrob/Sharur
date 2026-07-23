"""DuckDB storage backend (Part 8.2)."""

from __future__ import annotations

import re
from pathlib import Path
from typing import TYPE_CHECKING

import duckdb

from sharur.core.types import Annotation, Protein
from sharur.storage.migrations import run_migrations

from .schema import SCHEMA


if TYPE_CHECKING:
    import pandas as pd


class DuckDBStore:
    """DuckDB storage backend."""

    def __init__(
        self,
        db_path: Path | None = None,
        *,
        read_only: bool = False,
        threads: int | None = None,
        memory_limit: str | None = None,
        temp_directory: Path | str | None = None,
    ):
        """
        Initialize store.

        If db_path is None, uses in-memory database.
        """
        self.db_path = Path(db_path) if db_path else None
        self.read_only = read_only
        if threads is not None and (
            isinstance(threads, bool) or not isinstance(threads, int) or threads < 1
        ):
            raise ValueError("threads must be a positive integer")
        if memory_limit is not None:
            normalized_memory = memory_limit.strip()
            if not re.fullmatch(
                r"(?:\d+(?:\.\d+)?)\s*(?:B|KB|MB|GB|TB|KiB|MiB|GiB|TiB)",
                normalized_memory,
                flags=re.IGNORECASE,
            ):
                raise ValueError(
                    "memory_limit must include a supported unit, for example '4GB'"
                )
            self.memory_limit = normalized_memory
        else:
            self.memory_limit = None
        self.threads = threads
        self.temp_directory = (
            Path(temp_directory).expanduser().resolve()
            if temp_directory is not None
            else None
        )
        self._conn: duckdb.DuckDBPyConnection | None = None

    # ------------------------------------------------------------------ #
    # Connection management
    # ------------------------------------------------------------------ #
    @property
    def conn(self) -> duckdb.DuckDBPyConnection:
        """Lazy connection initialization."""
        if self._conn is None:
            if self.db_path:
                self._conn = duckdb.connect(
                    str(self.db_path),
                    read_only=self.read_only,
                )
            else:
                self._conn = duckdb.connect(":memory:")

            self._configure_resources()
            if not self.read_only:
                self._initialize_schema()
        return self._conn

    def _configure_resources(self) -> None:
        """Apply per-process limits before any analytical work begins."""

        assert self._conn is not None
        if self.threads is not None:
            self._conn.execute("SET threads = ?", [self.threads])
        if self.memory_limit is not None:
            self._conn.execute("SET memory_limit = ?", [self.memory_limit])
        if self.temp_directory is not None:
            self.temp_directory.mkdir(parents=True, exist_ok=True)
            self._conn.execute("SET temp_directory = ?", [str(self.temp_directory)])

    @property
    def resource_budget(self) -> dict[str, object | None]:
        """Return the explicitly configured per-connection resource budget."""

        return {
            "threads": self.threads,
            "memory_limit": self.memory_limit,
            "temp_directory": (
                str(self.temp_directory) if self.temp_directory is not None else None
            ),
        }

    def _initialize_schema(self) -> None:
        """Create tables if they don't exist, then run pending migrations."""
        self._conn.execute(SCHEMA)
        run_migrations(self._conn)

    def close(self) -> None:
        """Close the lazy connection without opening one solely to close it."""
        if self._conn is not None:
            self._conn.close()
            self._conn = None

    def __enter__(self) -> DuckDBStore:
        return self

    def __exit__(self, *_exc_info: object) -> None:
        self.close()

    # ------------------------------------------------------------------ #
    # Generic execution helpers
    # ------------------------------------------------------------------ #
    def execute(self, query: str, params: dict | list | tuple | None = None) -> list[tuple]:
        """Execute query and return results."""
        if params:
            return self.conn.execute(query, params).fetchall()
        return self.conn.execute(query).fetchall()

    def execute_df(self, query: str, params: dict | list | tuple | None = None) -> pd.DataFrame:
        """Execute query and return DataFrame."""
        if params:
            return self.conn.execute(query, params).df()
        return self.conn.execute(query).df()

    # ------------------------------------------------------------------ #
    # Convenience methods
    # ------------------------------------------------------------------ #
    def get_protein(self, protein_id: str) -> Protein | None:
        """Get single protein by ID."""
        row = self.conn.execute(
            """
            SELECT protein_id, contig_id, bin_id, start, end_coord, strand, sequence
            FROM proteins
            WHERE protein_id = ?
            """,
            [protein_id],
        ).fetchone()

        if not row:
            return None
        return Protein(
            protein_id=row[0],
            contig_id=row[1],
            bin_id=row[2],
            start=row[3],
            end=row[4],
            strand=row[5],
            sequence=row[6],
        )

    def get_proteins(self, protein_ids: list[str]) -> list[Protein]:
        """Get multiple proteins by ID."""
        if not protein_ids:
            return []

        placeholders = ",".join(["?"] * len(protein_ids))
        rows = self.conn.execute(
            f"""
            SELECT protein_id, contig_id, bin_id, start, end_coord, strand, sequence
            FROM proteins
            WHERE protein_id IN ({placeholders})
            """,
            protein_ids,
        ).fetchall()

        return [
            Protein(
                protein_id=r[0],
                contig_id=r[1],
                bin_id=r[2],
                start=r[3],
                end=r[4],
                strand=r[5],
                sequence=r[6],
            )
            for r in rows
        ]

    def get_proteins_in_window(
        self,
        contig_id: str,
        start: int,
        end: int,
        *,
        bin_id: str | None = None,
    ) -> list[Protein]:
        """Get proteins overlapping a coordinate window.

        ``contig_id`` is not guaranteed to be globally unique in legacy
        datasets. Supply ``bin_id`` when known. If it is omitted and the
        protein table associates the contig label with multiple bins, fail
        closed instead of returning a cross-genome mixture.
        """
        if bin_id is None:
            bins = self.conn.execute(
                """
                SELECT DISTINCT bin_id
                FROM proteins
                WHERE contig_id = ?
                ORDER BY bin_id
                """,
                [contig_id],
            ).fetchall()
            if len(bins) > 1:
                raise ValueError(
                    f"Contig ID {contig_id!r} occurs in multiple bins; "
                    "pass bin_id explicitly"
                )
            if bins:
                bin_id = bins[0][0]

        rows = self.conn.execute(
            """
            SELECT protein_id, contig_id, bin_id, start, end_coord, strand, sequence
            FROM proteins
            WHERE contig_id = ?
              AND (? IS NULL OR bin_id = ?)
              AND NOT (end_coord <= ? OR start >= ?)
            ORDER BY start
            """,
            [contig_id, bin_id, bin_id, start, end],
        ).fetchall()

        return [
            Protein(
                protein_id=r[0],
                contig_id=r[1],
                bin_id=r[2],
                start=r[3],
                end=r[4],
                strand=r[5],
                sequence=r[6],
            )
            for r in rows
        ]

    def get_annotations(self, protein_id: str) -> list[Annotation]:
        """Get all annotations for a protein."""
        rows = self.conn.execute(
            """
            SELECT protein_id, source, accession, name, description,
                   evalue, score, start_aa, end_aa
            FROM annotations
            WHERE protein_id = ?
            ORDER BY start_aa NULLS FIRST
            """,
            [protein_id],
        ).fetchall()

        return [
            Annotation(
                protein_id=r[0],
                source=r[1],
                accession=r[2],
                name=r[3],
                description=r[4],
                evalue=r[5],
                score=r[6],
                start_aa=r[7],
                end_aa=r[8],
            )
            for r in rows
        ]

    def search_annotations(
        self,
        source: str | None = None,
        accession: str | None = None,
        name_pattern: str | None = None,
    ) -> list[tuple[str, Annotation]]:
        """Search annotations, returns (protein_id, annotation) pairs."""
        clauses: list[str] = []
        params: list = []

        if source:
            clauses.append("source = ?")
            params.append(source)
        if accession:
            clauses.append("accession = ?")
            params.append(accession)
        if name_pattern:
            clauses.append("(LOWER(name) LIKE LOWER(?) OR LOWER(description) LIKE LOWER(?))")
            like = f"%{name_pattern}%"
            params.extend([like, like])

        where = f"WHERE {' AND '.join(clauses)}" if clauses else ""
        query = f"""
            SELECT protein_id, source, accession, name, description,
                   evalue, score, start_aa, end_aa
            FROM annotations
            {where}
        """

        rows = self.conn.execute(query, params if params else None).fetchall()

        return [
            (
                r[0],
                Annotation(
                    protein_id=r[0],
                    source=r[1],
                    accession=r[2],
                    name=r[3],
                    description=r[4],
                    evalue=r[5],
                    score=r[6],
                    start_aa=r[7],
                    end_aa=r[8],
                ),
            )
            for r in rows
        ]


__all__ = ["DuckDBStore"]
