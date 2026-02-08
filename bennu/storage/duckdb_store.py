"""DuckDB storage backend (Part 8.2)."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import duckdb
import pandas as pd

from bennu.core.types import Annotation, Protein

from .schema import SCHEMA


class DuckDBStore:
    """DuckDB storage backend."""

    def __init__(self, db_path: Optional[Path] = None):
        """
        Initialize store.

        If db_path is None, uses in-memory database.
        """
        self.db_path = Path(db_path) if db_path else None
        self._conn: Optional[duckdb.DuckDBPyConnection] = None

    # ------------------------------------------------------------------ #
    # Connection management
    # ------------------------------------------------------------------ #
    @property
    def conn(self) -> duckdb.DuckDBPyConnection:
        """Lazy connection initialization."""
        if self._conn is None:
            self._conn = duckdb.connect(str(self.db_path) if self.db_path else ":memory:")
            self._initialize_schema()
        return self._conn

    def _initialize_schema(self) -> None:
        """Create tables if they don't exist, then run pending migrations."""
        self._conn.execute(SCHEMA)
        from bennu.storage.migrations import get_current_version, run_migrations

        current = get_current_version(self._conn)
        if current == 0:
            run_migrations(self._conn)

    # ------------------------------------------------------------------ #
    # Generic execution helpers
    # ------------------------------------------------------------------ #
    def execute(self, query: str, params: Optional[dict | list | tuple] = None) -> list[tuple]:
        """Execute query and return results."""
        if params:
            return self.conn.execute(query, params).fetchall()
        return self.conn.execute(query).fetchall()

    def execute_df(self, query: str, params: Optional[dict | list | tuple] = None) -> pd.DataFrame:
        """Execute query and return DataFrame."""
        if params:
            return self.conn.execute(query, params).df()
        return self.conn.execute(query).df()

    # ------------------------------------------------------------------ #
    # Convenience methods
    # ------------------------------------------------------------------ #
    def get_protein(self, protein_id: str) -> Optional[Protein]:
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

    def get_proteins_in_window(self, contig_id: str, start: int, end: int) -> list[Protein]:
        """Get proteins overlapping coordinate window."""
        rows = self.conn.execute(
            """
            SELECT protein_id, contig_id, bin_id, start, end_coord, strand, sequence
            FROM proteins
            WHERE contig_id = ?
              AND NOT (end_coord <= ? OR start >= ?)
            ORDER BY start
            """,
            [contig_id, start, end],
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
        source: Optional[str] = None,
        accession: Optional[str] = None,
        name_pattern: Optional[str] = None,
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
