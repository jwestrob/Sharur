"""Read-only access to normalized, run-scoped ELSA result sidecars."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

import duckdb


if TYPE_CHECKING:
    from collections.abc import Sequence


EXPECTED_SCHEMA_VERSION = 1
DEFAULT_SIDECAR_NAME = "synteny.duckdb"
_REQUIRED_RELATIONS = (
    "elsa_schema_version",
    "elsa_runs",
    "elsa_dataset_state",
    "elsa_genes",
    "elsa_blocks",
    "elsa_anchor_pairs",
    "elsa_clusters",
    "elsa_cluster_loci",
    "elsa_cluster_members",
    "current_elsa_run",
)


class SyntenySidecarError(RuntimeError):
    """Raised when a discovered ELSA sidecar violates its storage contract."""


class SyntenyDatasetMismatchError(SyntenySidecarError):
    """Raised when a sidecar belongs to a different Sharur dataset state."""


@dataclass(frozen=True)
class SyntenySidecarInspection:
    """Typed, non-mutating inspection of one ELSA sidecar."""

    path: Path
    state: str
    schema_version: int | None
    active_run_id: str | None
    dataset_id: str | None
    status: str | None
    gene_count: int
    block_count: int
    cluster_count: int
    singleton_count: int
    anchor_pair_count: int
    locus_count: int
    member_count: int
    error: str | None = None


@dataclass(frozen=True)
class SyntenyDatasetIdentityInspection:
    """Comparison between one sidecar run and the live Sharur dataset seal."""

    state: str
    compatible: bool
    sidecar_dataset_id: str | None
    seal_path: Path
    seal_dataset_id: str | None
    error: str | None = None


def discover_synteny_sidecar(
    db_path: str | Path,
    *,
    explicit_path: str | Path | None = None,
) -> Path | None:
    """Discover the canonical dataset-local ELSA result sidecar."""

    database = Path(db_path).expanduser().resolve()
    candidates = (
        [Path(explicit_path).expanduser().resolve()]
        if explicit_path is not None
        else [
            database.parent / DEFAULT_SIDECAR_NAME,
            database.parent / "synteny" / DEFAULT_SIDECAR_NAME,
        ]
    )
    return next((candidate for candidate in candidates if candidate.is_file()), None)


def inspect_synteny_sidecar(path: str | Path) -> SyntenySidecarInspection:
    """Inspect schema, active run, and counts without mutating the sidecar."""

    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file():
        return SyntenySidecarInspection(
            path=resolved,
            state="unavailable",
            schema_version=None,
            active_run_id=None,
            dataset_id=None,
            status=None,
            gene_count=0,
            block_count=0,
            cluster_count=0,
            singleton_count=0,
            anchor_pair_count=0,
            locus_count=0,
            member_count=0,
            error="sidecar does not exist",
        )
    try:
        connection = duckdb.connect(str(resolved), read_only=True)
        try:
            for relation in _REQUIRED_RELATIONS:
                escaped = relation.replace('"', '""')
                connection.execute(f'SELECT * FROM "{escaped}" LIMIT 0')
            version_row = connection.execute(
                "SELECT MAX(version) FROM elsa_schema_version"
            ).fetchone()
            version = int(version_row[0]) if version_row and version_row[0] is not None else None
            run = connection.execute(
                """
                SELECT
                    run_id, dataset_id, status, gene_count, block_count,
                    cluster_count, singleton_count, anchor_pair_count,
                    locus_count, member_count
                FROM current_elsa_run
                """
            ).fetchone()
            if run is None:
                return SyntenySidecarInspection(
                    path=resolved,
                    state="stale",
                    schema_version=version,
                    active_run_id=None,
                    dataset_id=None,
                    status=None,
                    gene_count=0,
                    block_count=0,
                    cluster_count=0,
                    singleton_count=0,
                    anchor_pair_count=0,
                    locus_count=0,
                    member_count=0,
                    error="sidecar has no active run",
                )
            state = (
                "available"
                if version == EXPECTED_SCHEMA_VERSION and run[2] == "ready"
                else "stale"
            )
            return SyntenySidecarInspection(
                path=resolved,
                state=state,
                schema_version=version,
                active_run_id=str(run[0]),
                dataset_id=str(run[1]) if run[1] is not None else None,
                status=str(run[2]) if run[2] is not None else None,
                gene_count=int(run[3] or 0),
                block_count=int(run[4] or 0),
                cluster_count=int(run[5] or 0),
                singleton_count=int(run[6] or 0),
                anchor_pair_count=int(run[7] or 0),
                locus_count=int(run[8] or 0),
                member_count=int(run[9] or 0),
            )
        finally:
            connection.close()
    except Exception as exc:
        return SyntenySidecarInspection(
            path=resolved,
            state="failed",
            schema_version=None,
            active_run_id=None,
            dataset_id=None,
            status=None,
            gene_count=0,
            block_count=0,
            cluster_count=0,
            singleton_count=0,
            anchor_pair_count=0,
            locus_count=0,
            member_count=0,
            error=f"{type(exc).__name__}: {exc}",
        )


def inspect_synteny_dataset_identity(
    db_path: str | Path,
    sidecar: str | Path | SyntenySidecarInspection,
) -> SyntenyDatasetIdentityInspection:
    """Compare sidecar identity to the live dataset seal, when available."""

    database = Path(db_path).expanduser().resolve()
    inspection = (
        sidecar
        if isinstance(sidecar, SyntenySidecarInspection)
        else inspect_synteny_sidecar(sidecar)
    )
    seal_path = database.parent / "dataset.seal.json"
    if not seal_path.is_file():
        return SyntenyDatasetIdentityInspection(
            state="seal_unavailable",
            compatible=True,
            sidecar_dataset_id=inspection.dataset_id,
            seal_path=seal_path,
            seal_dataset_id=None,
        )
    try:
        payload = json.loads(seal_path.read_text())
        seal_dataset_id = payload.get("dataset_id")
        if not isinstance(seal_dataset_id, str) or not seal_dataset_id:
            raise ValueError("dataset seal has no non-empty dataset_id")
    except Exception as exc:
        return SyntenyDatasetIdentityInspection(
            state="seal_unreadable",
            compatible=False,
            sidecar_dataset_id=inspection.dataset_id,
            seal_path=seal_path,
            seal_dataset_id=None,
            error=f"{type(exc).__name__}: {exc}",
        )
    if inspection.dataset_id is None:
        return SyntenyDatasetIdentityInspection(
            state="sidecar_dataset_id_missing",
            compatible=False,
            sidecar_dataset_id=None,
            seal_path=seal_path,
            seal_dataset_id=seal_dataset_id,
        )
    matches = inspection.dataset_id == seal_dataset_id
    return SyntenyDatasetIdentityInspection(
        state="match" if matches else "mismatch",
        compatible=matches,
        sidecar_dataset_id=inspection.dataset_id,
        seal_path=seal_path,
        seal_dataset_id=seal_dataset_id,
    )


def _fetch_dicts(
    connection: duckdb.DuckDBPyConnection,
    query: str,
    params: Sequence[Any] | None = None,
) -> list[dict[str, Any]]:
    cursor = connection.execute(query, list(params)) if params else connection.execute(query)
    columns = [str(column[0]) for column in cursor.description]
    return [dict(zip(columns, row, strict=True)) for row in cursor.fetchall()]


class SyntenyStore:
    """Connected read-only query interface for one normalized ELSA sidecar."""

    def __init__(
        self,
        path: str | Path,
        *,
        core_db_path: str | Path | None = None,
        allow_stale: bool = False,
    ):
        self.path = Path(path).expanduser().resolve()
        self.core_db_path = (
            Path(core_db_path).expanduser().resolve()
            if core_db_path is not None
            else None
        )
        self.allow_stale = allow_stale
        self._connection: duckdb.DuckDBPyConnection | None = None

    @property
    def connection(self) -> duckdb.DuckDBPyConnection:
        if self._connection is None:
            inspection = inspect_synteny_sidecar(self.path)
            if inspection.state != "available":
                raise SyntenySidecarError(
                    f"ELSA sidecar is {inspection.state}: "
                    f"{inspection.error or inspection.status or self.path}"
                )
            if self.core_db_path is not None:
                identity = inspect_synteny_dataset_identity(
                    self.core_db_path,
                    inspection,
                )
                if not identity.compatible and not self.allow_stale:
                    raise SyntenyDatasetMismatchError(
                        "ELSA sidecar dataset identity is "
                        f"{identity.state}: sidecar={identity.sidecar_dataset_id!r}, "
                        f"live={identity.seal_dataset_id!r}. "
                        "Refresh the ELSA run or opt into historical results "
                        "with allow_stale_synteny=True."
                    )
            self._connection = duckdb.connect(str(self.path), read_only=True)
        return self._connection

    def close(self) -> None:
        if self._connection is not None:
            self._connection.close()
            self._connection = None

    def __enter__(self) -> SyntenyStore:
        return self

    def __exit__(self, *_exc_info: object) -> None:
        self.close()

    def active_run_id(self) -> str:
        row = self.connection.execute(
            "SELECT run_id FROM current_elsa_run"
        ).fetchone()
        if row is None:
            raise SyntenySidecarError("ELSA sidecar has no active run")
        return str(row[0])

    def resolve_run_id(self, run_id: str | None) -> str:
        resolved = run_id or self.active_run_id()
        row = self.connection.execute(
            "SELECT status FROM elsa_runs WHERE run_id = ?",
            [resolved],
        ).fetchone()
        if row is None:
            raise KeyError(f"ELSA run {resolved!r} was not found")
        if row[0] != "ready":
            raise SyntenySidecarError(
                f"ELSA run {resolved!r} has status {row[0]!r}"
            )
        return resolved

    def runs(self) -> list[dict[str, Any]]:
        return _fetch_dicts(
            self.connection,
            """
            SELECT
                run_id, run_label, created_at, status, is_active, dataset_id,
                elsa_version, elsa_commit, mapping_version, gene_count,
                block_count, cluster_count, singleton_count, anchor_pair_count,
                locus_count, member_count, validation_json
            FROM elsa_runs
            ORDER BY is_active DESC, created_at DESC, run_id
            """,
        )

    def protein_memberships(
        self,
        protein_ids: Sequence[str],
        *,
        run_id: str | None = None,
        limit: int = 50,
        offset: int = 0,
    ) -> tuple[list[dict[str, Any]], int, str]:
        """Return exact protein↔cluster-locus memberships."""

        if limit < 1 or limit > 1_000:
            raise ValueError("limit must be between 1 and 1000")
        if offset < 0:
            raise ValueError("offset must be non-negative")
        targets = list(dict.fromkeys(str(value) for value in protein_ids if str(value)))
        if not targets:
            return [], 0, self.resolve_run_id(run_id)
        resolved_run = self.resolve_run_id(run_id)
        base = """
            FROM elsa_cluster_members AS members
            JOIN elsa_clusters AS clusters
              ON clusters.run_id = members.run_id
             AND clusters.cluster_key = members.cluster_key
            JOIN elsa_cluster_loci AS loci
              ON loci.run_id = members.run_id
             AND loci.locus_key = members.locus_key
            WHERE members.run_id = ?
              AND members.protein_id IN (
                  SELECT UNNEST(?::VARCHAR[])
              )
        """
        total = int(
            self.connection.execute(
                f"SELECT COUNT(*) {base}",
                [resolved_run, targets],
            ).fetchone()[0]
        )
        rows = _fetch_dicts(
            self.connection,
            f"""
            SELECT
                members.run_id,
                members.protein_id,
                members.member_role,
                members.cluster_key,
                clusters.source_cluster_id,
                clusters.cluster_kind,
                clusters.size AS block_count,
                clusters.genome_support,
                clusters.locus_count,
                clusters.member_count,
                clusters.anchor_member_count,
                members.locus_key,
                loci.genome_id,
                loci.contig_id,
                loci.start_position_index,
                loci.end_position_index,
                loci.start_bp,
                loci.end_bp,
                loci.block_support
            {base}
            ORDER BY
                members.protein_id,
                CASE members.member_role WHEN 'anchor' THEN 0 ELSE 1 END,
                clusters.genome_support DESC,
                clusters.size DESC,
                members.cluster_key,
                members.locus_key
            LIMIT ? OFFSET ?
            """,
            [resolved_run, targets, limit, offset],
        )
        return rows, total, resolved_run

    def anchor_blocks(
        self,
        protein_id: str,
        *,
        run_id: str | None = None,
        limit: int = 50,
    ) -> tuple[list[dict[str, Any]], int, str]:
        """Return blocks where ``protein_id`` is an exact anchor."""

        if limit < 1 or limit > 1_000:
            raise ValueError("limit must be between 1 and 1000")
        resolved_run = self.resolve_run_id(run_id)
        query = """
            WITH hits AS (
                SELECT
                    pairs.run_id, pairs.block_id, pairs.cluster_key,
                    pairs.pair_order, 'query' AS anchor_side,
                    pairs.query_protein_id AS protein_id,
                    pairs.target_protein_id AS partner_protein_id
                FROM elsa_anchor_pairs AS pairs
                WHERE pairs.run_id = ? AND pairs.query_protein_id = ?
                UNION ALL
                SELECT
                    pairs.run_id, pairs.block_id, pairs.cluster_key,
                    pairs.pair_order, 'target' AS anchor_side,
                    pairs.target_protein_id AS protein_id,
                    pairs.query_protein_id AS partner_protein_id
                FROM elsa_anchor_pairs AS pairs
                WHERE pairs.run_id = ? AND pairs.target_protein_id = ?
            )
        """
        total = int(
            self.connection.execute(
                query + " SELECT COUNT(*) FROM hits",
                [resolved_run, protein_id, resolved_run, protein_id],
            ).fetchone()[0]
        )
        rows = _fetch_dicts(
            self.connection,
            query
            + """
            SELECT
                hits.run_id, hits.block_id, hits.cluster_key, hits.pair_order,
                hits.anchor_side, hits.protein_id, hits.partner_protein_id,
                blocks.source_cluster_id, blocks.query_genome,
                blocks.target_genome, blocks.query_contig,
                blocks.target_contig, blocks.orientation, blocks.n_anchors,
                blocks.chain_score
            FROM hits
            JOIN elsa_blocks AS blocks
              ON blocks.run_id = hits.run_id AND blocks.block_id = hits.block_id
            ORDER BY blocks.n_anchors DESC, blocks.chain_score DESC, hits.block_id
            LIMIT ?
            """,
            [resolved_run, protein_id, resolved_run, protein_id, limit],
        )
        return rows, total, resolved_run

    def resolve_cluster_key(
        self,
        cluster: str | int,
        *,
        run_id: str,
    ) -> str:
        """Resolve a canonical key or a run-scoped numeric source ID."""

        if isinstance(cluster, int) or str(cluster).isdigit():
            source_id = int(cluster)
            rows = self.connection.execute(
                """
                SELECT cluster_key
                FROM elsa_clusters
                WHERE run_id = ? AND source_cluster_id = ?
                """,
                [run_id, source_id],
            ).fetchall()
        else:
            rows = self.connection.execute(
                """
                SELECT cluster_key
                FROM elsa_clusters
                WHERE run_id = ? AND cluster_key = ?
                """,
                [run_id, str(cluster)],
            ).fetchall()
        if not rows:
            raise KeyError(
                f"ELSA cluster {cluster!r} was not found in run {run_id!r}"
            )
        if len(rows) > 1:
            raise SyntenySidecarError(
                f"ELSA cluster reference {cluster!r} is ambiguous in run {run_id!r}"
            )
        return str(rows[0][0])

    def cluster(
        self,
        cluster: str | int,
        *,
        run_id: str | None = None,
        member_limit: int = 100,
        member_offset: int = 0,
    ) -> dict[str, Any]:
        """Return one cluster summary with bounded loci and member records."""

        if member_limit < 1 or member_limit > 5_000:
            raise ValueError("member_limit must be between 1 and 5000")
        if member_offset < 0:
            raise ValueError("member_offset must be non-negative")
        resolved_run = self.resolve_run_id(run_id)
        cluster_key = self.resolve_cluster_key(cluster, run_id=resolved_run)
        summary_rows = _fetch_dicts(
            self.connection,
            """
            SELECT *
            FROM elsa_clusters
            WHERE run_id = ? AND cluster_key = ?
            """,
            [resolved_run, cluster_key],
        )
        loci = _fetch_dicts(
            self.connection,
            """
            SELECT *
            FROM elsa_cluster_loci
            WHERE run_id = ? AND cluster_key = ?
            ORDER BY genome_id, contig_id, start_position_index, locus_key
            """,
            [resolved_run, cluster_key],
        )
        members = _fetch_dicts(
            self.connection,
            """
            SELECT
                protein_id, member_role, locus_key, genome_id, contig_id,
                position_index, locus_order, start_bp, end_bp, strand
            FROM elsa_cluster_members
            WHERE run_id = ? AND cluster_key = ?
            ORDER BY
                genome_id, contig_id, position_index, protein_id
            LIMIT ? OFFSET ?
            """,
            [resolved_run, cluster_key, member_limit, member_offset],
        )
        payload = summary_rows[0]
        payload["loci"] = loci
        payload["members"] = members
        payload["member_offset"] = member_offset
        payload["members_truncated"] = member_offset + len(members) < payload["member_count"]
        validation = self.connection.execute(
            "SELECT validation_json FROM elsa_runs WHERE run_id = ?",
            [resolved_run],
        ).fetchone()[0]
        payload["run_validation"] = json.loads(validation) if validation else {}
        return payload


__all__ = [
    "DEFAULT_SIDECAR_NAME",
    "EXPECTED_SCHEMA_VERSION",
    "SyntenyDatasetIdentityInspection",
    "SyntenyDatasetMismatchError",
    "SyntenySidecarError",
    "SyntenySidecarInspection",
    "SyntenyStore",
    "discover_synteny_sidecar",
    "inspect_synteny_dataset_identity",
    "inspect_synteny_sidecar",
]
