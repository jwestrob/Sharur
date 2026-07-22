"""Token-bounded, exact ELSA synteny operators."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

from sharur.operators.base import OperatorContext, SharurResult
from sharur.synteny import SyntenyStore


if TYPE_CHECKING:
    from collections.abc import Sequence

    from sharur.storage.duckdb_store import DuckDBStore


def _membership_markdown(
    rows: list[dict[str, Any]],
    *,
    protein_ids: Sequence[str],
    run_id: str,
    total: int,
) -> str:
    lines = [
        "# ELSA exact synteny membership",
        "",
        f"- Run: `{run_id}`",
        f"- Proteins queried: {len(protein_ids)}",
        f"- Memberships: {total}",
    ]
    if not rows:
        return "\n".join(lines)
    lines.extend(
        [
            "",
            "| Protein | Role | Cluster | Blocks | Genomes | Locus |",
            "|---|---|---|---:|---:|---|",
        ]
    )
    for row in rows:
        source = row["source_cluster_id"]
        cluster = (
            f"{row['cluster_key']} (source {source})"
            if source is not None
            else row["cluster_key"]
        )
        lines.append(
            f"| `{row['protein_id']}` | {row['member_role']} | "
            f"`{cluster}` | {row['block_count']} | "
            f"{row['genome_support']} | `{row['locus_key']}` |"
        )
    if len(rows) < total:
        lines.extend(
            [
                "",
                f"Returned {len(rows)} of {total} memberships; use offset for the next page.",
            ]
        )
    return "\n".join(lines)


def synteny_for_proteins(
    core_store: DuckDBStore,
    sidecar_path: str | Path,
    protein_ids: Sequence[str],
    *,
    run_id: str | None = None,
    limit: int = 50,
    offset: int = 0,
    allow_stale: bool = False,
) -> SharurResult:
    """Return exact cluster-locus memberships for one or more proteins."""

    targets = list(dict.fromkeys(str(value) for value in protein_ids if str(value)))
    params = {
        "protein_ids": targets,
        "run_id": run_id,
        "limit": limit,
        "offset": offset,
        "allow_stale": allow_stale,
        "sidecar": str(Path(sidecar_path).expanduser().resolve()),
    }
    with OperatorContext("synteny_for_proteins", params, store=core_store) as context:
        with SyntenyStore(
            sidecar_path,
            core_db_path=core_store.db_path,
            allow_stale=allow_stale,
        ) as synteny:
            rows, total, resolved_run = synteny.protein_memberships(
                targets,
                run_id=run_id,
                limit=limit,
                offset=offset,
            )
        return context.make_result(
            data=_membership_markdown(
                rows,
                protein_ids=targets,
                run_id=resolved_run,
                total=total,
            ),
            rows=len(rows),
            total_rows=total,
            truncated=len(rows) < total,
            raw=rows,
            index_used="elsa_cluster_members(run_id, protein_id)",
            dataset_version=f"{context._trace_versions()[0]}|elsa_run={resolved_run}",
        )


def synteny_anchor_blocks(
    core_store: DuckDBStore,
    sidecar_path: str | Path,
    protein_id: str,
    *,
    run_id: str | None = None,
    limit: int = 50,
    allow_stale: bool = False,
) -> SharurResult:
    """Return exact orientation-resolved anchor pairs for one protein."""

    params = {
        "protein_id": protein_id,
        "run_id": run_id,
        "limit": limit,
        "allow_stale": allow_stale,
        "sidecar": str(Path(sidecar_path).expanduser().resolve()),
    }
    with OperatorContext("synteny_anchor_blocks", params, store=core_store) as context:
        with SyntenyStore(
            sidecar_path,
            core_db_path=core_store.db_path,
            allow_stale=allow_stale,
        ) as synteny:
            rows, total, resolved_run = synteny.anchor_blocks(
                protein_id,
                run_id=run_id,
                limit=limit,
            )
        lines = [
            f"# ELSA exact anchor blocks: {protein_id}",
            "",
            f"- Run: `{resolved_run}`",
            f"- Anchor pairs: {total}",
        ]
        if rows:
            lines.extend(
                [
                    "",
                    "| Block | Cluster | Side | Partner | Orientation | Anchors | Score |",
                    "|---:|---|---|---|---:|---:|---:|",
                ]
            )
            for row in rows:
                lines.append(
                    f"| {row['block_id']} | `{row['cluster_key']}` | "
                    f"{row['anchor_side']} | `{row['partner_protein_id']}` | "
                    f"{row['orientation']} | {row['n_anchors']} | "
                    f"{row['chain_score']:.4g} |"
                )
        return context.make_result(
            data="\n".join(lines),
            rows=len(rows),
            total_rows=total,
            truncated=len(rows) < total,
            raw=rows,
            index_used="elsa_anchor_pairs(run_id, protein_id)",
            dataset_version=f"{context._trace_versions()[0]}|elsa_run={resolved_run}",
        )


def get_synteny_cluster(
    core_store: DuckDBStore,
    sidecar_path: str | Path,
    cluster: str | int,
    *,
    run_id: str | None = None,
    member_limit: int = 100,
    member_offset: int = 0,
    allow_stale: bool = False,
) -> SharurResult:
    """Return one run-scoped cluster with bounded exact members."""

    params = {
        "cluster": cluster,
        "run_id": run_id,
        "member_limit": member_limit,
        "member_offset": member_offset,
        "allow_stale": allow_stale,
        "sidecar": str(Path(sidecar_path).expanduser().resolve()),
    }
    with OperatorContext("get_synteny_cluster", params, store=core_store) as context:
        with SyntenyStore(
            sidecar_path,
            core_db_path=core_store.db_path,
            allow_stale=allow_stale,
        ) as synteny:
            payload = synteny.cluster(
                cluster,
                run_id=run_id,
                member_limit=member_limit,
                member_offset=member_offset,
            )
        lines = [
            f"# ELSA cluster: {payload['cluster_key']}",
            "",
            f"- Run: `{payload['run_id']}`",
            f"- Source cluster ID: `{payload['source_cluster_id']}`",
            f"- Kind: `{payload['cluster_kind']}`",
            f"- Blocks: {payload['size']}",
            f"- Genome support: {payload['genome_support']}",
            f"- Loci: {payload['locus_count']}",
            f"- Members: {payload['member_count']} "
            f"({payload['anchor_member_count']} anchors)",
            "",
            "| Genome | Contig | Rank interval | Genes | Block support |",
            "|---|---|---:|---:|---:|",
        ]
        for locus in payload["loci"]:
            lines.append(
                f"| `{locus['genome_id']}` | `{locus['contig_id']}` | "
                f"{locus['start_position_index']}-{locus['end_position_index']} | "
                f"{locus['n_genes']} | {locus['block_support']} |"
            )
        if payload["members_truncated"]:
            lines.extend(
                [
                    "",
                    f"Members are paginated from offset {member_offset}; "
                    f"returned {len(payload['members'])}.",
                ]
            )
        return context.make_result(
            data="\n".join(lines),
            rows=len(payload["members"]),
            total_rows=int(payload["member_count"]),
            truncated=bool(payload["members_truncated"]),
            raw=payload,
            index_used="elsa_cluster_members(run_id, cluster_key)",
            dataset_version=(
                f"{context._trace_versions()[0]}|elsa_run={payload['run_id']}"
            ),
        )


__all__ = [
    "get_synteny_cluster",
    "synteny_anchor_blocks",
    "synteny_for_proteins",
]
