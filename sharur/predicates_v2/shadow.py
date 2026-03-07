"""
Shadow diff: run V1 and V2 on the same data, compare outputs.

Produces both a human-readable report and machine-readable JSONL
showing per-protein differences between V1 and V2 predicate systems.
"""

from __future__ import annotations

import json
from collections import Counter, defaultdict
from typing import TYPE_CHECKING, Any, Optional

from sharur.predicates.generator import (
    AnnotationRecord,
    PredicateGenerator,
    ProteinRecord,
)
from sharur.predicates_v2.aggregator import aggregate_atoms
from sharur.predicates_v2.compat import semantic_state_to_predicates
from sharur.predicates_v2.composites import evaluate_composites, load_composites
from sharur.predicates_v2.generator import AtomGenerator
from sharur.predicates_v2.model import SemanticAtom

if TYPE_CHECKING:
    from sharur.storage.duckdb_store import DuckDBStore


def shadow_diff(
    proteins: list[ProteinRecord],
    annotations_by_protein: dict[str, list[AnnotationRecord]],
) -> dict[str, Any]:
    """Run V1 and V2 on the same data and produce a diff.

    Args:
        proteins: List of protein records.
        annotations_by_protein: Mapping protein_id -> annotations.

    Returns:
        Dict with keys:
            - per_protein: list of per-protein diffs
            - summary: aggregate statistics
    """
    v1_gen = PredicateGenerator(
        expand_hierarchy=True,
        include_direct_access=True,
        predict_topology=False,  # Skip topology for consistency
    )
    v2_gen = AtomGenerator(
        expand_hierarchy=True,
        predict_topology=False,
    )
    composites = load_composites()

    per_protein_diffs = []
    total_v1_preds = Counter()
    total_v2_preds = Counter()
    match_count = 0
    diff_count = 0
    all_v2_atoms: list[SemanticAtom] = []

    for protein in proteins:
        pid = protein.protein_id
        annotations = annotations_by_protein.get(pid, [])

        # V1: flat predicates
        v1_preds = set(v1_gen.generate_for_protein(protein, annotations))

        # V2: atoms -> state -> flat predicates (via compat layer)
        v2_atoms = v2_gen.generate_atoms(protein, annotations)
        all_v2_atoms.extend(v2_atoms)
        v2_state = aggregate_atoms(pid, v2_atoms)
        composite_hits = evaluate_composites(v2_atoms, composites, v2_state.topology)
        v2_state.composite_predicates = composite_hits
        v2_preds_flat = set(semantic_state_to_predicates(v2_state, atoms=v2_atoms))

        # Compute diff
        # Filter out direct-access predicates from V1 (pfam:*, kegg:*, etc.)
        # since V2 doesn't produce these
        v1_semantic = {p for p in v1_preds if ":" not in p}

        only_v1 = v1_semantic - v2_preds_flat
        only_v2 = v2_preds_flat - v1_semantic
        common = v1_semantic & v2_preds_flat

        # Count totals
        for p in v1_semantic:
            total_v1_preds[p] += 1
        for p in v2_preds_flat:
            total_v2_preds[p] += 1

        if only_v1 or only_v2:
            diff_count += 1
            per_protein_diffs.append({
                "protein_id": pid,
                "v1_only": sorted(only_v1),
                "v2_only": sorted(only_v2),
                "common_count": len(common),
                "v2_composites": composite_hits,
                "v2_unresolved_count": len(v2_state.unresolved_accessions),
            })
        else:
            match_count += 1

    total_proteins = len(proteins)
    match_rate = (match_count / total_proteins * 100) if total_proteins > 0 else 0

    # Top V2 atoms by frequency with facet distribution
    atom_facet_counter: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    atom_relation_counter: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for atom in all_v2_atoms:
        atom_facet_counter[atom.atom_id][atom.facet.value] += 1
        atom_relation_counter[atom.atom_id][atom.relation.value] += 1

    top_atoms = sorted(
        atom_facet_counter.items(),
        key=lambda x: sum(x[1].values()),
        reverse=True,
    )[:20]

    # Count unresolved
    unresolved_count = sum(
        1 for a in all_v2_atoms if a.relation.value == "unresolved"
    )
    unique_unresolved = len({
        a.source_accession for a in all_v2_atoms if a.relation.value == "unresolved"
    })

    summary = {
        "total_proteins": total_proteins,
        "match_count": match_count,
        "diff_count": diff_count,
        "match_rate_pct": round(match_rate, 1),
        "v1_unique_predicates": len(total_v1_preds),
        "v2_unique_predicates": len(total_v2_preds),
        "v2_total_atoms": len(all_v2_atoms),
        "v2_unresolved_atoms": unresolved_count,
        "v2_unique_unresolved": unique_unresolved,
        "top_20_atoms": [
            {
                "atom_id": atom_id,
                "total_count": sum(facets.values()),
                "facets": dict(facets),
                "relations": dict(atom_relation_counter[atom_id]),
            }
            for atom_id, facets in top_atoms
        ],
    }

    return {
        "per_protein": per_protein_diffs,
        "summary": summary,
    }


def format_shadow_diff_report(diff: dict[str, Any]) -> str:
    """Format shadow diff as human-readable Markdown report.

    Args:
        diff: Output from shadow_diff().

    Returns:
        Markdown-formatted report string.
    """
    summary = diff["summary"]
    per_protein = diff["per_protein"]

    lines = []
    lines.append("# Predicate V2 Shadow Diff Report")
    lines.append("")
    lines.append("## Summary")
    lines.append("")
    lines.append(f"- Total proteins: {summary['total_proteins']}")
    lines.append(f"- Exact match: {summary['match_count']} ({summary['match_rate_pct']}%)")
    lines.append(f"- Differences: {summary['diff_count']}")
    lines.append(f"- V1 unique predicates: {summary['v1_unique_predicates']}")
    lines.append(f"- V2 unique predicates: {summary['v2_unique_predicates']}")
    lines.append(f"- V2 total atoms: {summary['v2_total_atoms']}")
    lines.append(f"- V2 unresolved atoms: {summary['v2_unresolved_atoms']}")
    lines.append(f"- V2 unique unresolved accessions: {summary['v2_unique_unresolved']}")
    lines.append("")

    lines.append("## Top 20 Atoms by Frequency")
    lines.append("")
    lines.append("| Atom | Count | Facet | Relations |")
    lines.append("|------|-------|-------|-----------|")
    for atom_info in summary["top_20_atoms"]:
        facet_str = ", ".join(f"{k}:{v}" for k, v in atom_info["facets"].items())
        rel_str = ", ".join(f"{k}:{v}" for k, v in atom_info["relations"].items())
        lines.append(
            f"| {atom_info['atom_id']} | {atom_info['total_count']} "
            f"| {facet_str} | {rel_str} |"
        )
    lines.append("")

    if per_protein:
        lines.append("## Per-Protein Differences (first 50)")
        lines.append("")
        for diff_entry in per_protein[:50]:
            pid = diff_entry["protein_id"]
            lines.append(f"### {pid}")
            if diff_entry["v1_only"]:
                lines.append(f"  V1 only: {', '.join(diff_entry['v1_only'])}")
            if diff_entry["v2_only"]:
                lines.append(f"  V2 only: {', '.join(diff_entry['v2_only'])}")
            if diff_entry["v2_composites"]:
                lines.append(f"  Composites: {', '.join(diff_entry['v2_composites'])}")
            lines.append("")

    return "\n".join(lines)


def format_shadow_diff_jsonl(diff: dict[str, Any]) -> str:
    """Format per-protein diffs as JSONL.

    Args:
        diff: Output from shadow_diff().

    Returns:
        JSONL string (one JSON object per line).
    """
    lines = []
    for entry in diff["per_protein"]:
        lines.append(json.dumps(entry))
    return "\n".join(lines) + "\n" if lines else ""


def run_shadow_diff_on_store(
    store: "DuckDBStore",
    protein_ids: Optional[list[str]] = None,
    output_report: Optional[str] = None,
    output_jsonl: Optional[str] = None,
) -> dict[str, Any]:
    """Run shadow diff on a DuckDB store.

    Convenience function that reads data from the store and runs
    the shadow diff.

    Args:
        store: DuckDB store.
        protein_ids: Optional subset.
        output_report: Path to write Markdown report.
        output_jsonl: Path to write JSONL diffs.

    Returns:
        Shadow diff results dict.
    """
    from sharur.predicates.generator import _get_contig_gc_stats

    gc_stats = _get_contig_gc_stats(store)

    # Get proteins
    if protein_ids:
        placeholders = ",".join(["?"] * len(protein_ids))
        rows = store.execute(
            f"SELECT protein_id, sequence_length, gc_content, contig_id "
            f"FROM proteins WHERE protein_id IN ({placeholders})",
            protein_ids,
        )
    else:
        rows = store.execute(
            "SELECT protein_id, sequence_length, gc_content, contig_id FROM proteins"
        )

    proteins = []
    for row in rows:
        gc_mean, gc_std = gc_stats.get(row[3], (None, None)) if row[3] else (None, None)
        proteins.append(ProteinRecord(
            protein_id=row[0],
            sequence_length=row[1],
            gc_content=row[2],
            contig_gc_mean=gc_mean,
            contig_gc_std=gc_std,
        ))

    # Get annotations
    if protein_ids:
        placeholders = ",".join(["?"] * len(protein_ids))
        ann_rows = store.execute(
            f"SELECT protein_id, source, accession, name, description, evalue, score "
            f"FROM annotations WHERE protein_id IN ({placeholders})",
            protein_ids,
        )
    else:
        ann_rows = store.execute(
            "SELECT protein_id, source, accession, name, description, evalue, score "
            "FROM annotations"
        )

    annotations_by_protein: dict[str, list[AnnotationRecord]] = {}
    for row in ann_rows:
        pid = row[0]
        if pid not in annotations_by_protein:
            annotations_by_protein[pid] = []
        annotations_by_protein[pid].append(
            AnnotationRecord(
                source=row[1],
                accession=row[2],
                name=row[3],
                description=row[4],
                evalue=row[5],
                score=row[6],
            )
        )

    # Run diff
    result = shadow_diff(proteins, annotations_by_protein)

    # Write outputs
    if output_report:
        with open(output_report, "w") as f:
            f.write(format_shadow_diff_report(result))

    if output_jsonl:
        with open(output_jsonl, "w") as f:
            f.write(format_shadow_diff_jsonl(result))

    return result


__all__ = [
    "format_shadow_diff_jsonl",
    "format_shadow_diff_report",
    "run_shadow_diff_on_store",
    "shadow_diff",
]
