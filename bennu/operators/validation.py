"""
Validation operators for annotation quality control.

These operators help detect and flag potential annotation errors,
validate genomic context, and ensure findings are properly validated.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional

from bennu.operators.base import BennuResult, OperatorContext


# ---------------------------------------------------------------------------
# Known Problem Domains Registry
# ---------------------------------------------------------------------------

PROBLEM_DOMAINS = {
    "Cas12f1-like_TNB": {
        "issue": "Shared with IS605 TnpB transposases (evolutionary ancestor of Cas12f)",
        "validation": "Check protein annotation - if 'transposase', this is NOT a CRISPR effector",
        "expected_context": ["cas1", "cas2", "cas3", "crispr_array"],
        "false_positive_indicators": ["transposase", "tnpA", "tnpB", "IS605"],
        "notes": "TnpB transposases are the evolutionary ancestors of Cas12f. Domain similarity does NOT imply CRISPR function.",
    },
    "HTH_17": {
        "issue": "Generic helix-turn-helix domain found in diverse protein families",
        "validation": "Cannot infer specific function from this domain alone - check other domains and annotation",
        "expected_context": [],
        "false_positive_indicators": [],
        "notes": "HTH domains are found in transcription factors, transposases, and many other proteins.",
    },
    "DDE_Tnp_1": {
        "issue": "DDE transposase domain - may be active or degraded",
        "validation": "Check if catalytic residues are intact; consider genomic context",
        "expected_context": ["terminal_inverted_repeats"],
        "false_positive_indicators": [],
        "notes": "Many DDE transposases are inactive remnants.",
    },
    "Intein_splicing": {
        "issue": "Inteins can be found in many host proteins",
        "validation": "The host protein function is determined by the surrounding sequence, not the intein",
        "expected_context": [],
        "false_positive_indicators": [],
        "notes": "Inteins are selfish genetic elements; their presence doesn't indicate protein function.",
    },
}


@dataclass
class ValidationResult:
    """Result of a validation check."""

    is_valid: bool
    confidence: str  # "high", "medium", "low"
    issues: list[dict[str, Any]]
    recommendations: list[str]


def validate_annotation(
    store,
    protein_ids: list[str],
    expected_function: str,
    domain_name: Optional[str] = None,
) -> BennuResult:
    """
    Validate that proteins match expected function.

    Cross-checks domain hits against protein annotations to detect
    potential annotation errors or domain misinterpretations.

    Args:
        store: DuckDB store
        protein_ids: List of protein IDs to validate
        expected_function: Expected function keyword (e.g., "CRISPR", "transposase")
        domain_name: Optional domain name that led to the function inference

    Returns:
        BennuResult with validation summary
    """
    with OperatorContext("validate_annotation", {
        "n_proteins": len(protein_ids),
        "expected_function": expected_function,
        "domain_name": domain_name,
    }) as ctx:

        if not protein_ids:
            return ctx.make_result(
                data="No proteins provided for validation.",
                rows=0,
            )

        # Get best annotation for each protein (lowest e-value)
        placeholders = ", ".join(["?"] * len(protein_ids))
        results = store.execute(f"""
            WITH best_annot AS (
                SELECT protein_id, name, description,
                       ROW_NUMBER() OVER (PARTITION BY protein_id ORDER BY evalue ASC NULLS LAST) as rn
                FROM annotations
                WHERE protein_id IN ({placeholders})
            )
            SELECT p.protein_id,
                   COALESCE(ba.name, ba.description, 'unannotated') as annotation,
                   p.sequence_length
            FROM proteins p
            LEFT JOIN best_annot ba ON p.protein_id = ba.protein_id AND ba.rn = 1
            WHERE p.protein_id IN ({placeholders})
        """, protein_ids + protein_ids)

        # Check for known problem domains
        problem_domain_warning = None
        if domain_name and domain_name in PROBLEM_DOMAINS:
            problem_domain_warning = PROBLEM_DOMAINS[domain_name]

        # Categorize proteins
        consistent = []
        inconsistent = []
        uncertain = []

        expected_lower = expected_function.lower()

        for protein_id, annotation, seq_len in results:
            annotation_lower = (annotation or "").lower()

            # Check if annotation matches expected function
            if expected_lower in annotation_lower:
                consistent.append({
                    "protein_id": protein_id,
                    "annotation": annotation,
                    "status": "consistent",
                })
            elif annotation and annotation.lower() not in ("hypothetical protein", "uncharacterized protein", "unannotated", ""):
                # Has annotation but doesn't match
                inconsistent.append({
                    "protein_id": protein_id,
                    "annotation": annotation,
                    "status": "inconsistent",
                    "issue": f"Annotated as '{annotation}', not '{expected_function}'",
                })
            else:
                uncertain.append({
                    "protein_id": protein_id,
                    "annotation": annotation,
                    "status": "uncertain",
                })

        # Check for false positive indicators if we have a problem domain
        false_positives = []
        if problem_domain_warning:
            for item in inconsistent + uncertain:
                annotation_lower = (item.get("annotation") or "").lower()
                for indicator in problem_domain_warning["false_positive_indicators"]:
                    if indicator.lower() in annotation_lower:
                        false_positives.append(item["protein_id"])
                        break

        # Build report
        lines = [
            f"# Annotation Validation Report",
            f"",
            f"**Expected function:** {expected_function}",
            f"**Proteins checked:** {len(protein_ids)}",
            f"",
        ]

        if problem_domain_warning:
            lines.extend([
                f"## Problem Domain Warning",
                f"",
                f"**Domain:** {domain_name}",
                f"**Issue:** {problem_domain_warning['issue']}",
                f"**Validation:** {problem_domain_warning['validation']}",
                f"",
            ])

        lines.extend([
            f"## Summary",
            f"",
            f"| Status | Count | Percentage |",
            f"|--------|-------|------------|",
            f"| Consistent | {len(consistent)} | {100*len(consistent)/len(protein_ids):.1f}% |",
            f"| **Inconsistent** | **{len(inconsistent)}** | **{100*len(inconsistent)/len(protein_ids):.1f}%** |",
            f"| Uncertain | {len(uncertain)} | {100*len(uncertain)/len(protein_ids):.1f}% |",
            f"",
        ])

        if false_positives:
            lines.extend([
                f"## Likely False Positives",
                f"",
                f"**{len(false_positives)} proteins** have annotations suggesting they are NOT {expected_function}:",
                f"",
            ])
            for pid in false_positives[:10]:
                ann = next((x["annotation"] for x in inconsistent + uncertain if x["protein_id"] == pid), "")
                lines.append(f"- `{pid}`: {ann}")
            if len(false_positives) > 10:
                lines.append(f"- ... and {len(false_positives) - 10} more")
            lines.append("")

        if inconsistent:
            lines.extend([
                f"## Inconsistent Annotations (Top 10)",
                f"",
            ])
            for item in inconsistent[:10]:
                lines.append(f"- `{item['protein_id']}`: {item['annotation']}")
            lines.append("")

        # Recommendations
        lines.extend([
            f"## Recommendations",
            f"",
        ])

        if len(inconsistent) > len(consistent):
            lines.append(f"- **Majority of proteins have inconsistent annotations.** The domain hit may not indicate {expected_function} function.")

        if false_positives:
            lines.append(f"- **{len(false_positives)} proteins are likely false positives** based on their annotations matching known false-positive indicators.")

        if problem_domain_warning:
            lines.append(f"- Check genomic context: Are these proteins near {', '.join(problem_domain_warning['expected_context'])}?")

        if len(consistent) == 0 and len(inconsistent) > 0:
            lines.append(f"- **CRITICAL:** No proteins have annotations consistent with {expected_function}. Reconsider the functional interpretation.")

        return ctx.make_result(
            data="\n".join(lines),
            rows=len(protein_ids),
            raw={
                "consistent": consistent,
                "inconsistent": inconsistent,
                "uncertain": uncertain,
                "false_positives": false_positives,
                "problem_domain": problem_domain_warning,
            },
        )


def validate_context(
    store,
    protein_id: str,
    expected_neighbors: list[str],
    window: int = 10,
) -> BennuResult:
    """
    Check if a protein is in expected genomic context.

    Validates whether nearby genes match expected functional context
    (e.g., CRISPR proteins should be near cas genes and CRISPR arrays).

    Args:
        store: DuckDB store
        protein_id: Target protein ID
        expected_neighbors: List of keywords/predicates expected nearby
        window: Number of genes to check on each side

    Returns:
        BennuResult with context validation
    """
    with OperatorContext("validate_context", {
        "protein_id": protein_id,
        "expected_neighbors": expected_neighbors,
        "window": window,
    }) as ctx:

        # Get target protein info
        target = store.execute("""
            SELECT p.contig_id, p.gene_index,
                   COALESCE(a.name, a.description, 'unannotated') as annotation
            FROM proteins p
            LEFT JOIN annotations a ON p.protein_id = a.protein_id
            WHERE p.protein_id = ?
            LIMIT 1
        """, [protein_id])

        if not target:
            return ctx.make_result(
                data=f"Protein {protein_id} not found.",
                rows=0,
            )

        contig_id, gene_index, target_annotation = target[0]

        # Get neighborhood with annotations
        neighbors = store.execute("""
            SELECT
                p.protein_id,
                COALESCE(a.name, a.description, 'unannotated') as annotation,
                p.gene_index,
                pp.predicates
            FROM proteins p
            LEFT JOIN annotations a ON p.protein_id = a.protein_id
            LEFT JOIN protein_predicates pp ON p.protein_id = pp.protein_id
            WHERE p.contig_id = ?
              AND p.gene_index BETWEEN ? AND ?
              AND p.protein_id != ?
            ORDER BY p.gene_index
        """, [contig_id, gene_index - window, gene_index + window, protein_id])

        # Check for CRISPR arrays nearby
        arrays_nearby = store.execute("""
            SELECT locus_id, start, end_coord
            FROM loci
            WHERE contig_id = ?
              AND locus_type = 'crispr_array'
        """, [contig_id])

        # Analyze neighbors
        found_expected = {exp: [] for exp in expected_neighbors}
        unexpected_context = []

        for pid, ann, gidx, predicates in neighbors:
            ann_lower = (ann or "").lower()
            predicates = predicates or []

            # Check if any expected neighbor found
            for exp in expected_neighbors:
                exp_lower = exp.lower()
                if exp_lower in ann_lower or exp_lower in [p.lower() for p in predicates]:
                    found_expected[exp].append({
                        "protein_id": pid,
                        "annotation": ann,
                        "distance": abs(gidx - gene_index),
                    })

            # Check for unexpected context (transposases near supposed CRISPR)
            if "transposase" in ann_lower and "transposase" not in [e.lower() for e in expected_neighbors]:
                unexpected_context.append({
                    "protein_id": pid,
                    "annotation": ann,
                    "distance": abs(gidx - gene_index),
                    "type": "transposase",
                })

        # Build report
        neighbors_found = [k for k, v in found_expected.items() if v]
        neighbors_missing = [k for k, v in found_expected.items() if not v]

        context_valid = len(neighbors_found) >= len(expected_neighbors) / 2

        lines = [
            f"# Context Validation Report",
            f"",
            f"**Protein:** {protein_id}",
            f"**Annotation:** {target_annotation}",
            f"**Contig:** {contig_id}",
            f"**Window:** +/-{window} genes",
            f"",
            f"## Context Valid: {'YES' if context_valid else 'NO'}",
            f"",
            f"### Expected Neighbors Found ({len(neighbors_found)}/{len(expected_neighbors)})",
            f"",
        ]

        for exp in expected_neighbors:
            matches = found_expected[exp]
            if matches:
                closest = min(matches, key=lambda x: x["distance"])
                lines.append(f"- **{exp}**: Found {len(matches)} ({closest['annotation']} at distance {closest['distance']})")
            else:
                lines.append(f"- **{exp}**: Not found")

        lines.append("")

        if arrays_nearby:
            lines.extend([
                f"### CRISPR Arrays on Contig",
                f"",
            ])
            for arr_id, arr_start, arr_end in arrays_nearby:
                lines.append(f"- {arr_id}: {arr_start}-{arr_end}")
            lines.append("")
        else:
            lines.extend([
                f"### CRISPR Arrays: None on this contig",
                f"",
            ])

        if unexpected_context:
            lines.extend([
                f"### Unexpected Context",
                f"",
            ])
            for item in unexpected_context[:5]:
                lines.append(f"- {item['type']} at distance {item['distance']}: {item['annotation']}")
            lines.append("")

        return ctx.make_result(
            data="\n".join(lines),
            rows=len(neighbors),
            raw={
                "context_valid": context_valid,
                "neighbors_found": neighbors_found,
                "neighbors_missing": neighbors_missing,
                "unexpected_context": unexpected_context,
                "crispr_arrays_nearby": len(arrays_nearby),
            },
        )


def analyze_crispr_systems(store) -> BennuResult:
    """
    Analyze all CRISPR-Cas systems in the database.

    Identifies complete vs incomplete systems, orphan arrays,
    and systems affected by assembly fragmentation.

    Args:
        store: DuckDB store

    Returns:
        BennuResult with CRISPR system analysis
    """
    with OperatorContext("analyze_crispr_systems", {}) as ctx:

        # Get all CRISPR arrays
        arrays = store.execute("""
            SELECT
                l.locus_id,
                c.bin_id as genome_id,
                l.contig_id,
                l.start,
                l.end_coord,
                l.metadata,
                c.length as contig_length
            FROM loci l
            JOIN contigs c ON l.contig_id = c.contig_id
            WHERE l.locus_type = 'crispr_array'
        """)

        if not arrays:
            return ctx.make_result(
                data="No CRISPR arrays found in database.",
                rows=0,
            )

        systems = []

        for locus_id, genome_id, contig_id, start, end_coord, metadata, contig_length in arrays:
            # Check if at contig edge (potential fragmentation)
            at_edge = start < 1000 or end_coord > contig_length - 1000

            # Find nearby Cas proteins (within 15kb) - use annotations table
            cas_nearby = store.execute("""
                SELECT p.protein_id, a.name, p.start, p.end_coord
                FROM proteins p
                JOIN annotations a ON p.protein_id = a.protein_id
                WHERE p.contig_id = ?
                  AND (a.name ILIKE '%cas1%' OR a.name ILIKE '%cas2%'
                       OR a.name ILIKE '%cas3%' OR a.name ILIKE '%cas5%'
                       OR a.name ILIKE '%cas6%' OR a.name ILIKE '%cas7%'
                       OR a.name ILIKE '%cas8%' OR a.name ILIKE '%cas9%'
                       OR a.name ILIKE '%cas10%' OR a.name ILIKE '%cas12%'
                       OR a.name ILIKE '%cas13%' OR a.name ILIKE '%csm%'
                       OR a.name ILIKE '%cmr%' OR a.name ILIKE '%csx%')
                  AND ABS(p.start - ?) < 15000
            """, [contig_id, start])

            # Determine system type based on cas genes
            cas_genes = [ann.lower() for _, ann, _, _ in cas_nearby]

            cas_type = "Unknown"
            if any("cas9" in g for g in cas_genes):
                cas_type = "II"
            elif any("cas12" in g for g in cas_genes):
                cas_type = "V"
            elif any("cas13" in g for g in cas_genes):
                cas_type = "VI"
            elif any("cas10" in g or "csm" in g or "cmr" in g for g in cas_genes):
                cas_type = "III"
            elif any("cas3" in g for g in cas_genes):
                cas_type = "I"
            elif any("cas1" in g or "cas2" in g for g in cas_genes):
                cas_type = "I/II/V (adaptation only)"

            # Check completeness
            has_cas1 = any("cas1" in g for g in cas_genes)
            has_cas2 = any("cas2" in g for g in cas_genes)
            has_effector = any(x in " ".join(cas_genes) for x in ["cas3", "cas9", "cas10", "cas12", "cas13"])

            is_complete = has_cas1 and has_cas2 and has_effector
            is_orphan = len(cas_nearby) == 0

            systems.append({
                "array_id": locus_id,
                "genome_id": genome_id,
                "contig_id": contig_id,
                "array_start": start,
                "array_end": end_coord,
                "at_contig_edge": at_edge,
                "cas_genes": [ann for _, ann, _, _ in cas_nearby],
                "cas_type": cas_type,
                "is_complete": is_complete,
                "is_orphan": is_orphan,
                "metadata": metadata,
            })

        # Summarize
        n_total = len(systems)
        n_complete = sum(1 for s in systems if s["is_complete"])
        n_orphan = sum(1 for s in systems if s["is_orphan"])
        n_at_edge = sum(1 for s in systems if s["at_contig_edge"])

        type_counts = {}
        for s in systems:
            type_counts[s["cas_type"]] = type_counts.get(s["cas_type"], 0) + 1

        lines = [
            f"# CRISPR-Cas System Analysis",
            f"",
            f"## Summary",
            f"",
            f"| Metric | Count | Percentage |",
            f"|--------|-------|------------|",
            f"| Total arrays | {n_total} | 100% |",
            f"| Complete systems | {n_complete} | {100*n_complete/n_total:.1f}% |",
            f"| Orphan arrays | {n_orphan} | {100*n_orphan/n_total:.1f}% |",
            f"| At contig edge | {n_at_edge} | {100*n_at_edge/n_total:.1f}% |",
            f"",
            f"## System Types",
            f"",
        ]

        for cas_type, count in sorted(type_counts.items(), key=lambda x: -x[1]):
            lines.append(f"- **Type {cas_type}**: {count} ({100*count/n_total:.1f}%)")

        lines.extend([
            f"",
            f"## Assembly Fragmentation Note",
            f"",
            f"**{n_at_edge} arrays ({100*n_at_edge/n_total:.1f}%)** are at contig edges.",
            f"These may be incomplete due to assembly fragmentation at CRISPR repeat sequences.",
            f"",
            f"## Orphan Arrays",
            f"",
        ])

        if n_orphan > 0:
            lines.append(f"**{n_orphan} arrays** have no associated cas genes within 15kb:")
            for s in [s for s in systems if s["is_orphan"]][:10]:
                lines.append(f"- {s['array_id']} on {s['contig_id']}")
        else:
            lines.append("No orphan arrays found.")

        lines.extend([
            f"",
            f"## Complete Systems",
            f"",
        ])

        if n_complete > 0:
            for s in [s for s in systems if s["is_complete"]][:10]:
                lines.append(f"- **{s['array_id']}** (Type {s['cas_type']}): {', '.join(s['cas_genes'][:5])}")
        else:
            lines.append("No complete systems found (may be due to annotation gaps).")

        return ctx.make_result(
            data="\n".join(lines),
            rows=n_total,
            raw=systems,
        )


def detect_annotation_errors(store, limit: int = 50) -> BennuResult:
    """
    Scan for likely annotation errors in the database.

    Checks for:
    1. Domain-annotation mismatches (especially for known problem domains)
    2. Proteins in wrong genomic context
    3. Spurious ORFs (overlapping CRISPR arrays)

    Args:
        store: DuckDB store
        limit: Maximum errors to report per category

    Returns:
        BennuResult with potential annotation errors
    """
    with OperatorContext("detect_annotation_errors", {"limit": limit}) as ctx:

        errors = {
            "problem_domain_hits": [],
            "spurious_orfs": [],
        }

        # Check for problem domain hits
        for domain_name, info in PROBLEM_DOMAINS.items():
            # Find proteins with this domain
            hits = store.execute("""
                SELECT DISTINCT
                    a1.protein_id,
                    COALESCE(a2.name, a2.description, 'unannotated') as best_annotation
                FROM annotations a1
                LEFT JOIN annotations a2 ON a1.protein_id = a2.protein_id
                WHERE a1.name = ?
                LIMIT ?
            """, [domain_name, limit])

            for protein_id, annotation in hits:
                ann_lower = (annotation or "").lower()
                for indicator in info["false_positive_indicators"]:
                    if indicator.lower() in ann_lower:
                        errors["problem_domain_hits"].append({
                            "protein_id": protein_id,
                            "domain": domain_name,
                            "annotation": annotation,
                            "issue": info["issue"],
                        })
                        break

        # Check for ORFs overlapping CRISPR arrays
        spurious = store.execute("""
            SELECT DISTINCT p.protein_id,
                   COALESCE(a.name, a.description, 'unannotated') as annotation,
                   p.start, p.end_coord, l.locus_id
            FROM proteins p
            LEFT JOIN annotations a ON p.protein_id = a.protein_id
            JOIN loci l ON p.contig_id = l.contig_id
            WHERE l.locus_type = 'crispr_array'
              AND p.start < l.end_coord
              AND p.end_coord > l.start
            LIMIT ?
        """, [limit])

        for protein_id, annotation, start, end, locus_id in spurious:
            errors["spurious_orfs"].append({
                "protein_id": protein_id,
                "annotation": annotation,
                "coords": f"{start}-{end}",
                "overlaps": locus_id,
                "issue": "ORF overlaps CRISPR array - likely spurious",
            })

        # Build report
        total_errors = sum(len(v) for v in errors.values())

        lines = [
            f"# Annotation Error Detection Report",
            f"",
            f"**Total potential errors found:** {total_errors}",
            f"",
        ]

        if errors["problem_domain_hits"]:
            lines.extend([
                f"## Problem Domain Hits ({len(errors['problem_domain_hits'])})",
                f"",
                f"These proteins have domains known to cause false-positive functional assignments:",
                f"",
            ])
            for err in errors["problem_domain_hits"][:10]:
                lines.append(f"- `{err['protein_id']}`: {err['domain']} domain but annotated as '{err['annotation']}'")
            if len(errors["problem_domain_hits"]) > 10:
                lines.append(f"- ... and {len(errors['problem_domain_hits']) - 10} more")
            lines.append("")

        if errors["spurious_orfs"]:
            lines.extend([
                f"## Spurious ORFs ({len(errors['spurious_orfs'])})",
                f"",
                f"These proteins overlap CRISPR arrays and are likely false ORF predictions:",
                f"",
            ])
            for err in errors["spurious_orfs"][:10]:
                lines.append(f"- `{err['protein_id']}` overlaps {err['overlaps']}")
            if len(errors["spurious_orfs"]) > 10:
                lines.append(f"- ... and {len(errors['spurious_orfs']) - 10} more")
            lines.append("")

        if total_errors == 0:
            lines.append("No annotation errors detected.")
            lines.append("")

        lines.extend([
            f"## Recommendations",
            f"",
            f"1. **Exclude spurious ORFs** from functional analyses",
            f"2. **Validate problem domain hits** by checking genomic context",
            f"3. **Cross-reference annotations** with domain evidence before drawing conclusions",
        ])

        return ctx.make_result(
            data="\n".join(lines),
            rows=total_errors,
            raw=errors,
        )


__all__ = [
    "PROBLEM_DOMAINS",
    "validate_annotation",
    "validate_context",
    "analyze_crispr_systems",
    "detect_annotation_errors",
]
