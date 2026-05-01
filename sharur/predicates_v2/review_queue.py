"""
Review queue: surfaces unmapped accessions for curation.

Scans all unresolved atoms, aggregates counts per accession,
prioritizes by frequency, and outputs a TSV for human review.
"""

from __future__ import annotations

import csv
import io
from typing import Any

from sharur.predicates_v2.model import ClaimRelation, SemanticAtom


# Keyword -> suggested facet mapping for auto-suggesting facets
_FACET_KEYWORDS: dict[str, str] = {
    "kinase": "activity",
    "protease": "activity",
    "hydrolase": "activity",
    "transferase": "activity",
    "oxidoreductase": "activity",
    "synthase": "activity",
    "reductase": "activity",
    "dehydrogenase": "activity",
    "lyase": "activity",
    "ligase": "activity",
    "phosphatase": "activity",
    "polymerase": "activity",
    "nuclease": "activity",
    "esterase": "activity",
    "transporter": "role",
    "transport": "role",
    "regulator": "role",
    "regulation": "role",
    "defense": "role",
    "toxin": "role",
    "resistance": "role",
    "membrane": "localization",
    "periplasmic": "localization",
    "outer membrane": "localization",
    "secreted": "localization",
    "repeat": "architecture",
    "domain": "architecture",
    "coiled": "architecture",
    "barrel": "architecture",
}


def build_review_queue(
    all_atoms: list[SemanticAtom],
    protein_genomes: dict[str, str] | None = None,
    total_proteins: int = 0,
) -> list[dict[str, Any]]:
    """Build a prioritized review queue from unresolved atoms.

    Args:
        all_atoms: All atoms across all proteins.
        protein_genomes: Optional mapping protein_id -> genome_id
            for computing genome spread.
        total_proteins: Total number of proteins in the dataset.

    Returns:
        List of review queue entries, sorted by priority (descending).
    """
    # Collect unresolved accessions
    unresolved: dict[tuple[str, str], dict[str, Any]] = {}

    for atom in all_atoms:
        if atom.relation != ClaimRelation.unresolved:
            continue

        acc = atom.source_accession
        source_db = atom.source_db
        key = (source_db, acc)
        if key not in unresolved:
            unresolved[key] = {
                "accession": acc,
                "source_db": source_db,
                "protein_ids": set(),
                "genome_ids": set(),
                "example_protein_ids": [],
            }

        entry = unresolved[key]
        entry["protein_ids"].add(atom.protein_id)

        # Track genomes if mapping provided
        if protein_genomes and atom.protein_id in protein_genomes:
            entry["genome_ids"].add(protein_genomes[atom.protein_id])

        # Keep up to 5 example protein IDs
        if (
            len(entry["example_protein_ids"]) < 5
            and atom.protein_id not in entry["example_protein_ids"]
        ):
            entry["example_protein_ids"].append(atom.protein_id)

    # Build output rows
    rows = []
    for (_source_db, acc), entry in unresolved.items():
        n_proteins = len(entry["protein_ids"])
        n_genomes = len(entry["genome_ids"])
        pct_proteome = (n_proteins / total_proteins * 100) if total_proteins > 0 else 0.0

        # Suggest a facet based on accession name
        suggested_facet = suggest_facet(acc)

        rows.append({
            "accession": acc,
            "source_db": entry["source_db"],
            "n_proteins": n_proteins,
            "n_genomes": n_genomes,
            "pct_proteome": round(pct_proteome, 2),
            "example_protein_ids": ";".join(entry["example_protein_ids"]),
            "suggested_facet": suggested_facet,
            # Priority score: proteins * genomes
            "_priority": n_proteins * max(n_genomes, 1),
        })

    # Sort by priority descending
    rows.sort(key=lambda r: r["_priority"], reverse=True)

    return rows


def suggest_facet(accession: str) -> str:
    """Suggest a facet based on accession name keywords.

    Returns empty string if no suggestion.
    """
    lower = accession.lower()
    for keyword, facet in _FACET_KEYWORDS.items():
        if keyword in lower:
            return facet
    return ""


_suggest_facet = suggest_facet


def format_review_queue_tsv(rows: list[dict[str, Any]]) -> str:
    """Format review queue as TSV string.

    Args:
        rows: Review queue entries from build_review_queue().

    Returns:
        TSV-formatted string with header.
    """
    output = io.StringIO()
    writer = csv.writer(output, delimiter="\t")

    # Header
    writer.writerow([
        "accession",
        "source_db",
        "n_proteins",
        "n_genomes",
        "pct_proteome",
        "example_protein_ids",
        "suggested_facet",
    ])

    for row in rows:
        writer.writerow([
            row["accession"],
            row["source_db"],
            row["n_proteins"],
            row["n_genomes"],
            row["pct_proteome"],
            row["example_protein_ids"],
            row["suggested_facet"],
        ])

    return output.getvalue()


__all__ = [
    "build_review_queue",
    "format_review_queue_tsv",
    "suggest_facet",
]
