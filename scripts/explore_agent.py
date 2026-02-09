#!/usr/bin/env python3
"""
Agent-driven exploration script for Sharur.

This script provides a structured exploration workflow that an AI agent
can use to discover interesting loci in metagenomic data.

Usage:
    python scripts/explore_agent.py                    # General exploration
    python scripts/explore_agent.py --focus metabolism # Focus area
    python scripts/explore_agent.py --contig <id>     # Browse specific contig

The agent can import and call these functions directly:
    from scripts.explore_agent import overview, search, neighborhood, protein
"""

import argparse
import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from sharur.operators import Sharur

# Default database path
DEFAULT_DB = "data/sharur.duckdb"

# Global Sharur instance (initialized lazily)
_bennu = None


def get_bennu(db_path: str = DEFAULT_DB) -> Sharur:
    """Get or create Sharur instance."""
    global _bennu
    if _bennu is None:
        _bennu = Sharur(db_path)
    return _bennu


def overview(db_path: str = DEFAULT_DB):
    """Get dataset overview."""
    b = get_bennu(db_path)
    result = b.overview()
    print(result.data)
    return result


def search(has=None, lacks=None, limit=20, db_path: str = DEFAULT_DB):
    """Search proteins by predicates."""
    b = get_bennu(db_path)
    result = b.search_by_predicates(has=has, lacks=lacks, limit=limit)
    print(result.data)
    return result


def neighborhood(protein_id, window=10, db_path: str = DEFAULT_DB):
    """Get genomic neighborhood around a protein."""
    b = get_bennu(db_path)
    result = b.get_neighborhood(protein_id, window=window)
    print(result.data)
    return result


def protein(protein_id, db_path: str = DEFAULT_DB):
    """Get protein details."""
    b = get_bennu(db_path)
    result = b.get_protein(protein_id, verbosity=2)
    print(result.data)
    return result


def genomes(limit=10, db_path: str = DEFAULT_DB):
    """List genomes."""
    b = get_bennu(db_path)
    result = b.list_genomes(limit=limit)
    print(result.data)
    return result


def proteins_on_contig(contig_id, limit=100, db_path: str = DEFAULT_DB):
    """List proteins on a specific contig."""
    b = get_bennu(db_path)
    result = b.list_proteins(contig_id=contig_id, limit=limit)
    print(result.data)
    return result


def contigs(genome_id=None, limit=20, db_path: str = DEFAULT_DB):
    """List contigs, optionally filtered by genome."""
    b = get_bennu(db_path)
    result = b.list_contigs(genome_id=genome_id, limit=limit)
    print(result.data)
    return result


# Focus-specific search configurations
FOCUS_SEARCHES = {
    "metabolism": [
        {"has": ["hydrogenase"], "name": "Hydrogenases"},
        {"has": ["nitrogenase"], "name": "Nitrogenases"},
        {"has": ["methanogenesis"], "name": "Methanogenesis enzymes"},
        {"has": ["sulfur_metabolism"], "name": "Sulfur metabolism"},
        {"has": ["oxidoreductase", "unannotated"], "name": "Novel oxidoreductases"},
        {"has": ["iron_sulfur", "hypothetical"], "name": "Hypothetical Fe-S proteins"},
        {"has": ["molybdenum_binding"], "name": "Molybdenum enzymes"},
        {"has": ["nickel_binding"], "name": "Nickel enzymes"},
    ],
    "phage": [
        {"has": ["phage_related"], "name": "Phage-related proteins"},
        {"has": ["integrase"], "name": "Integrases"},
        {"has": ["giant", "phage_related"], "name": "Giant phage proteins"},
        {"has": ["unannotated", "large"], "lacks": ["transposase"], "name": "Large unannotated (potential phage)"},
    ],
    "defense": [
        {"has": ["crispr_associated"], "name": "CRISPR-associated"},
        {"has": ["restriction_modification"], "name": "Restriction-modification"},
        {"has": ["toxin_antitoxin"], "name": "Toxin-antitoxin systems"},
        {"has": ["abortive_infection"], "name": "Abortive infection"},
    ],
    "novel": [
        {"has": ["giant", "unannotated"], "name": "Giant unannotated"},
        {"has": ["massive"], "name": "Massive proteins (>2000 aa)"},
        {"has": ["multi_domain", "hypothetical"], "name": "Multi-domain hypotheticals"},
        {"has": ["transmembrane_predicted", "unannotated"], "name": "Novel membrane proteins"},
        {"has": ["metal_binding", "hypothetical"], "name": "Hypothetical metal-binding"},
    ],
    "general": [
        {"has": ["giant", "unannotated"], "name": "Giant unannotated proteins"},
        {"has": ["massive"], "name": "Massive proteins (>2000 aa)"},
        {"has": ["multi_domain", "hypothetical"], "name": "Multi-domain hypotheticals"},
        {"has": ["nickel_binding"], "name": "Nickel-binding proteins"},
        {"has": ["molybdenum_binding"], "name": "Molybdenum-binding proteins"},
        {"has": ["crispr_associated"], "name": "CRISPR-associated proteins"},
        {"has": ["phage_related", "giant"], "name": "Giant phage proteins"},
        {"has": ["transmembrane_predicted", "unannotated"], "name": "Novel membrane proteins"},
    ],
}


def run_focused_exploration(focus: str = "general", db_path: str = DEFAULT_DB):
    """Run exploration with a specific focus area."""
    b = get_bennu(db_path)

    print("=" * 70)
    print(f"BENNU EXPLORATION: {focus.upper()}")
    print("=" * 70)

    # Overview
    print("\n## DATASET OVERVIEW\n")
    overview(db_path)

    # Get searches for this focus
    searches = FOCUS_SEARCHES.get(focus, FOCUS_SEARCHES["general"])

    findings = []
    for search_def in searches:
        name = search_def["name"]
        has = search_def.get("has", [])
        lacks = search_def.get("lacks", [])

        print(f"\n## {name.upper()}\n")
        result = b.search_by_predicates(has=has, lacks=lacks, limit=5)

        if result.meta.rows > 0:
            print(result.data)
            findings.append({
                "name": name,
                "count": result.meta.total_rows or result.meta.rows,
                "has": has,
                "lacks": lacks,
            })

            # Show neighborhood for first interesting hit
            if result.meta.rows > 0:
                # Try to get first protein ID and show neighborhood
                try:
                    # Parse the result to get first protein
                    lines = result.data.split("\n")
                    for line in lines[2:]:  # Skip header lines
                        parts = line.split()
                        if parts and "_" in parts[0]:
                            protein_id = parts[0].strip()
                            print(f"\n### Neighborhood of {protein_id}\n")
                            neighborhood(protein_id, window=5, db_path=db_path)
                            break
                except Exception:
                    pass

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY OF FINDINGS")
    print("=" * 70)
    for f in findings:
        print(f"- {f['name']}: {f['count']} proteins")

    return findings


def browse_contig(contig_id: str, window_size: int = 15, db_path: str = DEFAULT_DB):
    """Browse a contig window-by-window."""
    b = get_bennu(db_path)

    print("=" * 70)
    print(f"BROWSING CONTIG: {contig_id}")
    print("=" * 70)

    # Get all proteins on contig
    result = b.list_proteins(contig_id=contig_id, limit=500)
    print(f"\n## All proteins on contig\n")
    print(result.data)

    # Parse protein count
    protein_count = result.meta.rows
    print(f"\nTotal proteins: {protein_count}")
    print(f"Window size: {window_size}")
    print("\nUse neighborhood() to examine specific regions.")

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Agent-driven exploration of Sharur metagenomic data"
    )
    parser.add_argument(
        "--focus",
        choices=["metabolism", "phage", "defense", "novel", "general"],
        default="general",
        help="Focus area for exploration"
    )
    parser.add_argument(
        "--contig",
        type=str,
        help="Browse a specific contig"
    )
    parser.add_argument(
        "--db",
        type=str,
        default=DEFAULT_DB,
        help="Path to Sharur database"
    )

    args = parser.parse_args()

    if args.contig:
        browse_contig(args.contig, db_path=args.db)
    else:
        run_focused_exploration(args.focus, db_path=args.db)


if __name__ == "__main__":
    main()
