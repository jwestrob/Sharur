#!/usr/bin/env python3
"""
Hydrogenase Classification Pipeline

Two-stage classification:
1. HMM search with HydDB profiles (via Astra) - identifies hydrogenase class (NiFe, FeFe, Fe)
2. DIAMOND search against HydDB reference - assigns subgroup (e.g., Group_1a, Group_4e)

All DIAMOND-classified hits receive subgroup predicates. PFAM domains are
recorded as metadata but do NOT gate predicate assignment — this avoids
silently dropping Group 4 NiFe hydrogenases (Hyf/Hyc/Mbh/Ech), which
diverged too far from the NiFeSe_Hases domain (PF00374).

Hits that lack PFAM corroboration or show Complex I markers are tagged
with `hyddb_needs_curation` so agents can prioritize neighborhood-based
validation during analysis. See .claude/skills/survey.md for the
hydrogenase curation protocol.

Usage:
    python scripts/classify_hydrogenases.py --db data/my_dataset/sharur.duckdb

References:
    Søndergaard D, et al. (2016) HydDB: A web tool for hydrogenase classification
    and analysis. Sci Rep 6:34212. doi:10.1038/srep34212
"""

import argparse
import subprocess
import tempfile
from pathlib import Path
import duckdb
import pandas as pd


# HydDB subgroup mappings for predicates
# Predicate names must match vocabulary.py
NIFE_SUBGROUPS = {
    # Group 1: Membrane-bound H2 uptake
    "Group_1a": ("nife_group1", "uptake_hydrogenase", "membrane_bound_hydrogenase"),
    "Group_1b": ("nife_group1", "uptake_hydrogenase", "membrane_bound_hydrogenase"),
    "Group_1c": ("nife_group1", "uptake_hydrogenase", "membrane_bound_hydrogenase"),
    "Group_1d": ("nife_group1", "uptake_hydrogenase", "membrane_bound_hydrogenase"),
    "Group_1e": ("nife_group1", "uptake_hydrogenase"),
    "Group_1f": ("nife_group1", "uptake_hydrogenase"),
    "Group_1g": ("nife_group1", "uptake_hydrogenase"),
    "Group_1h": ("nife_group1", "uptake_hydrogenase"),  # Actinobacterial
    "Group_1i": ("nife_group1", "uptake_hydrogenase"),
    "Group_1j": ("nife_group1", "uptake_hydrogenase"),
    "Group_1k": ("nife_group1", "uptake_hydrogenase"),
    "Group_1l": ("nife_group1", "uptake_hydrogenase"),
    # Group 2: Cytoplasmic H2 sensors and uptake
    "Group_2a": ("nife_group2", "cytoplasmic_hydrogenase", "h2_sensor"),
    "Group_2b": ("nife_group2", "cytoplasmic_hydrogenase", "h2_sensor"),
    "Group_2c": ("nife_group2", "cytoplasmic_hydrogenase"),
    "Group_2d": ("nife_group2", "cytoplasmic_hydrogenase"),
    "Group_2e": ("nife_group2", "cytoplasmic_hydrogenase"),
    # Group 3: Bidirectional/F420-reducing
    "Group_3a": ("nife_group3", "bidirectional_hydrogenase", "f420_reducing"),
    "Group_3b": ("nife_group3", "bidirectional_hydrogenase", "nadp_reducing"),
    "Group_3c": ("nife_group3", "bidirectional_hydrogenase", "methyl_viologen_reducing"),
    "Group_3d": ("nife_group3", "bidirectional_hydrogenase", "nadp_reducing"),
    # Group 4: H2-evolving, energy-conserving
    "Group_4a": ("nife_group4", "h2_evolving", "energy_conserving_hydrogenase", "formate_coupled"),
    "Group_4b": ("nife_group4", "h2_evolving", "energy_conserving_hydrogenase", "formate_coupled"),
    "Group_4c": ("nife_group4", "h2_evolving", "energy_conserving_hydrogenase", "co_coupled"),
    "Group_4d": ("nife_group4", "h2_evolving", "energy_conserving_hydrogenase", "co_coupled"),
    "Group_4e": ("nife_group4", "h2_evolving", "energy_conserving_hydrogenase", "mbh_hydrogenase"),
    "Group_4f": ("nife_group4", "h2_evolving", "energy_conserving_hydrogenase"),
    "Group_4g": ("nife_group4", "h2_evolving", "energy_conserving_hydrogenase", "ech_hydrogenase"),
    "Group_4h": ("nife_group4", "h2_evolving", "energy_conserving_hydrogenase"),
    "Group_4i": ("nife_group4", "h2_evolving", "energy_conserving_hydrogenase"),
}

FEFE_SUBGROUPS = {
    # Group A: Monomeric, cytoplasmic
    "Group_A1": ("fefe_groupA", "monomeric_fefe", "fermentative_hydrogenase"),
    "Group_A2": ("fefe_groupA", "monomeric_fefe", "fermentative_hydrogenase"),
    "Group_A3": ("fefe_groupA", "monomeric_fefe", "fermentative_hydrogenase"),
    "Group_A4": ("fefe_groupA", "monomeric_fefe"),
    # Group B: Bifurcating
    "Group_B": ("fefe_groupB", "bifurcating_hydrogenase"),
    # Group C: Sensory/regulatory
    "Group_C1": ("fefe_groupC", "sensory_hydrogenase"),
    "Group_C2": ("fefe_groupC", "sensory_hydrogenase"),
    "Group_C3": ("fefe_groupC", "sensory_hydrogenase"),
}


def get_hyddb_reference_path() -> Path:
    """Get path to HydDB reference files."""
    # Check standard locations
    candidates = [
        Path(__file__).parent.parent / "data/reference/hyddb",
        Path("data/reference/hyddb"),
        Path.home() / ".sharur/hyddb",
    ]
    for path in candidates:
        if (path / "HydDB_all.dmnd").exists():
            return path
    raise FileNotFoundError(
        "HydDB reference not found. Run: diamond makedb --in HydDB_all_hydrogenases.faa --db HydDB_all"
    )


def classify_with_diamond(
    protein_ids: list[str],
    sequences: dict[str, str],
    hyddb_path: Path,
    threads: int = 4,
) -> dict[str, dict]:
    """
    Classify hydrogenases by DIAMOND search against HydDB reference.

    Returns dict mapping protein_id to classification info.
    """
    if not protein_ids:
        return {}

    # Write sequences to temp file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.faa', delete=False) as f:
        for pid in protein_ids:
            if pid in sequences:
                f.write(f">{pid}\n{sequences[pid]}\n")
        query_file = f.name

    # Run DIAMOND
    db_path = hyddb_path / "HydDB_all.dmnd"

    try:
        result = subprocess.run(
            [
                "diamond", "blastp",
                "--db", str(db_path),
                "--query", query_file,
                "--outfmt", "6", "qseqid", "sseqid", "pident", "evalue", "bitscore",
                "--max-target-seqs", "1",
                "--threads", str(threads),
                "--ultra-sensitive",
            ],
            capture_output=True,
            text=True,
            timeout=300,
        )
    except subprocess.TimeoutExpired:
        print("Warning: DIAMOND search timed out")
        return {}
    finally:
        Path(query_file).unlink(missing_ok=True)

    if result.returncode != 0:
        print(f"Warning: DIAMOND error: {result.stderr}")
        return {}

    # Parse results
    classifications = {}
    for line in result.stdout.strip().split('\n'):
        if not line:
            continue
        parts = line.split('\t')
        if len(parts) >= 5:
            qid, sid, pident, evalue, bitscore = parts[:5]

            # Parse subject ID to get classification
            # Format: WP_xxx|Organism_name|[Type]_Group_Xx
            sid_parts = sid.split('|')
            if len(sid_parts) >= 3:
                classification = sid_parts[2]  # e.g., "[NiFe]_Group_1a"

                # Extract type and subgroup
                if "[NiFe]_" in classification:
                    hyd_type = "NiFe"
                    subgroup = classification.replace("[NiFe]_", "")
                elif "[FeFe]_" in classification:
                    hyd_type = "FeFe"
                    subgroup = classification.replace("[FeFe]_", "")
                elif "[Fe]" in classification:
                    hyd_type = "Fe"
                    subgroup = "Fe_only"
                else:
                    continue

                classifications[qid] = {
                    "type": hyd_type,
                    "subgroup": subgroup,
                    "reference_id": sid_parts[0],
                    "organism": sid_parts[1] if len(sid_parts) > 1 else "",
                    "pident": float(pident),
                    "evalue": float(evalue),
                    "bitscore": float(bitscore),
                }

    return classifications


def get_subgroup_predicates(hyd_type: str, subgroup: str) -> list[str]:
    """Get semantic predicates for a hydrogenase subgroup."""
    predicates = []

    if hyd_type == "NiFe":
        if subgroup in NIFE_SUBGROUPS:
            predicates.extend(NIFE_SUBGROUPS[subgroup])
        # Add generic group predicate
        group_num = subgroup.split('_')[1][0] if '_' in subgroup else None
        if group_num:
            predicates.append(f"nife_group{group_num}")

    elif hyd_type == "FeFe":
        if subgroup in FEFE_SUBGROUPS:
            predicates.extend(FEFE_SUBGROUPS[subgroup])

    elif hyd_type == "Fe":
        predicates.extend(["fe_only_hydrogenase", "methanogen_hydrogenase"])

    return list(set(predicates))


def classify_hydrogenases(
    db_path: str,
    threads: int = 4,
    update_predicates: bool = True,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Classify hydrogenases in a Sharur database.

    Args:
        db_path: Path to sharur.duckdb
        threads: Number of threads for DIAMOND
        update_predicates: Whether to update predicates in database
        verbose: Print progress

    Returns:
        DataFrame with classification results
    """
    db = duckdb.connect(db_path)

    # Find HydDB reference
    try:
        hyddb_path = get_hyddb_reference_path()
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return pd.DataFrame()

    if verbose:
        print(f"Using HydDB reference: {hyddb_path}")

    # Get proteins with HydDB annotations
    hyddb_proteins = db.execute("""
        SELECT DISTINCT a.protein_id, a.accession as hmm_type, a.score as hmm_score
        FROM annotations a
        WHERE LOWER(a.source) = 'hyddb'
    """).fetchdf()

    if hyddb_proteins.empty:
        print("No HydDB annotations found. Run Astra with HydDB first.")
        return pd.DataFrame()

    if verbose:
        print(f"Found {len(hyddb_proteins)} proteins with HydDB annotations")

    # Get sequences
    sequences = {}
    for row in db.execute("""
        SELECT protein_id, sequence FROM proteins WHERE sequence IS NOT NULL
    """).fetchall():
        sequences[row[0]] = row[1]

    # Get PFAM annotations for validation
    pfam_domains = {}
    for row in db.execute("""
        SELECT protein_id, accession FROM annotations WHERE LOWER(source) = 'pfam'
    """).fetchall():
        if row[0] not in pfam_domains:
            pfam_domains[row[0]] = set()
        pfam_domains[row[0]].add(row[1])

    # Classify with DIAMOND
    protein_ids = hyddb_proteins['protein_id'].tolist()
    if verbose:
        print(f"Running DIAMOND classification for {len(protein_ids)} proteins...")

    diamond_results = classify_with_diamond(protein_ids, sequences, hyddb_path, threads)

    if verbose:
        print(f"DIAMOND classified {len(diamond_results)} proteins")

    # Build results — classify all hits, flag those needing curation
    results = []
    for _, row in hyddb_proteins.iterrows():
        pid = row['protein_id']
        hmm_type = row['hmm_type']
        hmm_score = row['hmm_score']

        # Get DIAMOND classification
        diamond = diamond_results.get(pid, {})
        subgroup = diamond.get('subgroup', 'unknown')
        pident = diamond.get('pident', 0)

        # Check PFAM context (informational, not gating)
        domains = pfam_domains.get(pid, set())
        has_nifese = 'PF00374' in domains  # NiFeSe_Hases
        has_fe_hyd = 'PF02906' in domains or 'PF02256' in domains  # Fe_hyd domains
        has_complex1 = 'PF00346' in domains or 'PF00329' in domains  # Complex I

        # Determine confidence level — all hits get classified,
        # but some are flagged for agent-level neighborhood curation
        if hmm_type == 'NiFe':
            if has_nifese:
                needs_curation = False
                pfam_note = "Corroborated by NiFeSe_Hases (PF00374)"
            elif has_complex1 and not has_nifese:
                needs_curation = True
                pfam_note = "Likely Complex I — has PF00346/PF00329 without PF00374"
            else:
                needs_curation = True
                pfam_note = "No PF00374 — likely Group 4 or needs neighborhood check"
        elif hmm_type == 'FeFe':
            if has_fe_hyd:
                needs_curation = False
                pfam_note = "Corroborated by Fe_hyd domain (PF02906/PF02256)"
            else:
                needs_curation = True
                pfam_note = "No Fe_hyd domain — needs neighborhood check"
        else:  # Fe_only
            needs_curation = False
            pfam_note = "Fe-only (rare, trusted)"

        results.append({
            'protein_id': pid,
            'hmm_type': hmm_type,
            'hmm_score': hmm_score,
            'subgroup': subgroup,
            'diamond_pident': pident,
            'needs_curation': needs_curation,
            'pfam_note': pfam_note,
            'has_nifese_hases': has_nifese,
            'has_fe_hyd': has_fe_hyd,
            'has_complex1': has_complex1,
        })

    results_df = pd.DataFrame(results)

    # Print summary
    if verbose:
        print("\n=== CLASSIFICATION SUMMARY ===")
        print(f"\nBy HMM type:")
        print(results_df.groupby('hmm_type')['protein_id'].count())
        print(f"\nCuration status:")
        print(results_df.groupby(['hmm_type', 'needs_curation'])['protein_id'].count())
        print(f"\nSubgroups (all classified):")
        print(results_df.groupby(['hmm_type', 'subgroup'])['protein_id'].count().head(30))
        n_curation = results_df['needs_curation'].sum()
        if n_curation > 0:
            print(f"\n{n_curation} hits flagged for agent neighborhood curation")

    # Update predicates if requested — ALL hits get subgroup predicates
    if update_predicates:
        if verbose:
            print("\nUpdating predicates...")

        for _, row in results_df.iterrows():
            pid = row['protein_id']
            subgroup = row['subgroup']
            hmm_type = row['hmm_type']

            # Get subgroup predicates
            new_predicates = get_subgroup_predicates(hmm_type, subgroup)
            if not new_predicates:
                continue

            # Add subgroup direct access predicate
            new_predicates.append(f"hyddb_subgroup:{subgroup}")

            # Flag hits that need agent-level neighborhood curation
            if row['needs_curation']:
                new_predicates.append("hyddb_needs_curation")

            # Get current predicates
            current = db.execute("""
                SELECT predicates FROM protein_predicates WHERE protein_id = ?
            """, [pid]).fetchone()

            if current:
                current_preds = set(current[0])
                current_preds.update(new_predicates)
                db.execute("""
                    UPDATE protein_predicates
                    SET predicates = ?, updated_at = CURRENT_TIMESTAMP
                    WHERE protein_id = ?
                """, [list(current_preds), pid])

        db.commit()
        if verbose:
            print("Predicates updated")

    # Save detailed results
    output_path = Path(db_path).parent / "hydrogenase_classification.tsv"
    results_df.to_csv(output_path, sep='\t', index=False)
    if verbose:
        print(f"\nDetailed results saved to: {output_path}")

    db.close()
    return results_df


def main():
    parser = argparse.ArgumentParser(
        description="Classify hydrogenases using HydDB",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--db", required=True,
        help="Path to sharur.duckdb"
    )
    parser.add_argument(
        "--threads", type=int, default=4,
        help="Number of threads for DIAMOND (default: 4)"
    )
    parser.add_argument(
        "--no-update", action="store_true",
        help="Don't update predicates in database"
    )
    parser.add_argument(
        "--quiet", action="store_true",
        help="Suppress progress output"
    )

    args = parser.parse_args()

    results = classify_hydrogenases(
        db_path=args.db,
        threads=args.threads,
        update_predicates=not args.no_update,
        verbose=not args.quiet,
    )

    if results.empty:
        print("No hydrogenases classified")
        return 1

    return 0


if __name__ == "__main__":
    exit(main())
