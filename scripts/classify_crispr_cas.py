#!/usr/bin/env python3
"""
CRISPR-Cas Classification Pipeline

Classifies CRISPR-Cas systems by type and subtype using CRISPRCasFinder HMMs.

Classification hierarchy:
- Class 1 (multi-subunit): Types I, III, IV
- Class 2 (single-effector): Types II, V, VI

Validation:
- Proteins near CRISPR arrays (<20kb) get higher confidence
- Core proteins (Cas1/Cas2) are marked as adaptation module
- Signature proteins determine type (Cas3=I, Cas9=II, Cas10=III, Csf=IV, Cas12=V, Cas13=VI)

Usage:
    python scripts/classify_crispr_cas.py --db data/my_dataset/sharur.duckdb

Requires:
    - CRISPRCasFinder results in {db_dir}/crisprcasfinder_results/
    - CRISPR arrays in loci table (from MinCED)
"""

import argparse
import re
from pathlib import Path
import duckdb
import pandas as pd


# TIGR ID to CRISPR type mapping from NCBI CDD/TIGRFAMs
# Source: https://ftp.ncbi.nih.gov/pub/wolf/_suppl/CRISPRclass/crisprPro.html
TIGR_TO_CRISPR = {
    # Cas1 - adaptation module (core, found in most types)
    "TIGR00287": ("Core", None, "adaptation"),  # cas1
    "TIGR03641": ("I", "B", "adaptation"),       # cas1_HMARI
    "TIGR03640": ("I", "C", "adaptation"),       # cas1_DVULG
    "TIGR03638": ("I", "E", "adaptation"),       # cas1_ECOLI
    "TIGR03637": ("I", "F", "adaptation"),       # cas1_YPEST
    "TIGR03639": ("I", "A", "adaptation"),       # cas1 Type I-A variant

    # Cas2 - adaptation module (core)
    "TIGR01573": ("Core", None, "adaptation"),   # cas2
    "TIGR01877": ("Core", None, "adaptation"),   # cas2 variant
    "TIGR01875": ("I", "E", "adaptation"),       # cas2 I-E
    "TIGR01876": ("I", "F", "adaptation"),       # cas2 I-F

    # Cas3 - Type I signature (helicase-nuclease)
    "TIGR01587": ("I", None, "effector"),        # cas3_core (helicase domain)
    "TIGR01596": ("I", None, "effector"),        # cas3_HD (nuclease domain)
    "TIGR02621": ("I", None, "effector"),        # cas3_GSU0051
    "TIGR03158": ("I", "D", "effector"),         # cas3_cyano
    "TIGR02562": ("I", "F", "effector"),         # cas3_yersinia
    "TIGR01595": ("I", None, "effector"),        # cas3 variant

    # Cas4 - accessory (found in I-A, I-B, I-C, I-D, II-B, V)
    "TIGR00372": ("Core", None, "accessory"),    # cas4
    "TIGR02574": ("Core", None, "accessory"),    # cas4 variant

    # Type I-A backbone (Csa)
    "TIGR01896": ("I", "A", "backbone"),         # cas_AF1879
    "TIGR01894": ("I", "A", "backbone"),         # csa2
    "TIGR01895": ("I", "A", "backbone"),         # csa4
    "TIGR01899": ("I", "A", "backbone"),         # csa5
    "TIGR01908": ("I", "A", "backbone"),         # cas8a1

    # Type I-B backbone
    "TIGR02556": ("I", "B", "backbone"),         # cas8b
    "TIGR02557": ("I", "B", "backbone"),         # cas5 I-B
    "TIGR02558": ("I", "B", "backbone"),         # cas7 I-B
    "TIGR02165": ("I", "B", "backbone"),         # cas6 I-B

    # Type I-C backbone
    "TIGR02577": ("I", "C", "backbone"),         # cas8c
    "TIGR02578": ("I", "C", "backbone"),         # cas5 I-C
    "TIGR02580": ("I", "C", "backbone"),         # cas7 I-C

    # Type I-D backbone
    "TIGR02585": ("I", "D", "backbone"),         # cas10d
    "TIGR02589": ("I", "D", "backbone"),         # cas7 I-D
    "TIGR02590": ("I", "D", "backbone"),         # cas5 I-D
    "TIGR02591": ("I", "D", "backbone"),         # cas6 I-D
    "TIGR02592": ("I", "D", "backbone"),         # csc1
    "TIGR02593": ("I", "D", "backbone"),         # csc2

    # Type I-E backbone (CRISPR/Cascade)
    "TIGR01873": ("I", "E", "backbone"),         # CT1978/cse2
    "TIGR01863": ("I", "E", "backbone"),         # cas5 I-E
    "TIGR01862": ("I", "E", "backbone"),         # cas6 I-E
    "TIGR02547": ("I", "E", "backbone"),         # cas7 I-E
    "TIGR01870": ("I", "E", "backbone"),         # cse1/cas8e
    "TIGR01871": ("I", "E", "backbone"),         # cse4

    # Type I-F backbone
    "TIGR01903": ("I", "F", "backbone"),         # csy1/cas8f
    "TIGR02548": ("I", "F", "backbone"),         # csy2/cas5 I-F
    "TIGR02549": ("I", "F", "backbone"),         # csy3/cas7 I-F
    "TIGR01865": ("I", "F", "backbone"),         # csy4/cas6f

    # Type II (Cas9-based)
    "TIGR01866": ("II", None, "backbone"),       # csn1 (II-associated)
    "TIGR01867": ("II", "A", "backbone"),        # csn2 II-A
    "TIGR02584": ("II", "C", "backbone"),        # cas9 II-C variant

    # Type III-A backbone (Csm)
    "TIGR02570": ("III", "A", "backbone"),       # cas10 III-A
    "TIGR01903": ("III", "A", "backbone"),       # csm2
    "TIGR01898": ("III", "A", "backbone"),       # csm3
    "TIGR01899": ("III", "A", "backbone"),       # csm4
    "TIGR01900": ("III", "A", "backbone"),       # csm5
    "TIGR01901": ("III", "A", "backbone"),       # csm6

    # Type III-B backbone (Cmr)
    "TIGR02582": ("III", "B", "backbone"),       # cas10 III-B / cmr2
    "TIGR03174": ("III", "B", "backbone"),       # cmr1
    "TIGR01873": ("III", "B", "backbone"),       # cmr3
    "TIGR01874": ("III", "B", "backbone"),       # cmr4
    "TIGR01878": ("III", "B", "backbone"),       # cmr5
    "TIGR01879": ("III", "B", "backbone"),       # cmr6

    # Type III-D backbone
    "TIGR02674": ("III", "D", "backbone"),       # csx10

    # Type IV
    "TIGR04327": ("IV", None, "backbone"),       # csf1
    "TIGR04328": ("IV", None, "backbone"),       # csf2
    "TIGR04329": ("IV", None, "backbone"),       # csf3
    "TIGR04106": ("IV", None, "backbone"),       # csf4

    # Type V (Cas12)
    "TIGR04330": ("V", "A", "effector"),         # cpf1/cas12a
    "C2c3": ("V", "C", "effector"),              # cas12c

    # Additional Type I-E backbone
    "TIGR01869": ("I", "E", "backbone"),         # Cas7/Cse4/CasC
    "TIGR01868": ("I", "E", "backbone"),         # Cas5/CasD
    "TIGR01907": ("I", "E", "backbone"),         # Cas6/Cse3/CasE

    # Type I-F additional
    "TIGR02564": ("I", "F", "backbone"),         # Csy1
    "TIGR02565": ("I", "F", "backbone"),         # Csy2
    "TIGR02566": ("I", "F", "backbone"),         # Csy3

    # MYXAN subtype (Type I variant)
    "TIGR02807": ("I", None, "backbone"),        # Cas6 MYXAN subtype
    "TIGR03485": ("I", None, "backbone"),        # Cas8a1/Csx13 MYXAN

    # Dpsyc subtype (Type I variant)
    "TIGR04113": ("I", None, "backbone"),        # Csx17 Dpsyc subtype

    # Type IV additional
    "TIGR03115": ("IV", None, "backbone"),       # Csf2 AFERR
    "TIGR03114": ("IV", None, "backbone"),       # Csf3 AFERR
    "TIGR03116": ("IV", None, "backbone"),       # Csf1 AFERR

    # General CRISPR-associated (untyped)
    "TIGR03486": ("Core", None, "accessory"),    # cas_Cse1-associated
    "TIGR03487": ("Core", None, "accessory"),    # CRISPR-associated protein
    "TIGR03489": ("Core", None, "accessory"),    # CRISPR-associated protein
    "TIGR03983": ("Core", None, "accessory"),    # CRISPR-assoc. protein
    "TIGR03117": ("Core", None, "accessory"),    # CRISPR-assoc. DxTHG
    "TIGR02672": ("III", None, "backbone"),      # CRISPR Type III-associated
}


def _parse_roman_type(code: str) -> str:
    """Convert Roman numeral portion to type string."""
    code = code.upper()
    if code.startswith("VI"):
        return "VI"
    if code.startswith("IV"):
        return "IV"
    if code.startswith("III"):
        return "III"
    if code.startswith("II"):
        return "II"
    if code.startswith("V"):
        return "V"
    if code.startswith("I"):
        return "I"
    return "Unknown"


def _parse_subtype_letter(code: str) -> str | None:
    """Extract subtype letter (A-F) from subtype code."""
    code = code.upper()
    # Remove Roman numerals
    for numeral in ["VI", "IV", "III", "II", "V", "I"]:
        if code.startswith(numeral):
            remainder = code[len(numeral):]
            if remainder and remainder[0] in "ABCDEF":
                return remainder[0]
            break
    return None


def _infer_protein_class(hmm: str) -> str:
    """Infer protein class from HMM name."""
    hmm = hmm.lower()
    # Effector proteins
    if any(x in hmm for x in ['cas3', 'cas9', 'cas10', 'cas12', 'cas13', 'cpf1']):
        return "effector"
    # Adaptation module
    if any(x in hmm for x in ['cas1', 'cas2']):
        return "adaptation"
    # Accessory
    if 'cas4' in hmm:
        return "accessory"
    # Default to backbone (structural components)
    return "backbone"


# HMM name to subtype mapping
def parse_subtype_from_hmm(hmm_name: str) -> tuple[str, str, str]:
    """
    Parse CRISPR-Cas type and subtype from HMM name.

    Returns: (type, subtype, protein_class)
        type: I, II, III, IV, V, VI, or Core
        subtype: A, B, C, D, E, F, or None
        protein_class: effector, backbone, adaptation, accessory
    """
    # Check TIGR mapping first (most reliable)
    if hmm_name in TIGR_TO_CRISPR:
        return TIGR_TO_CRISPR[hmm_name]

    hmm = hmm_name.lower()

    # Multi-type HMMs (e.g., cas2_I_II_III_V_maka) are CORE proteins, not type-specific
    # These are adaptation/accessory proteins found across multiple system types
    if re.search(r'cas[12]_[iv_]+_maka', hmm) or re.search(r'cas[12]_.*_[iv]_.*_[iv]', hmm):
        return ('Core', None, 'adaptation')
    if re.search(r'cas4_[iv_]+_maka', hmm) or re.search(r'cas4_.*_[iv]_.*_[iv]', hmm):
        return ('Core', None, 'accessory')

    # Check for "maka" naming convention from CRISPRCasFinder
    # Format: protein_SUBTYPE_maka_N (e.g., csm3_IIIAD_maka_5 → III-A)
    # But NOT multi-type patterns like cas2_I_II_III_V_maka
    maka_match = re.search(r'_([iv]+[a-d]?)_maka', hmm)
    if maka_match:
        subtype_code = maka_match.group(1).upper()
        # Skip if this looks like a multi-type pattern
        if '_' not in hmm.split('_maka')[0].split('_')[-2]:
            cas_type = _parse_roman_type(subtype_code)
            subtype = _parse_subtype_letter(subtype_code)
            protein_class = _infer_protein_class(hmm)
            return (cas_type, subtype, protein_class)

    # Check for numeric suffix (e.g., cas9_maka_4_II → Type II)
    maka_suffix = re.search(r'maka_\d+_([iv]+)$', hmm)
    if maka_suffix:
        type_code = maka_suffix.group(1).upper()
        cas_type = _parse_roman_type(type_code)
        protein_class = _infer_protein_class(hmm)
        return (cas_type, None, protein_class)

    # Core adaptation proteins (general)
    if any(x in hmm for x in ['cas1', 'cas2']):
        return ('Core', None, 'adaptation')

    # Type II - Cas9
    if 'cas9' in hmm:
        if 'iib' in hmm or '_ii_b' in hmm or 'ii-b' in hmm:
            return ('II', 'B', 'effector')
        if 'iic' in hmm or '_ii_c' in hmm or 'ii-c' in hmm:
            return ('II', 'C', 'effector')
        return ('II', None, 'effector')

    # Type V - Cas12 effectors (the actual signature proteins)
    if 'cas12' in hmm or 'cpf1' in hmm:
        if 'va' in hmm or '_v_a' in hmm or 'v-a' in hmm or '12a' in hmm:
            return ('V', 'A', 'effector')
        if 'vb' in hmm or '_v_b' in hmm or 'v-b' in hmm or '12b' in hmm:
            return ('V', 'B', 'effector')
        if 'vc' in hmm or '_v_c' in hmm or 'v-c' in hmm or '12c' in hmm:
            return ('V', 'C', 'effector')
        return ('V', None, 'effector')
    # C2c1 = Cas12b, C2c3 = Cas12c
    if hmm == 'c2c1':
        return ('V', 'B', 'effector')
    if hmm == 'c2c3':
        return ('V', 'C', 'effector')

    # Type VI - Cas13
    if 'cas13' in hmm:
        return ('VI', None, 'effector')

    # Type III - Cas10 and Csm/Cmr
    if 'cas10' in hmm:
        if 'iiia' in hmm:
            return ('III', 'A', 'effector')
        if 'iiib' in hmm:
            return ('III', 'B', 'effector')
        if 'iiic' in hmm:
            return ('III', 'C', 'effector')
        if 'iiid' in hmm:
            return ('III', 'D', 'effector')
        return ('III', None, 'effector')

    if 'csm' in hmm or 'iiia' in hmm:
        return ('III', 'A', 'backbone')

    if 'cmr' in hmm or 'iiib' in hmm:
        return ('III', 'B', 'backbone')

    if 'iiic' in hmm:
        return ('III', 'C', 'backbone')

    if 'iiid' in hmm or 'csx10' in hmm:
        return ('III', 'D', 'backbone')

    # Type IV - Csf
    if 'csf' in hmm or '_iv' in hmm:
        return ('IV', None, 'backbone')

    # Type I - Cas3 and variants
    if 'cas3' in hmm:
        return ('I', None, 'effector')

    # Type I subtypes from Cas5, Cas7, Cas8
    for cas in ['cas5', 'cas6', 'cas7', 'cas8', 'cse', 'csy', 'csc']:
        if cas in hmm:
            if '_ia' in hmm or 'i-a' in hmm or '_i_a' in hmm:
                return ('I', 'A', 'backbone')
            if '_ib' in hmm or 'i-b' in hmm or '_i_b' in hmm:
                return ('I', 'B', 'backbone')
            if '_ic' in hmm or 'i-c' in hmm or '_i_c' in hmm:
                return ('I', 'C', 'backbone')
            if '_id' in hmm or 'i-d' in hmm or '_i_d' in hmm:
                return ('I', 'D', 'backbone')
            if '_ie' in hmm or 'i-e' in hmm or '_i_e' in hmm:
                return ('I', 'E', 'backbone')
            if '_if' in hmm or 'i-f' in hmm or '_i_f' in hmm:
                return ('I', 'F', 'backbone')
            return ('I', None, 'backbone')

    # Cas4 - associated with I, II, V
    if 'cas4' in hmm:
        return ('Core', None, 'accessory')

    return ('Unknown', None, 'unknown')


def cluster_into_loci(
    results_df: pd.DataFrame,
    protein_coords: dict,
    max_gap_kb: int = 10,
) -> pd.DataFrame:
    """
    Cluster Cas proteins into loci based on genomic proximity.

    Proteins on the same contig within max_gap_kb are grouped into a locus.
    Each locus is assigned a type based on signature proteins.

    Returns DataFrame with locus assignments and locus-level classification.
    """
    # Add coordinates to results
    results_df = results_df.copy()
    results_df['contig'] = results_df['protein_id'].map(
        lambda x: protein_coords.get(x, (None, None, None))[0]
    )
    results_df['start'] = results_df['protein_id'].map(
        lambda x: protein_coords.get(x, (None, None, None))[1]
    )
    results_df['end'] = results_df['protein_id'].map(
        lambda x: protein_coords.get(x, (None, None, None))[2]
    )

    # Remove proteins without coordinates
    results_df = results_df.dropna(subset=['contig', 'start'])

    # Cluster by genome and contig
    loci = []
    locus_id = 0

    for (genome, contig), group in results_df.groupby(['genome', 'contig']):
        # Sort by position
        group = group.sort_values('start')

        current_locus = []
        prev_end = None

        for _, row in group.iterrows():
            if prev_end is None or (row['start'] - prev_end) <= max_gap_kb * 1000:
                current_locus.append(row)
            else:
                # Gap too large - save current locus and start new one
                if current_locus:
                    locus_id += 1
                    for r in current_locus:
                        loci.append({**r.to_dict(), 'locus_id': locus_id})
                current_locus = [row]
            prev_end = max(prev_end or 0, row['end'])

        # Save final locus
        if current_locus:
            locus_id += 1
            for r in current_locus:
                loci.append({**r.to_dict(), 'locus_id': locus_id})

    return pd.DataFrame(loci)


def classify_loci(loci_df: pd.DataFrame) -> pd.DataFrame:
    """
    Classify each locus by type based on signature proteins.

    Signature proteins:
    - Type I: Cas3 (effector)
    - Type II: Cas9 (effector)
    - Type III: Cas10 (effector)
    - Type IV: Csf proteins
    - Type V: Cas12/Cpf1 (effector)
    - Type VI: Cas13 (effector)
    """
    locus_summaries = []

    for locus_id, group in loci_df.groupby('locus_id'):
        genome = group['genome'].iloc[0]
        contig = group['contig'].iloc[0]
        start = group['start'].min()
        end = group['end'].max()
        n_proteins = len(group)
        near_array = group['near_array'].any()

        # Collect types and subtypes from proteins
        types = group['cas_type'].value_counts()
        subtypes = group[group['subtype'].notna()].groupby(['cas_type', 'subtype']).size()

        # Determine locus type by signature proteins (effectors)
        effectors = group[group['protein_class'] == 'effector']

        if len(effectors) > 0:
            # Use effector type
            locus_type = effectors['cas_type'].mode().iloc[0] if len(effectors) > 0 else 'Unknown'
            locus_subtype = None

            # Get most specific subtype from effectors
            effector_subtypes = effectors[effectors['subtype'].notna()]
            if len(effector_subtypes) > 0:
                locus_subtype = effector_subtypes['subtype'].mode().iloc[0]
        else:
            # No effector - use most common backbone type
            backbone = group[group['protein_class'] == 'backbone']
            if len(backbone) > 0:
                locus_type = backbone['cas_type'].mode().iloc[0]
                backbone_subtypes = backbone[backbone['subtype'].notna()]
                locus_subtype = backbone_subtypes['subtype'].mode().iloc[0] if len(backbone_subtypes) > 0 else None
            else:
                locus_type = 'Core'
                locus_subtype = None

        # Skip Core-only loci (just Cas1/Cas2 without a system)
        if locus_type == 'Core':
            # Check if there's any non-Core type
            non_core = types[types.index != 'Core']
            if len(non_core) > 0:
                locus_type = non_core.idxmax()

        locus_summaries.append({
            'locus_id': locus_id,
            'genome': genome,
            'contig': contig,
            'start': start,
            'end': end,
            'length_kb': (end - start) / 1000,
            'n_proteins': n_proteins,
            'locus_type': locus_type,
            'locus_subtype': locus_subtype,
            'near_array': near_array,
            'protein_ids': ','.join(group['protein_id'].tolist()),
            'hmm_names': ','.join(group['hmm_name'].unique().tolist()),
        })

    return pd.DataFrame(locus_summaries)


def get_predicates_for_type(cas_type: str, subtype: str | None, protein_class: str) -> list[str]:
    """Get predicates for a CRISPR-Cas classification."""
    predicates = ['crispr_associated']

    if cas_type == 'Core':
        if protein_class == 'adaptation':
            predicates.append('crispr_adaptation')
        else:
            predicates.append('crispr_accessory')
        return predicates

    if cas_type == 'Unknown':
        return predicates

    # Add class
    if cas_type in ('I', 'III', 'IV'):
        predicates.append('crispr_class1')
    elif cas_type in ('II', 'V', 'VI'):
        predicates.append('crispr_class2')

    # Add type
    type_pred = f'crispr_type_{cas_type.lower()}'
    predicates.append(type_pred)

    # Add subtype if known
    if subtype:
        subtype_pred = f'crispr_type_{cas_type.lower()}_{subtype.lower()}'
        predicates.append(subtype_pred)

    # Add effector/nuclease markers
    if protein_class == 'effector':
        predicates.append('cas_nuclease')

    return predicates


def classify_crispr_cas(
    db_path: str,
    results_dir: str | None = None,
    update_predicates: bool = True,
    array_distance_kb: int = 20,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Classify CRISPR-Cas proteins and update predicates.

    Args:
        db_path: Path to sharur.duckdb
        results_dir: Path to CRISPRCasFinder results (default: {db_dir}/crisprcasfinder_results)
        update_predicates: Whether to update predicates in database
        array_distance_kb: Max distance (kb) from CRISPR array for validation
        verbose: Print progress

    Returns:
        DataFrame with classification results
    """
    db = duckdb.connect(db_path)
    db_dir = Path(db_path).parent

    if results_dir is None:
        results_dir = db_dir / "crisprcasfinder_results"
    else:
        results_dir = Path(results_dir)

    results_file = results_dir / "CRISPRCasFinder_hits_df.tsv"
    if not results_file.exists():
        print(f"Error: CRISPRCasFinder results not found at {results_file}")
        print("Run: astra search --prot_in <proteins> --installed_hmms CRISPRCasFinder --outdir <output>")
        return pd.DataFrame()

    if verbose:
        print(f"Loading CRISPRCasFinder results from {results_file}")

    # Load results
    hits_df = pd.read_csv(results_file, sep='\t')

    if verbose:
        print(f"Found {len(hits_df)} hits for {hits_df['sequence_id'].nunique()} proteins")

    # Get protein->genome mapping
    protein_genome = dict(db.execute(
        "SELECT protein_id, bin_id FROM proteins"
    ).fetchall())

    # Get protein coordinates for array proximity check
    protein_coords = {}
    for row in db.execute("""
        SELECT protein_id, contig_id, start, end_coord FROM proteins
    """).fetchall():
        protein_coords[row[0]] = (row[1], row[2], row[3])

    # Get CRISPR array coordinates
    array_coords = []
    for row in db.execute("""
        SELECT contig_id, start, end_coord FROM loci WHERE locus_type = 'crispr_array'
    """).fetchall():
        array_coords.append((row[0], row[1], row[2]))

    if verbose:
        print(f"Found {len(array_coords)} CRISPR arrays for validation")

    # Classify each protein
    results = []
    for _, hit in hits_df.iterrows():
        protein_id = hit['sequence_id']
        hmm_name = hit['hmm_name']
        bitscore = hit['bitscore']
        evalue = hit['evalue']

        # Parse type from HMM
        cas_type, subtype, protein_class = parse_subtype_from_hmm(hmm_name)

        # Check proximity to CRISPR array
        near_array = False
        if protein_id in protein_coords:
            p_contig, p_start, p_end = protein_coords[protein_id]
            for a_contig, a_start, a_end in array_coords:
                if p_contig == a_contig:
                    distance = min(
                        abs(p_start - a_end),
                        abs(p_end - a_start)
                    )
                    if distance <= array_distance_kb * 1000:
                        near_array = True
                        break

        # Get genome
        genome = protein_genome.get(protein_id, 'unknown')

        results.append({
            'protein_id': protein_id,
            'genome': genome,
            'hmm_name': hmm_name,
            'cas_type': cas_type,
            'subtype': subtype,
            'protein_class': protein_class,
            'bitscore': bitscore,
            'evalue': evalue,
            'near_array': near_array,
        })

    results_df = pd.DataFrame(results)

    # Deduplicate - keep best hit per protein
    results_df = results_df.sort_values('bitscore', ascending=False)
    results_df = results_df.drop_duplicates(subset=['protein_id'], keep='first')

    # Cluster into loci
    if verbose:
        print("Clustering proteins into loci...")

    loci_df = cluster_into_loci(results_df, protein_coords, max_gap_kb=10)
    locus_summary = classify_loci(loci_df)

    # Print locus-level summary
    if verbose:
        print(f"\n=== CRISPR-Cas LOCUS CLASSIFICATION ===")

        # Complete systems (4+ proteins) are the meaningful unit
        complete = locus_summary[locus_summary['n_proteins'] >= 4]
        validated = locus_summary[locus_summary['near_array'] == True]
        complete_validated = complete[complete['near_array'] == True]

        print(f"\nOverview:")
        print(f"  Total protein clusters: {len(locus_summary)}")
        print(f"  Complete systems (≥4 proteins): {len(complete)}")
        print(f"  Array-validated loci: {len(validated)}")
        print(f"  Complete + validated: {len(complete_validated)}")

        # Complete systems by type
        print(f"\nComplete CRISPR-Cas systems (≥4 proteins):")
        if len(complete) > 0:
            type_counts = complete.groupby('locus_type').agg({
                'locus_id': 'count',
                'genome': 'nunique',
            }).rename(columns={'locus_id': 'n_loci'})
            for t, row in type_counts.iterrows():
                print(f"  Type {t}: {int(row['n_loci'])} systems in {int(row['genome'])} genomes")

        # Complete systems by subtype
        print(f"\nBy subtype (complete systems only):")
        if len(complete) > 0:
            complete_typed = complete[~complete['locus_type'].isin(['Core', 'Unknown'])]
            subtype_counts = complete_typed.groupby(['locus_type', 'locus_subtype']).size()
            for (t, s), n in sorted(subtype_counts.items()):
                label = f"{t}-{s}" if s else f"{t} (unsubtyped)"
                print(f"  {label}: {n}")

        # Array-validated systems
        print(f"\nArray-validated systems:")
        if len(validated) > 0:
            val_counts = validated.groupby('locus_type').size()
            for t, n in val_counts.items():
                print(f"  Type {t}: {n}")

        # Genome coverage
        print(f"\nGenome coverage:")
        genomes_with_complete = complete['genome'].nunique()
        genomes_with_any = locus_summary['genome'].nunique()
        print(f"  Genomes with complete systems: {genomes_with_complete}")
        print(f"  Genomes with any Cas proteins: {genomes_with_any}")

        # Note about fragmentation
        orphans = len(locus_summary[locus_summary['n_proteins'] == 1])
        print(f"\nNote: {orphans} isolated Cas proteins (MAG fragmentation or orphans)")

    # Update predicates
    if update_predicates:
        if verbose:
            print("\nUpdating predicates...")

        updated = 0
        for _, row in results_df.iterrows():
            protein_id = row['protein_id']
            cas_type = row['cas_type']
            subtype = row['subtype']
            protein_class = row['protein_class']

            # Get predicates for this classification
            new_predicates = get_predicates_for_type(cas_type, subtype, protein_class)

            # Add direct access predicate for subtype
            if cas_type not in ('Core', 'Unknown'):
                if subtype:
                    new_predicates.append(f"crispr_subtype:{cas_type}-{subtype}")
                else:
                    new_predicates.append(f"crispr_subtype:{cas_type}")

            # Get current predicates
            current = db.execute("""
                SELECT predicates FROM protein_predicates WHERE protein_id = ?
            """, [protein_id]).fetchone()

            if current:
                current_preds = set(current[0])
                current_preds.update(new_predicates)
                db.execute("""
                    UPDATE protein_predicates
                    SET predicates = ?, updated_at = CURRENT_TIMESTAMP
                    WHERE protein_id = ?
                """, [list(current_preds), protein_id])
                updated += 1

        db.commit()
        if verbose:
            print(f"Updated predicates for {updated} proteins")

    # Save detailed results
    protein_output = db_dir / "crispr_cas_proteins.tsv"
    results_df.to_csv(protein_output, sep='\t', index=False)

    locus_output = db_dir / "crispr_cas_loci.tsv"
    locus_summary.to_csv(locus_output, sep='\t', index=False)

    if verbose:
        print(f"\nProtein-level results: {protein_output}")
        print(f"Locus-level results: {locus_output}")

    db.close()
    return locus_summary


def main():
    parser = argparse.ArgumentParser(
        description="Classify CRISPR-Cas systems",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--db", required=True,
        help="Path to sharur.duckdb"
    )
    parser.add_argument(
        "--results-dir",
        help="Path to CRISPRCasFinder results (default: auto-detect)"
    )
    parser.add_argument(
        "--array-distance", type=int, default=20,
        help="Max distance (kb) from CRISPR array for validation (default: 20)"
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

    results = classify_crispr_cas(
        db_path=args.db,
        results_dir=args.results_dir,
        update_predicates=not args.no_update,
        array_distance_kb=args.array_distance,
        verbose=not args.quiet,
    )

    if results.empty:
        print("No CRISPR-Cas proteins classified")
        return 1

    return 0


if __name__ == "__main__":
    exit(main())
