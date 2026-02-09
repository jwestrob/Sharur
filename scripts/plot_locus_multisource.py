#!/usr/bin/env python3
"""
Generate publication-quality gene neighborhood diagrams with multi-source annotation.

Annotation priority: Foldseek structural > DefenseFinder > PADLOC > PFAM/KEGG/VOGdb

Usage:
    python plot_locus_multisource.py --db data/dataset/sharur.duckdb --protein PROTEIN_ID --window 10 --output figure.png
"""

import argparse
import re
from pathlib import Path
from typing import Optional, Dict, List, Tuple

import duckdb
import pandas as pd
from dna_features_viewer import GraphicFeature, GraphicRecord
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# Color scheme by annotation source
# Keys are canonical display names; _SOURCE_NORMALIZE maps DB values to these keys
COLORS = {
    'Foldseek': '#FF6B6B',       # Red - Foldseek structural homology
    'DefenseFinder': '#4ECDC4',  # Teal - DefenseFinder defense systems
    'PADLOC': '#95E1D3',         # Light teal - PADLOC defense systems
    'PFAM': '#3498DB',           # Blue - PFAM domains
    'KEGG': '#9B59B6',           # Purple - KEGG orthologs
    'VOGdb': '#F39C12',          # Orange - VOGdb viral genes
    'CAZy': '#2ECC71',           # Green - CAZy carbohydrate-active enzymes
    'HydDB': '#E74C3C',          # Red - HydDB hydrogenase classification
    'none': '#BDC3C7',           # Light gray - unannotated
}

# Map database source values (any case) to canonical display names
_SOURCE_NORMALIZE = {
    'pfam': 'PFAM',
    'kegg': 'KEGG',
    'vogdb': 'VOGdb',
    'cazy': 'CAZy',
    'hyddb': 'HydDB',
    'foldseek': 'Foldseek',
    'defensefinder': 'DefenseFinder',
    'padloc': 'PADLOC',
}


def normalize_source(source: str) -> str:
    """Normalize a source string to its canonical display name."""
    return _SOURCE_NORMALIZE.get(source.lower(), source)


def clean_pdb_description(desc: str) -> str:
    """Clean up PDB descriptions for use as gene labels."""
    if not desc:
        return "?"

    # Store original for debugging
    original = desc

    # 1. Remove PDB IDs at the start (many variations)
    # Examples: "6abc_A", "Type-assembly1_A-2", "6fiw-assembly1_A-2", "7F0N-ASSEMBLY1_A"
    desc = re.sub(r'^\w{4}-?assembly\d+_[A-Z\d]+(-\d+)?\s+', '', desc, flags=re.IGNORECASE)
    desc = re.sub(r'^\w{4}_[A-Z\d]+\s+', '', desc)
    desc = re.sub(r'^[\w\-]+assembly\d+_[A-Z\d]+-?\d*\s+', '', desc)

    # 2. Remove experimental method prefixes (case-insensitive)
    # Do this in multiple passes to catch nested patterns
    desc = re.sub(r'^x-?ray\s+(crystal\s+)?', '', desc, flags=re.IGNORECASE)  # Remove x-ray first
    experimental_methods = [
        r'crystal\s+structure\s+of\s+(the\s+)?',
        r'cryo-?em\s+structure\s+of\s+(the\s+)?',
        r'electron\s+microscopy\s+structure\s+of\s+(the\s+)?',
        r'nmr\s+(solution\s+)?structure\s+of\s+(the\s+)?',
        r'neutron\s+structure\s+of\s+(the\s+)?',
    ]
    for pattern in experimental_methods:
        desc = re.sub(pattern, '', desc, flags=re.IGNORECASE)

    # 3. Remove model/molecular descriptors
    desc = re.sub(r'^atomic\s+model\s+of\s+(the\s+)?', '', desc, flags=re.IGNORECASE)
    desc = re.sub(r'^molecular\s+model\s+of\s+(the\s+)?', '', desc, flags=re.IGNORECASE)
    desc = re.sub(r'^model\s+of\s+(the\s+)?', '', desc, flags=re.IGNORECASE)
    desc = re.sub(r'^structure\s+of\s+(the\s+)?', '', desc, flags=re.IGNORECASE)

    # 4. Remove articles at the start
    desc = re.sub(r'^(the|an?)\s+', '', desc, flags=re.IGNORECASE)

    # 5. Remove state/form descriptors
    desc = re.sub(r'^(apo|holo)\s+(form\s+of\s+)?', '', desc, flags=re.IGNORECASE)
    # Only remove DNA/RNA/ligand etc when they're explicitly bound/free/complexed
    desc = re.sub(r'^(ligand|substrate|inhibitor|product)-(bound|free|complexed)\s+', '', desc, flags=re.IGNORECASE)
    desc = re.sub(r'^(dna|rna)-(bound|free|complexed)\s+', '', desc, flags=re.IGNORECASE)

    # 6. Extract first meaningful phrase (stop at context clauses)
    # Stop at: "from [organism]", "in complex with", "bound to", "with", commas
    match = re.search(
        r'^([^,]+?)(?:\s+from\s+|\s+in\s+complex\s+with\s+|\s+complexed\s+with\s+|\s+bound\s+to\s+|\s+with\s+|\s+containing\s+|,|$)',
        desc,
        flags=re.IGNORECASE
    )
    if match:
        desc = match.group(1).strip()

    # 7. Clean up remaining articles after all other processing
    desc = re.sub(r'^(the|an?)\s+', '', desc, flags=re.IGNORECASE)

    # 8. Remove trailing state descriptors
    desc = re.sub(r'\s+(apo|holo)\s*$', '', desc, flags=re.IGNORECASE)

    # 9. Clean up whitespace
    desc = ' '.join(desc.split())

    # 10. Capitalize first letter
    if desc:
        desc = desc[0].upper() + desc[1:]

    # If we ended up with nothing or just punctuation, return "?"
    if not desc or not re.search(r'[a-zA-Z]', desc):
        return "?"

    return desc  # No character limit


def get_color_for_source(source: str) -> str:
    """Assign color based on annotation source."""
    canonical = normalize_source(source)
    return COLORS.get(canonical, COLORS['none'])


def load_defensefinder(db_path: Path) -> Dict[str, str]:
    """Load DefenseFinder annotations."""
    df_path = db_path.parent / "defensefinder_results" / "DefenseFinder_hits_df.tsv"
    if not df_path.exists():
        # Try alternate location
        df_path = db_path.parent / "annotations" / "defensefinder.tsv"

    if not df_path.exists():
        return {}

    df = pd.read_csv(df_path, sep='\t')
    annotations = {}

    for _, row in df.iterrows():
        protein_id = row['sequence_id']
        hmm_name = row.get('hmm_name', '')

        # Parse system__component format
        if '__' in hmm_name:
            parts = hmm_name.split('__')
            label = f"{parts[0]} {parts[-1]}"
        else:
            label = hmm_name.replace('_', ' ')

        # Only keep best hit per protein (lowest e-value, already sorted)
        if protein_id not in annotations:
            annotations[protein_id] = label

    return annotations


def load_padloc(db_path: Path) -> Dict[str, str]:
    """Load PADLOC annotations."""
    # Try multiple possible locations
    possible_paths = [
        db_path.parent / "padloc_results" / "padloc_hits_df.tsv",
        db_path.parent / "data" / db_path.parent.name / "padloc_results" / "padloc_hits_df.tsv",
        db_path.parent / "annotations" / "padloc.tsv",
    ]

    padloc_path = None
    for p in possible_paths:
        if p.exists():
            padloc_path = p
            break

    if not padloc_path:
        return {}

    df = pd.read_csv(padloc_path, sep='\t')
    annotations = {}

    for _, row in df.iterrows():
        protein_id = row['sequence_id']
        hmm = row.get('hmm_name', '')

        # Clean up PADLOC HMM names
        label = hmm.replace('_', ' ')

        # Only keep best hit per protein (lowest e-value, already sorted)
        if protein_id not in annotations:
            annotations[protein_id] = label

    return annotations


def load_kegg_descriptions() -> Dict[str, str]:
    """Load KEGG KO descriptions from reference file."""
    ko_list_path = Path(__file__).parent.parent / "data" / "reference" / "ko_list"

    if not ko_list_path.exists():
        return {}

    ko_map = {}
    with open(ko_list_path, 'r') as f:
        # Skip header
        next(f)
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 13:
                ko_id = parts[0]  # knum column
                simplified_def = parts[12]  # simplified_definition column
                if simplified_def:
                    ko_map[ko_id] = simplified_def

    return ko_map


def get_best_annotation(
    protein_id: str,
    db: duckdb.DuckDBPyConnection,
    defensefinder: Dict[str, str],
    padloc: Dict[str, str],
    kegg_map: Dict[str, str],
) -> Tuple[str, str]:
    """
    Get best annotation for a protein using priority stack.

    Priority: Foldseek > DefenseFinder > PADLOC > PFAM/KEGG/VOGdb

    Returns: (label, source)
    """
    # 1. Foldseek structural homology (table may not exist)
    try:
        foldseek = db.execute(
            """SELECT pdb_id_chain, pdb_description, evalue
               FROM foldseek_hits
               WHERE protein_id = ? AND evalue < 0.01
               ORDER BY evalue LIMIT 1""",
            [protein_id]
        ).fetchone()

        if foldseek:
            desc = foldseek[1] or foldseek[0]
            return (clean_pdb_description(desc), 'Foldseek')
    except duckdb.CatalogException:
        pass  # foldseek_hits table doesn't exist

    # 2. DefenseFinder
    if protein_id in defensefinder:
        return (defensefinder[protein_id], 'DefenseFinder')

    # 3. PADLOC
    if protein_id in padloc:
        return (padloc[protein_id], 'PADLOC')

    # 4. PFAM/KEGG/VOGdb (best by e-value)
    seq_ann = db.execute(
        """SELECT name, source FROM annotations
           WHERE protein_id = ?
           ORDER BY evalue NULLS LAST LIMIT 1""",
        [protein_id]
    ).fetchone()

    if seq_ann and seq_ann[0] and seq_ann[0] != 'NO HITS':
        name = seq_ann[0].split('(')[0].strip()
        source = normalize_source(seq_ann[1])

        # Translate KEGG KO IDs to descriptions
        if source == 'KEGG' and name in kegg_map:
            name = kegg_map[name]

        return (name, source)

    return ('?', 'none')


def plot_locus(
    db_path: Path,
    protein_id: str,
    window: int = 10,
    output_path: Optional[Path] = None,
    title: Optional[str] = None,
    subtitle: Optional[str] = None,
):
    """
    Generate gene neighborhood diagram with multi-source annotation.

    Args:
        db_path: Path to sharur.duckdb
        protein_id: Center protein ID
        window: Genes on each side
        output_path: Output PNG path
        title: Custom title (default: auto-generated)
        subtitle: Custom subtitle (default: shows annotation sources)
    """
    db = duckdb.connect(str(db_path), read_only=True)

    # Load external annotation sources
    print("Loading DefenseFinder annotations...")
    defensefinder = load_defensefinder(db_path)
    print(f"  Found {len(defensefinder)} DefenseFinder hits")

    print("Loading PADLOC annotations...")
    padloc = load_padloc(db_path)
    print(f"  Found {len(padloc)} PADLOC hits")

    print("Loading KEGG KO descriptions...")
    kegg_map = load_kegg_descriptions()
    print(f"  Loaded {len(kegg_map)} KO descriptions")

    # Get anchor protein
    anchor = db.execute(
        """SELECT contig_id, gene_index, start, end_coord, strand
           FROM proteins WHERE protein_id = ?""",
        [protein_id]
    ).fetchone()

    if not anchor:
        raise ValueError(f"Protein {protein_id} not found")

    contig_id, gene_index, anchor_start, anchor_end, anchor_strand = anchor

    # Get neighborhood
    rows = db.execute(
        """SELECT protein_id, start, end_coord, strand, gene_index
           FROM proteins
           WHERE contig_id = ? AND gene_index BETWEEN ? AND ?
           ORDER BY gene_index""",
        [contig_id, gene_index - window, gene_index + window]
    ).fetchall()

    if not rows:
        raise ValueError(f"No neighborhood found for {protein_id}")

    # Calculate span
    min_coord = min(r[1] for r in rows)
    max_coord = max(r[2] for r in rows)

    # Check for CRISPR arrays in region
    crispr_arrays = db.execute(
        """SELECT locus_id, start, end_coord FROM loci
           WHERE contig_id = ? AND locus_type = 'crispr_array'
           AND start <= ? AND end_coord >= ?""",
        [contig_id, max_coord + 1000, min_coord - 1000]
    ).fetchall()

    # Extend to include arrays
    for _, arr_start, arr_end in crispr_arrays:
        min_coord = min(min_coord, arr_start)
        max_coord = max(max_coord, arr_end)

    # Build features
    features = []
    gene_positions = []  # For gene number labels

    # Add CRISPR arrays
    for locus_id, arr_start, arr_end in crispr_arrays:
        # Extract repeat/spacer count from locus metadata if available
        metadata = db.execute(
            "SELECT metadata FROM loci WHERE locus_id = ?", [locus_id]
        ).fetchone()

        repeats = '?'
        spacers = '?'
        if metadata and metadata[0]:
            import json
            try:
                meta = json.loads(metadata[0])
                repeats = meta.get('num_repeats', '?')
                spacers = meta.get('num_spacers', '?')
            except:
                pass

        features.append(
            GraphicFeature(
                start=arr_start - min_coord,
                end=arr_end - min_coord,
                strand=0,
                color='#FFE6E6',
                linecolor='#FF6B6B',
                linewidth=2,
                label=f'CRISPR\n({repeats}R/{spacers}S)',
            )
        )

    # Add genes
    sources_used = set()
    for pid, start, end, strand, gidx in rows:
        label, source = get_best_annotation(pid, db, defensefinder, padloc, kegg_map)
        sources_used.add(source)

        strand_int = 1 if strand == '+' else -1
        color = get_color_for_source(source)

        # Highlight query protein with star marker and darker red color
        if pid == protein_id:
            color = '#C0392B'  # Dark red to distinguish from Foldseek red
            # Add star above label to mark query protein
            label = f"★\n{label}"

        features.append(
            GraphicFeature(
                start=start - min_coord,
                end=end - min_coord,
                strand=strand_int,
                color=color,
                label=label,
            )
        )

        # Store gene position for labeling
        mid = (start + end) / 2 - min_coord
        gene_positions.append((gidx, mid))

    # Create plot
    record = GraphicRecord(sequence_length=max_coord - min_coord, features=features)

    fig, ax = plt.subplots(1, 1, figsize=(16, 6))
    record.plot(ax=ax, with_ruler=True, annotate_inline=False)

    # Add absolute genome coordinates
    ticks = ax.get_xticks()
    ax.set_xticklabels([f'{int(t + min_coord):,}' for t in ticks])

    # Add gene numbers BELOW the track (move down more to avoid overlap)
    # Use 1-indexed gene numbers for display (gene_index + 1)
    y_bottom = ax.get_ylim()[0]
    for gidx, mid in gene_positions:
        # Vertical line from bottom of track to label
        ax.plot([mid, mid], [y_bottom, y_bottom - 0.5],
                color='gray', linewidth=0.5, alpha=0.5, zorder=0)
        # Gene number label positioned lower (1-indexed for display)
        ax.text(mid, y_bottom - 0.7, f'g{gidx + 1}',
                ha='center', va='top', fontsize=7, color='gray')

    # Extend y-axis to accommodate gene labels
    ylim = ax.get_ylim()
    ax.set_ylim(ylim[0] - 1.2, ylim[1])

    # Title and subtitle
    if title is None:
        # Auto-generate title from context
        crispr_desc = f" — {len(crispr_arrays)} CRISPR array(s)" if crispr_arrays else ""
        title = f"Locus around {protein_id}{crispr_desc}"

    if subtitle is None:
        # Show full annotation priority order (always show all sources)
        subtitle = "Labels: Foldseek structural > DefenseFinder > PADLOC > PFAM/KEGG/VOGdb"

    ax.set_title(f"{title}\n{subtitle}", fontsize=11, pad=15)

    # Add color legend — only show sources actually used in this figure
    from matplotlib.patches import Rectangle

    # Canonical display order for legend
    _LEGEND_ORDER = ['Foldseek', 'DefenseFinder', 'PADLOC', 'PFAM', 'KEGG',
                     'VOGdb', 'CAZy', 'HydDB', 'none']
    # Normalize sources_used to canonical names
    sources_canonical = {normalize_source(s) for s in sources_used}

    legend_elements = []
    for src in _LEGEND_ORDER:
        if src in sources_canonical:
            display = 'Unannotated' if src == 'none' else src
            legend_elements.append((display, COLORS[src]))

    handles = []
    for label, color in legend_elements:
        handles.append(Rectangle((0, 0), 1, 1, fc=color, edgecolor='black', linewidth=0.5))

    ncol = min(len(legend_elements), 7)
    ax.legend(handles, [label for label, _ in legend_elements],
              loc='upper right', fontsize=8, ncol=ncol, framealpha=0.95,
              title='Annotation Source', title_fontsize=8)

    # Save
    if output_path is None:
        output_path = Path(f"locus_{protein_id}.png")

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

    print(f"\nSaved: {output_path}")
    print(f"Genes: {len(rows)}, CRISPR arrays: {len(crispr_arrays)}")
    print(f"Sources used: {', '.join(sorted(sources_used))}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--db', required=True, help='Path to sharur.duckdb')
    parser.add_argument('--protein', required=True, help='Center protein ID')
    parser.add_argument('--window', type=int, default=10, help='Genes on each side (default: 10)')
    parser.add_argument('--output', help='Output PNG path')
    parser.add_argument('--title', help='Custom title')
    parser.add_argument('--subtitle', help='Custom subtitle')

    args = parser.parse_args()

    plot_locus(
        db_path=Path(args.db),
        protein_id=args.protein,
        window=args.window,
        output_path=Path(args.output) if args.output else None,
        title=args.title,
        subtitle=args.subtitle,
    )


if __name__ == '__main__':
    main()
