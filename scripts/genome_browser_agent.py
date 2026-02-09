#!/usr/bin/env python3
"""
Genome Browser Agent - Systematically examines genome in quarters.

This script provides utilities for an agent to browse through the genome
25% at a time, recording interesting findings to a persistent file.

Usage:
    python scripts/genome_browser_agent.py --quarter 1  # Browse first 25%
    python scripts/genome_browser_agent.py --quarter 2  # Browse second 25%
    python scripts/genome_browser_agent.py --quarter 3  # Browse third 25%
    python scripts/genome_browser_agent.py --quarter 4  # Browse fourth 25%
    python scripts/genome_browser_agent.py --summary    # Generate summary
"""

import argparse
import json
import duckdb
from pathlib import Path
from datetime import datetime

DB_PATH = "data/thorarchaeota_production/sharur.duckdb"
FINDINGS_PATH = "data/thorarchaeota_production/genome_browser_findings.jsonl"
SUMMARY_PATH = "data/thorarchaeota_production/genome_browser_summary.md"


def get_quarter_proteins(quarter: int) -> list:
    """Get protein IDs for a specific quarter (1-4)."""
    db = duckdb.connect(DB_PATH, read_only=True)

    # Get all proteins ordered by gene_index
    all_proteins = db.execute("""
        SELECT protein_id, gene_index, sequence_length
        FROM proteins
        ORDER BY gene_index
    """).fetchall()

    db.close()

    total = len(all_proteins)
    chunk_size = total // 4

    start_idx = (quarter - 1) * chunk_size
    end_idx = start_idx + chunk_size if quarter < 4 else total

    return all_proteins[start_idx:end_idx], start_idx, end_idx, total


def get_protein_details(protein_id: str) -> dict:
    """Get full details for a protein including annotations and neighbors."""
    db = duckdb.connect(DB_PATH, read_only=True)

    # Basic info
    info = db.execute("""
        SELECT protein_id, contig_id, gene_index, sequence_length,
               start, end_coord, strand
        FROM proteins WHERE protein_id = ?
    """, [protein_id]).fetchone()

    if not info:
        db.close()
        return None

    result = {
        'protein_id': info[0],
        'contig_id': info[1],
        'gene_index': info[2],
        'length': info[3],
        'start': info[4],
        'end': info[5],
        'strand': info[6],
        'annotations': [],
        'neighbors': []
    }

    # Annotations
    annots = db.execute("""
        SELECT source, annotation_id, description, evalue
        FROM annotations WHERE protein_id = ?
    """, [protein_id]).fetchall()

    for a in annots:
        result['annotations'].append({
            'source': a[0],
            'id': a[1],
            'description': a[2],
            'evalue': a[3]
        })

    # Neighbors (5 upstream, 5 downstream)
    neighbors = db.execute("""
        SELECT p.protein_id, p.gene_index, p.sequence_length,
               GROUP_CONCAT(a.description, '; ') as annotations
        FROM proteins p
        LEFT JOIN annotations a ON p.protein_id = a.protein_id
        WHERE p.contig_id = ? AND p.gene_index BETWEEN ? AND ?
        GROUP BY p.protein_id, p.gene_index, p.sequence_length
        ORDER BY p.gene_index
    """, [result['contig_id'], result['gene_index'] - 5, result['gene_index'] + 5]).fetchall()

    for n in neighbors:
        result['neighbors'].append({
            'protein_id': n[0],
            'gene_index': n[1],
            'length': n[2],
            'annotations': n[3] or 'unannotated'
        })

    db.close()
    return result


def get_quarter_summary(quarter: int) -> dict:
    """Get summary statistics for a quarter."""
    proteins, start_idx, end_idx, total = get_quarter_proteins(quarter)

    db = duckdb.connect(DB_PATH, read_only=True)

    protein_ids = [p[0] for p in proteins]
    placeholders = ','.join(['?'] * len(protein_ids))

    # Count annotations
    annotated = db.execute(f"""
        SELECT COUNT(DISTINCT protein_id)
        FROM annotations
        WHERE protein_id IN ({placeholders})
    """, protein_ids).fetchone()[0]

    # Size distribution
    sizes = db.execute(f"""
        SELECT
            SUM(CASE WHEN sequence_length >= 1000 THEN 1 ELSE 0 END) as giant,
            SUM(CASE WHEN sequence_length >= 500 AND sequence_length < 1000 THEN 1 ELSE 0 END) as large,
            SUM(CASE WHEN sequence_length < 500 THEN 1 ELSE 0 END) as small
        FROM proteins
        WHERE protein_id IN ({placeholders})
    """, protein_ids).fetchone()

    # Top domains
    top_domains = db.execute(f"""
        SELECT annotation_id, COUNT(*) as cnt
        FROM annotations
        WHERE protein_id IN ({placeholders}) AND source = 'pfam'
        GROUP BY annotation_id
        ORDER BY cnt DESC
        LIMIT 10
    """, protein_ids).fetchall()

    db.close()

    return {
        'quarter': quarter,
        'start_gene': start_idx,
        'end_gene': end_idx,
        'total_proteins': len(proteins),
        'annotated': annotated,
        'unannotated': len(proteins) - annotated,
        'giant': sizes[0],
        'large': sizes[1],
        'small': sizes[2],
        'top_domains': [(d[0], d[1]) for d in top_domains]
    }


def record_finding(finding: dict):
    """Append a finding to the findings file."""
    finding['timestamp'] = datetime.now().isoformat()

    with open(FINDINGS_PATH, 'a') as f:
        f.write(json.dumps(finding) + '\n')

    print(f"Recorded finding: {finding.get('title', 'Untitled')}")


def get_all_findings() -> list:
    """Load all recorded findings."""
    findings = []
    if Path(FINDINGS_PATH).exists():
        with open(FINDINGS_PATH) as f:
            for line in f:
                if line.strip():
                    findings.append(json.loads(line))
    return findings


def browse_quarter(quarter: int):
    """Interactive browse of a quarter - outputs info for agent to analyze."""
    proteins, start_idx, end_idx, total = get_quarter_proteins(quarter)
    summary = get_quarter_summary(quarter)

    print(f"\n{'='*60}")
    print(f"GENOME BROWSER - QUARTER {quarter}/4")
    print(f"{'='*60}")
    print(f"Genes {start_idx}-{end_idx} of {total}")
    print(f"Proteins in this quarter: {len(proteins)}")
    print(f"Annotated: {summary['annotated']} ({summary['annotated']/len(proteins)*100:.1f}%)")
    print(f"Unannotated: {summary['unannotated']} ({summary['unannotated']/len(proteins)*100:.1f}%)")
    print(f"\nSize distribution:")
    print(f"  Giant (≥1000aa): {summary['giant']}")
    print(f"  Large (500-999aa): {summary['large']}")
    print(f"  Small (<500aa): {summary['small']}")
    print(f"\nTop PFAM domains in this quarter:")
    for domain, count in summary['top_domains'][:10]:
        print(f"  {domain}: {count}")

    # Identify interesting proteins
    print(f"\n{'='*60}")
    print("NOTABLE PROTEINS IN THIS QUARTER")
    print(f"{'='*60}")

    db = duckdb.connect(DB_PATH, read_only=True)
    protein_ids = [p[0] for p in proteins]
    placeholders = ','.join(['?'] * len(protein_ids))

    # Giant unannotated
    giants = db.execute(f"""
        SELECT p.protein_id, p.sequence_length, p.gene_index
        FROM proteins p
        WHERE p.protein_id IN ({placeholders})
        AND p.sequence_length >= 1000
        AND NOT EXISTS (SELECT 1 FROM annotations a WHERE a.protein_id = p.protein_id)
        ORDER BY p.sequence_length DESC
        LIMIT 5
    """, protein_ids).fetchall()

    if giants:
        print("\n--- Giant Unannotated Proteins ---")
        for g in giants:
            print(f"  {g[0]}: {g[1]} aa (gene #{g[2]})")

    # Interesting annotations
    interesting = db.execute(f"""
        SELECT p.protein_id, p.sequence_length, a.description
        FROM proteins p
        JOIN annotations a ON p.protein_id = a.protein_id
        WHERE p.protein_id IN ({placeholders})
        AND (a.description LIKE '%secretion%'
             OR a.description LIKE '%toxin%'
             OR a.description LIKE '%defense%'
             OR a.description LIKE '%CRISPR%'
             OR a.description LIKE '%transpos%'
             OR a.description LIKE '%phage%'
             OR a.description LIKE '%virus%'
             OR a.description LIKE '%flagell%'
             OR a.description LIKE '%pilus%'
             OR a.description LIKE '%adhesin%')
        LIMIT 10
    """, protein_ids).fetchall()

    if interesting:
        print("\n--- Interesting Annotated Proteins ---")
        for i in interesting:
            print(f"  {i[0]} ({i[1]}aa): {i[2][:60]}")

    # Clusters of unannotated
    print("\n--- Regions with Consecutive Unannotated Proteins ---")

    # Find clusters
    unannotated_genes = db.execute(f"""
        SELECT p.gene_index
        FROM proteins p
        WHERE p.protein_id IN ({placeholders})
        AND NOT EXISTS (SELECT 1 FROM annotations a WHERE a.protein_id = p.protein_id)
        ORDER BY p.gene_index
    """, protein_ids).fetchall()

    unannotated_indices = [u[0] for u in unannotated_genes]

    # Find consecutive runs
    if unannotated_indices:
        clusters = []
        current_cluster = [unannotated_indices[0]]
        for idx in unannotated_indices[1:]:
            if idx == current_cluster[-1] + 1:
                current_cluster.append(idx)
            else:
                if len(current_cluster) >= 5:
                    clusters.append(current_cluster)
                current_cluster = [idx]
        if len(current_cluster) >= 5:
            clusters.append(current_cluster)

        # Show top clusters
        clusters.sort(key=len, reverse=True)
        for cluster in clusters[:5]:
            print(f"  Genes {cluster[0]}-{cluster[-1]}: {len(cluster)} consecutive unannotated")

    db.close()

    print(f"\n{'='*60}")
    print("SAMPLE PROTEINS FOR DETAILED EXAMINATION")
    print(f"{'='*60}")

    # Sample 10 proteins from different parts of the quarter
    sample_indices = [0, len(proteins)//10, len(proteins)//5, len(proteins)//3,
                      len(proteins)//2, 2*len(proteins)//3, 3*len(proteins)//4,
                      4*len(proteins)//5, 9*len(proteins)//10, len(proteins)-1]

    for idx in sample_indices[:10]:
        if idx < len(proteins):
            pid = proteins[idx][0]
            details = get_protein_details(pid)
            if details:
                annot_str = ', '.join([str(a['id']) for a in details['annotations'][:3]]) or 'unannotated'
                print(f"\n[Gene #{details['gene_index']}] {pid}")
                print(f"  Length: {details['length']} aa | Annotations: {annot_str}")
                print(f"  Neighborhood:")
                for n in details['neighbors']:
                    marker = ">>>" if n['protein_id'] == pid else "   "
                    annot = n['annotations'][:40] if n['annotations'] != 'unannotated' else 'unannotated'
                    print(f"    {marker} [{n['gene_index']}] {n['length']}aa: {annot}")


def generate_summary():
    """Generate a summary of all findings."""
    findings = get_all_findings()

    summary = f"""# Genome Browser Findings Summary

**Generated:** {datetime.now().strftime("%Y-%m-%d %H:%M")}
**Total Findings:** {len(findings)}

---

"""

    # Group by quarter
    by_quarter = {}
    for f in findings:
        q = f.get('quarter', 'unknown')
        if q not in by_quarter:
            by_quarter[q] = []
        by_quarter[q].append(f)

    for quarter in sorted(by_quarter.keys()):
        quarter_findings = by_quarter[quarter]
        summary += f"## Quarter {quarter}\n\n"

        for f in quarter_findings:
            summary += f"### {f.get('title', 'Untitled')}\n\n"
            if f.get('proteins'):
                summary += f"**Proteins:** {', '.join(f['proteins'][:5])}"
                if len(f.get('proteins', [])) > 5:
                    summary += f" (+{len(f['proteins'])-5} more)"
                summary += "\n\n"
            if f.get('description'):
                summary += f"{f['description']}\n\n"
            if f.get('significance'):
                summary += f"**Significance:** {f['significance']}\n\n"
            summary += "---\n\n"

    with open(SUMMARY_PATH, 'w') as f:
        f.write(summary)

    print(f"Summary written to: {SUMMARY_PATH}")
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Genome Browser Agent')
    parser.add_argument('--quarter', type=int, choices=[1, 2, 3, 4], help='Browse specific quarter')
    parser.add_argument('--summary', action='store_true', help='Generate findings summary')
    parser.add_argument('--protein', type=str, help='Get details for specific protein')

    args = parser.parse_args()

    if args.quarter:
        browse_quarter(args.quarter)
    elif args.summary:
        generate_summary()
    elif args.protein:
        details = get_protein_details(args.protein)
        print(json.dumps(details, indent=2))
    else:
        parser.print_help()
