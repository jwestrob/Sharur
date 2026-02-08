#!/usr/bin/env python3
"""
Direct protein FASTA ingestion for Bennu.

For datasets where proteins are already called (no Prodigal needed).
Parses Prodigal-style FASTA headers with coordinates.

Example header:
>Complete_Asgard_Heimdall_2022_curated_1|Complete_Asgard_Thorarchaeota_Heimdall_2022_curated|Heimdall_2022 # 59 # 1510 # 1
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path
from typing import Iterator, Tuple

import duckdb
import typer
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn

console = Console()


def parse_protein_fasta(fasta_path: Path) -> Iterator[dict]:
    """
    Parse a protein FASTA with Prodigal-style headers.

    Expected header format:
    >{protein_id}|{contig_id}|{genome_id} # {start} # {end} # {strand}

    Also handles simpler formats like:
    >{protein_id} # {start} # {end} # {strand} # ID={id};...
    """
    open_fn = gzip.open if str(fasta_path).endswith('.gz') else open

    header = None
    seq_parts = []

    with open_fn(fasta_path, 'rt') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield _parse_header_and_seq(header, ''.join(seq_parts))
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)

        if header is not None:
            yield _parse_header_and_seq(header, ''.join(seq_parts))


def _parse_header_and_seq(header: str, sequence: str) -> dict:
    """Parse header to extract protein metadata."""
    # Remove trailing stop codon if present
    if sequence.endswith('*'):
        sequence = sequence[:-1]

    # Try to parse Prodigal-style: ID|contig|genome # start # end # strand
    # Or: ID # start # end # strand # ID=...
    parts = header.split(' # ')

    if len(parts) >= 4:
        # Has coordinate info
        id_part = parts[0]
        start = int(parts[1])
        end = int(parts[2])
        strand = '+' if parts[3] in ('1', '+') else '-'

        # Parse ID part - may contain pipe separators
        if '|' in id_part:
            id_parts = id_part.split('|')
            protein_id = id_parts[0]
            contig_id = id_parts[1] if len(id_parts) > 1 else protein_id.rsplit('_', 1)[0]
            bin_id = id_parts[2] if len(id_parts) > 2 else contig_id
        else:
            protein_id = id_part.strip()
            # Try to derive contig from protein ID
            contig_id = protein_id.rsplit('_', 1)[0]
            bin_id = contig_id
    else:
        # Simple header without coordinates
        protein_id = header.split()[0]
        contig_id = protein_id.rsplit('_', 1)[0]
        bin_id = contig_id
        start = 1
        end = len(sequence) * 3
        strand = '+'

    return {
        'protein_id': protein_id,
        'contig_id': contig_id,
        'bin_id': bin_id,
        'start': start,
        'end': end,
        'strand': strand,
        'sequence': sequence,
        'sequence_length': len(sequence),
    }


def ingest_proteins(
    fasta_path: Path,
    db_path: Path,
    genome_name: str = None,
    force: bool = False,
) -> dict:
    """
    Ingest protein FASTA into Bennu DuckDB.

    Returns stats dict.
    """
    from bennu.storage.schema import SCHEMA

    if db_path.exists() and not force:
        raise FileExistsError(f"Database exists: {db_path}. Use --force to overwrite.")

    if db_path.exists() and force:
        db_path.unlink()

    db_path.parent.mkdir(parents=True, exist_ok=True)

    conn = duckdb.connect(str(db_path))
    conn.execute(SCHEMA)

    # Collect all proteins first
    proteins = list(parse_protein_fasta(fasta_path))

    # Identify unique contigs and bins
    contigs = {}
    bins = {}

    for p in proteins:
        cid = p['contig_id']
        bid = p['bin_id']

        if cid not in contigs:
            contigs[cid] = {'bin_id': bid, 'proteins': []}
        contigs[cid]['proteins'].append(p)

        if bid not in bins:
            bins[bid] = {'contigs': set()}
        bins[bid]['contigs'].add(cid)

    # Insert bins
    bin_name = genome_name or list(bins.keys())[0]
    for bid in bins:
        conn.execute(
            """
            INSERT INTO bins (bin_id, completeness, contamination, taxonomy, n_contigs, total_length)
            VALUES (?, 95.0, 2.0, ?, ?, ?)
            """,
            [bid, f"Asgard archaeon ({bin_name})", len(bins[bid]['contigs']), 0]
        )

    # Insert contigs
    for cid, cdata in contigs.items():
        # Calculate contig length from max protein end
        max_end = max(p['end'] for p in cdata['proteins'])
        conn.execute(
            """
            INSERT INTO contigs (contig_id, bin_id, length, gc_content)
            VALUES (?, ?, ?, NULL)
            """,
            [cid, cdata['bin_id'], max_end]
        )

    # Insert proteins with gene_index
    for cid, cdata in contigs.items():
        # Sort proteins by start position
        sorted_prots = sorted(cdata['proteins'], key=lambda x: x['start'])
        for idx, p in enumerate(sorted_prots):
            conn.execute(
                """
                INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand, gene_index, sequence, sequence_length)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                [p['protein_id'], p['contig_id'], p['bin_id'], p['start'], p['end'], p['strand'], idx, p['sequence'], p['sequence_length']]
            )

    # Update bin total_length
    for bid in bins:
        total_len = conn.execute(
            "SELECT SUM(length) FROM contigs WHERE bin_id = ?", [bid]
        ).fetchone()[0]
        conn.execute("UPDATE bins SET total_length = ? WHERE bin_id = ?", [total_len, bid])

    conn.close()

    return {
        'proteins': len(proteins),
        'contigs': len(contigs),
        'bins': len(bins),
        'db_path': str(db_path),
    }


app = typer.Typer(no_args_is_help=True, add_completion=False)


@app.command()
def main(
    fasta: Path = typer.Argument(..., help="Protein FASTA file (.faa or .faa.gz)"),
    output: Path = typer.Option(Path("data/bennu.duckdb"), "--output", "-o", help="Output DuckDB path"),
    genome_name: str = typer.Option(None, "--name", "-n", help="Genome name (default: from FASTA)"),
    force: bool = typer.Option(False, "--force", "-f", help="Overwrite existing database"),
):
    """
    Ingest a protein FASTA directly into Bennu DuckDB.

    For use with pre-called proteins (no Prodigal step needed).
    Annotations must be added separately.
    """
    if not fasta.exists():
        console.print(f"[red]File not found: {fasta}[/red]")
        raise typer.Exit(1)

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console,
    ) as progress:
        progress.add_task("Ingesting proteins...", total=None)
        stats = ingest_proteins(fasta, output, genome_name, force)

    console.print(f"[green]✓ Ingested {stats['proteins']} proteins[/green]")
    console.print(f"  Contigs: {stats['contigs']}")
    console.print(f"  Bins: {stats['bins']}")
    console.print(f"  Database: {stats['db_path']}")


if __name__ == "__main__":
    app()
