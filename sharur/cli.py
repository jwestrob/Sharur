#!/usr/bin/env python3
"""
Sharur CLI - Metagenomic dataset exploration.

Commands:
- overview: Dataset summary
- genomes: List/filter genomes
- proteins: List/filter proteins
- neighborhood: View genomic context around a protein
- search: Search proteins by predicates

Example usage:
    sharur overview --db data/sharur.duckdb
    sharur genomes --taxonomy Archaea --db data/sharur.duckdb
    sharur proteins --min-length 2000 --db data/sharur.duckdb
    sharur neighborhood PROTEIN_ID --window 10 --db data/sharur.duckdb
    sharur search --has giant,unannotated --db data/sharur.duckdb
"""

from pathlib import Path
import sys
from typing import Optional

import typer

from sharur.operators import Sharur

# Default DB path
DEFAULT_DB = "data/sharur.duckdb"

# Main app with subcommands
app = typer.Typer(
    no_args_is_help=True,
    add_completion=False,
    help="Sharur - Metagenomic dataset exploration CLI",
    rich_markup_mode=None,  # Disable rich to work around Typer 0.15 bug
)


# ------------------------------------------------------------------ #
# Overview command
# ------------------------------------------------------------------ #


@app.command()
def overview(
    db: str = typer.Option(DEFAULT_DB, "--db", "-d", help="Path to DuckDB database"),
):
    """
    Show dataset overview with summary statistics.

    Displays genome/protein counts, annotation coverage,
    taxonomy distribution, and predicate summary.
    """
    db_path = Path(db)
    if not db_path.exists():
        typer.echo(f"DB not found: {db}", err=True)
        raise typer.Exit(1)

    b = Sharur(db_path)
    result = b.overview()
    typer.echo(result.data)


# ------------------------------------------------------------------ #
# Genomes command
# ------------------------------------------------------------------ #


@app.command()
def genomes(
    taxonomy: Optional[str] = typer.Option(None, "--taxonomy", "-t", help="Filter by taxonomy substring"),
    min_completeness: Optional[float] = typer.Option(None, "--min-comp", help="Minimum completeness %"),
    max_contamination: Optional[float] = typer.Option(None, "--max-contam", help="Maximum contamination %"),
    limit: int = typer.Option(20, "--limit", "-n", help="Maximum results"),
    db: str = typer.Option(DEFAULT_DB, "--db", "-d", help="Path to DuckDB database"),
):
    """
    List genomes (MAGs) with optional filtering.

    Examples:
        sharur genomes --taxonomy Archaea
        sharur genomes --min-comp 90 --max-contam 5
    """
    db_path = Path(db)
    if not db_path.exists():
        typer.echo(f"DB not found: {db}", err=True)
        raise typer.Exit(1)

    b = Sharur(db_path)
    result = b.list_genomes(
        taxonomy_filter=taxonomy,
        min_completeness=min_completeness,
        max_contamination=max_contamination,
        limit=limit,
    )
    typer.echo(result.data)

    if result.meta.truncated:
        typer.echo(f"\n[Showing {result.meta.rows} of {result.meta.total_rows} results]", err=True)


# ------------------------------------------------------------------ #
# Proteins command
# ------------------------------------------------------------------ #


@app.command()
def proteins(
    genome: Optional[str] = typer.Option(None, "--genome", "-g", help="Filter by genome (bin_id)"),
    contig: Optional[str] = typer.Option(None, "--contig", "-c", help="Filter by contig"),
    min_length: Optional[int] = typer.Option(None, "--min-length", help="Minimum length (aa)"),
    max_length: Optional[int] = typer.Option(None, "--max-length", help="Maximum length (aa)"),
    annotated: Optional[bool] = typer.Option(None, "--annotated/--unannotated", help="Filter by annotation status"),
    limit: int = typer.Option(50, "--limit", "-n", help="Maximum results"),
    db: str = typer.Option(DEFAULT_DB, "--db", "-d", help="Path to DuckDB database"),
):
    """
    List proteins with optional filtering.

    Examples:
        sharur proteins --genome bin_001
        sharur proteins --min-length 2000 --unannotated
    """
    db_path = Path(db)
    if not db_path.exists():
        typer.echo(f"DB not found: {db}", err=True)
        raise typer.Exit(1)

    b = Sharur(db_path)
    result = b.list_proteins(
        genome_id=genome,
        contig_id=contig,
        min_length=min_length,
        max_length=max_length,
        has_annotation=annotated,
        limit=limit,
    )
    typer.echo(result.data)

    if result.meta.truncated:
        typer.echo(f"\n[Showing {result.meta.rows} of {result.meta.total_rows} results]", err=True)


# ------------------------------------------------------------------ #
# Neighborhood command
# ------------------------------------------------------------------ #


@app.command()
def neighborhood(
    protein_id: str = typer.Argument(..., help="Protein ID as anchor"),
    window: int = typer.Option(10, "--window", "-w", help="Genes on each side"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Show predicate details"),
    db: str = typer.Option(DEFAULT_DB, "--db", "-d", help="Path to DuckDB database"),
):
    """
    Show genomic neighborhood around a protein.

    Displays proteins in context with coordinates, annotations,
    and predicates as an ASCII table.

    Example:
        sharur neighborhood PROTEIN_ID --window 15
    """
    db_path = Path(db)
    if not db_path.exists():
        typer.echo(f"DB not found: {db}", err=True)
        raise typer.Exit(1)

    b = Sharur(db_path)
    verbosity = 2 if verbose else 1
    result = b.get_neighborhood(entity_id=protein_id, window=window, verbosity=verbosity)
    typer.echo(result.data)


# ------------------------------------------------------------------ #
# Search command
# ------------------------------------------------------------------ #


@app.command()
def search(
    has: Optional[str] = typer.Option(None, "--has", help="Predicates that must be true (comma-separated)"),
    lacks: Optional[str] = typer.Option(None, "--lacks", help="Predicates that must be false (comma-separated)"),
    annotation: Optional[str] = typer.Option(None, "--annotation", "-a", help="Annotation pattern to match"),
    accession: Optional[str] = typer.Option(None, "--accession", help="Exact accession (e.g., PF00142)"),
    taxonomy: Optional[str] = typer.Option(None, "--taxonomy", "-t", help="Taxonomy filter"),
    limit: int = typer.Option(50, "--limit", "-n", help="Maximum results"),
    db: str = typer.Option(DEFAULT_DB, "--db", "-d", help="Path to DuckDB database"),
):
    """
    Search proteins by predicates or annotations.

    Predicate search uses set logic:
    - --has: Protein must have ALL specified predicates (AND)
    - --lacks: Protein must have NONE of specified predicates

    Examples:
        sharur search --has giant,unannotated
        sharur search --has confident_hit --lacks hypothetical
        sharur search --annotation "hydrogenase"
        sharur search --accession PF00142
    """
    db_path = Path(db)
    if not db_path.exists():
        typer.echo(f"DB not found: {db}", err=True)
        raise typer.Exit(1)

    b = Sharur(db_path)

    # Parse comma-separated predicates
    has_list = [p.strip() for p in has.split(",")] if has else None
    lacks_list = [p.strip() for p in lacks.split(",")] if lacks else None

    # Decide which search to use
    if has_list or lacks_list:
        result = b.search_by_predicates(has=has_list, lacks=lacks_list, limit=limit)
    elif annotation or accession or taxonomy:
        result = b.search_proteins(
            annotation_pattern=annotation,
            accession=accession,
            taxonomy_filter=taxonomy,
            limit=limit,
        )
    else:
        typer.echo("Please specify search criteria: --has, --lacks, --annotation, --accession, or --taxonomy", err=True)
        raise typer.Exit(1)

    typer.echo(result.data)

    if result.meta.truncated:
        typer.echo(f"\n[Showing {result.meta.rows} of {result.meta.total_rows} results]", err=True)


# ------------------------------------------------------------------ #
# Compute predicates command
# ------------------------------------------------------------------ #


@app.command(name="compute-predicates")
def compute_predicates(
    db: str = typer.Option(DEFAULT_DB, "--db", "-d", help="Path to DuckDB database"),
    protein_id: Optional[str] = typer.Option(None, "--protein", "-p", help="Compute for specific protein only"),
):
    """
    Compute and store predicates for proteins.

    This analyzes annotations and properties to compute semantic
    predicates like 'transporter', 'hydrogenase', 'giant', etc.

    Run this after loading annotations to enable predicate-based search.

    Examples:
        sharur compute-predicates --db data/sharur.duckdb
        sharur compute-predicates --protein PROT_001 --db data/sharur.duckdb
    """
    from sharur.predicates.generator import (
        generate_predicates_for_proteins,
        persist_generated_predicates,
    )
    from sharur.storage.duckdb_store import DuckDBStore

    db_path = Path(db)
    if not db_path.exists():
        typer.echo(f"DB not found: {db}", err=True)
        raise typer.Exit(1)

    store = DuckDBStore(db_path=db_path)

    if protein_id:
        typer.echo(f"Computing predicates for {protein_id}...", err=True)
        predicates = generate_predicates_for_proteins(store, protein_ids=[protein_id])
        if protein_id in predicates:
            preds = predicates[protein_id]
            typer.echo(f"\nPredicates for {protein_id}:")
            for pred in preds:
                typer.echo(f"  - {pred}")
            typer.echo(f"\nTotal: {len(preds)} predicates")
        else:
            typer.echo(f"No predicates generated for {protein_id}")
    else:
        typer.echo("Computing predicates for all proteins...", err=True)
        predicates = generate_predicates_for_proteins(store)
        count = persist_generated_predicates(store, predicates)
        typer.echo(f"\nComputed and stored predicates for {count} proteins")

        # Show summary
        all_preds = set()
        for preds in predicates.values():
            all_preds.update(preds)
        typer.echo(f"Unique predicates: {len(all_preds)}")


# ------------------------------------------------------------------ #
# Predicates command (list available predicates)
# ------------------------------------------------------------------ #


@app.command()
def predicates(
    category: Optional[str] = typer.Option(None, "--category", "-c", help="Filter by category"),
    show_hierarchy: bool = typer.Option(False, "--hierarchy", "-h", help="Show parent predicates"),
):
    """
    List available predicates for search.

    Predicates are boolean properties computed over proteins.
    Use with 'sharur search --has PREDICATE'.

    Categories include:
    - enzyme: Enzyme classes (oxidoreductase, hydrolase, etc.)
    - transport: Transporters and membrane proteins
    - regulation: Regulators and signaling
    - metabolism: Metabolic pathway enzymes
    - cazy: Carbohydrate-active enzymes
    - binding: Binding domains
    - envelope: Cell surface and envelope
    - mobile: Mobile elements and defense
    - stress: Stress response and resistance
    - size: Size-based (tiny, giant, massive)
    - annotation: Annotation status

    Examples:
        sharur predicates
        sharur predicates --category enzyme
        sharur predicates --category transport --hierarchy
    """
    from sharur.predicates.vocabulary import list_predicates, list_categories

    preds = list_predicates(category=category)

    if not preds:
        if category:
            typer.echo(f"No predicates found in category '{category}'.\n")
            typer.echo("Available categories:")
            for cat in list_categories():
                typer.echo(f"  - {cat}")
        return

    typer.echo("# Available Predicates\n")

    current_category = None
    for pred in preds:
        if pred.category != current_category:
            current_category = pred.category
            typer.echo(f"\n## {current_category.title()}\n")

        if show_hierarchy and pred.parent:
            typer.echo(f"**{pred.predicate_id}** ({pred.parent}): {pred.description}")
        else:
            typer.echo(f"**{pred.predicate_id}**: {pred.description}")

    typer.echo(f"\nTotal: {len(preds)} predicates")


# ------------------------------------------------------------------ #
# Main entry point
# ------------------------------------------------------------------ #


def main():
    app(prog_name="sharur")


if __name__ == "__main__":
    main()
