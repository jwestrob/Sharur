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

import json
import os
from enum import Enum
from pathlib import Path
from typing import Optional

import typer

from sharur import __version__
from sharur.operators import Sharur

# Default DB path
DEFAULT_DB = "data/sharur.duckdb"


class OutputFormat(str, Enum):
    markdown = "markdown"
    json = "json"
    jsonl = "jsonl"
    tsv = "tsv"


class BriefFormat(str, Enum):
    markdown = "markdown"
    json = "json"

# Main app with subcommands
app = typer.Typer(
    no_args_is_help=True,
    add_completion=False,
    help="Sharur - Metagenomic dataset exploration CLI",
    rich_markup_mode=None,  # Disable rich to work around Typer 0.15 bug
)


def _version_callback(value: bool) -> None:
    if value:
        typer.echo(f"sharur {__version__}")
        raise typer.Exit()


@app.callback()
def _main(
    version: bool = typer.Option(
        False,
        "--version",
        callback=_version_callback,
        is_eager=True,
        help="Show the Sharur version and exit.",
    ),
) -> None:
    """Sharur - Metagenomic dataset exploration CLI."""


def _emit_result(result, output_format: OutputFormat) -> None:
    """Render one operator result in a human- or machine-readable format."""
    if output_format == OutputFormat.markdown:
        typer.echo(result.data)
        return
    if output_format == OutputFormat.json:
        typer.echo(result.to_json(indent=2))
        return
    if output_format == OutputFormat.jsonl:
        for record in result.records:
            typer.echo(json.dumps(record, default=str))
        return

    typer.echo(result.to_dataframe().to_csv(sep="\t", index=False), nl=False)


# ------------------------------------------------------------------ #
# Overview command
# ------------------------------------------------------------------ #


@app.command()
def overview(
    db: str = typer.Option(DEFAULT_DB, "--db", "-d", help="Path to DuckDB database"),
    output_format: OutputFormat = typer.Option(
        OutputFormat.markdown,
        "--format",
        "-f",
        help="Output format: markdown, json, jsonl, or tsv",
    ),
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

    b = Sharur(db_path, read_only=True)
    result = b.overview()
    _emit_result(result, output_format)


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
    output_format: OutputFormat = typer.Option(
        OutputFormat.markdown,
        "--format",
        "-f",
    ),
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

    b = Sharur(db_path, read_only=True)
    result = b.list_genomes(
        taxonomy_filter=taxonomy,
        min_completeness=min_completeness,
        max_contamination=max_contamination,
        limit=limit,
    )
    _emit_result(result, output_format)

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
    output_format: OutputFormat = typer.Option(
        OutputFormat.markdown,
        "--format",
        "-f",
    ),
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

    b = Sharur(db_path, read_only=True)
    result = b.list_proteins(
        genome_id=genome,
        contig_id=contig,
        min_length=min_length,
        max_length=max_length,
        has_annotation=annotated,
        limit=limit,
    )
    _emit_result(result, output_format)

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
    output_format: OutputFormat = typer.Option(
        OutputFormat.markdown,
        "--format",
        "-f",
    ),
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

    b = Sharur(db_path, read_only=True)
    verbosity = 2 if verbose else 1
    result = b.get_neighborhood(entity_id=protein_id, window=window, verbosity=verbosity)
    _emit_result(result, output_format)


# ------------------------------------------------------------------ #
# First-class case inspection and context comparison
# ------------------------------------------------------------------ #


@app.command(name="inspect")
def inspect_entity(
    entity_id: str = typer.Argument(..., help="Protein, system, locus, contig, or bin ID"),
    entity_type: str | None = typer.Option(
        None,
        "--type",
        help="Disambiguate as protein, system, locus, contig, or bin.",
    ),
    bin_id: str | None = typer.Option(
        None,
        "--bin",
        help="Required only when a contig label occurs in multiple bins.",
    ),
    source_table: str | None = typer.Option(
        None,
        "--source-table",
        help="Structured caller table when a system/locus ID is ambiguous.",
    ),
    window: int = typer.Option(
        10,
        "--window",
        "-w",
        help="Default ORFs on each side.",
    ),
    upstream_orfs: int | None = typer.Option(
        None,
        "--upstream",
        help="Biological upstream ORFs; overrides --window on that side.",
    ),
    downstream_orfs: int | None = typer.Option(
        None,
        "--downstream",
        help="Biological downstream ORFs; overrides --window on that side.",
    ),
    assembly_evidence: Path | None = typer.Option(
        None,
        "--assembly-evidence",
        help="Optional assembly_evidence.duckdb sidecar.",
    ),
    include_sequences: bool = typer.Option(
        False,
        "--include-sequences",
        help="Embed context sequences in JSON output (can make output large).",
    ),
    plot: Path | None = typer.Option(
        None,
        "--plot",
        help="Render the resolved locus to this PNG/SVG path.",
    ),
    bundle: Path | None = typer.Option(
        None,
        "--bundle",
        help="Write a compact, replayable evidence-bundle directory.",
    ),
    bundle_sequences: bool = typer.Option(
        True,
        "--bundle-sequences/--no-bundle-sequences",
        help="Include anchor-component FASTA in --bundle.",
    ),
    overwrite: bool = typer.Option(
        False,
        "--overwrite",
        help="Replace an existing --bundle directory.",
    ),
    db: str = typer.Option(DEFAULT_DB, "--db", "-d", help="Path to DuckDB database"),
    output_format: BriefFormat = typer.Option(
        BriefFormat.markdown,
        "--format",
        "-f",
        help="Output format: markdown or json.",
    ),
):
    """Resolve an entity into a provenance-separated, strand-aware case."""
    db_path = Path(db)
    if not db_path.is_file():
        typer.echo(f"DB not found: {db}", err=True)
        raise typer.Exit(1)

    try:
        sharur = Sharur(
            db_path,
            read_only=True,
            assembly_evidence_path=assembly_evidence,
        )
        case = sharur.inspect(
            entity_id,
            entity_type=entity_type,
            bin_id=bin_id,
            source_table=source_table,
            window=window,
            upstream_orfs=upstream_orfs,
            downstream_orfs=downstream_orfs,
            include_sequences=include_sequences,
        )
        if plot is not None:
            rendered = case.plot(plot)
            typer.echo(f"Plot: {rendered}", err=True)
        if bundle is not None:
            bundle_path = case.export(
                bundle,
                include_sequences=bundle_sequences,
                include_plot=False,
                overwrite=overwrite,
            )
            typer.echo(f"Bundle: {bundle_path}", err=True)
    except (KeyError, ValueError, RuntimeError, FileExistsError) as exc:
        typer.echo(f"Could not inspect case: {exc}", err=True)
        raise typer.Exit(1) from exc

    if output_format == BriefFormat.json:
        typer.echo(case.to_json())
    else:
        typer.echo(case.to_markdown())


@app.command(name="compare-context")
def compare_context_command(
    entity_id: str = typer.Argument(..., help="System, locus, or protein case ID"),
    feature: list[str] = typer.Option(
        [],
        "--feature",
        help=(
            "Repeatable feature: pfam:ACCESSION, name:TEXT, predicate:ID, "
            "system:TYPE, locus:TYPE, or other_called_system."
        ),
    ),
    entity_type: str | None = typer.Option(None, "--type"),
    source_table: str | None = typer.Option(None, "--source-table"),
    bin_id: str | None = typer.Option(None, "--bin"),
    window: int = typer.Option(10, "--window", "-w", help="Default ORFs on each side."),
    upstream_orfs: int | None = typer.Option(
        None,
        "--upstream",
        help="Biological upstream ORFs; overrides --window.",
    ),
    downstream_orfs: int | None = typer.Option(
        None,
        "--downstream",
        help="Biological downstream ORFs; overrides --window.",
    ),
    foreground_id: list[str] = typer.Option(
        [],
        "--foreground-id",
        help="Explicit foreground entity ID; repeat for protein/custom cohorts.",
    ),
    background_id: list[str] = typer.Option(
        [],
        "--background-id",
        help="Explicit background entity ID; repeat for protein/custom cohorts.",
    ),
    combine: str = typer.Option(
        "all",
        "--combine",
        help="Combine features with all or any.",
    ),
    min_components: int = typer.Option(
        1,
        "--min-components",
        help="Exclude caller-emitted systems/loci with fewer components.",
    ),
    require_full_context: bool = typer.Option(
        True,
        "--require-full-context/--allow-edge-censored",
        help="Exclude contig-edge-censored neighborhoods.",
    ),
    deduplicate_by: str = typer.Option(
        "replicon",
        "--deduplicate-by",
        help="Independent unit: entity, replicon, or bin.",
    ),
    exclude_foreground_units: bool = typer.Option(
        True,
        "--exclude-foreground-units/--allow-foreground-overlap",
        help="Keep foreground-bearing units out of the background.",
    ),
    taxonomy: str | None = typer.Option(
        None,
        "--taxonomy",
        help="Require this taxonomy substring in both groups.",
    ),
    same_taxonomy_rank: str | None = typer.Option(
        None,
        "--same-taxonomy-rank",
        help="Match the case at domain/phylum/class/order/family/genus/species.",
    ),
    alternative: str = typer.Option(
        "greater",
        "--alternative",
        help="Fisher alternative: greater, less, or two-sided.",
    ),
    bundle: Path | None = typer.Option(
        None,
        "--bundle",
        help="Write the case, comparison matrix, recipe, and verifier.",
    ),
    overwrite: bool = typer.Option(False, "--overwrite"),
    db: str = typer.Option(DEFAULT_DB, "--db", "-d", help="Path to DuckDB database"),
    output_format: BriefFormat = typer.Option(
        BriefFormat.markdown,
        "--format",
        "-f",
    ),
):
    """Run an exact, reproducible foreground/background ORF-context comparison."""
    if not feature:
        typer.echo("At least one --feature is required.", err=True)
        raise typer.Exit(1)
    if combine not in {"all", "any"}:
        raise typer.BadParameter("--combine must be all or any")
    if deduplicate_by not in {"entity", "replicon", "bin"}:
        raise typer.BadParameter("--deduplicate-by must be entity, replicon, or bin")
    if alternative not in {"greater", "less", "two-sided"}:
        raise typer.BadParameter("--alternative must be greater, less, or two-sided")

    db_path = Path(db)
    if not db_path.is_file():
        typer.echo(f"DB not found: {db}", err=True)
        raise typer.Exit(1)
    try:
        sharur = Sharur(db_path, read_only=True)
        case = sharur.inspect(
            entity_id,
            entity_type=entity_type,
            bin_id=bin_id,
            source_table=source_table,
            window=window,
            upstream_orfs=upstream_orfs,
            downstream_orfs=downstream_orfs,
        )
        comparison = case.compare_context(
            features=feature,
            window=window,
            upstream_orfs=upstream_orfs,
            downstream_orfs=downstream_orfs,
            foreground_ids=foreground_id or None,
            background_ids=background_id or None,
            combine=combine,
            min_components=min_components,
            require_full_context=require_full_context,
            deduplicate_by=deduplicate_by,
            exclude_foreground_units=exclude_foreground_units,
            taxonomy_filter=taxonomy,
            same_taxonomy_rank=same_taxonomy_rank,
            alternative=alternative,
        )
        if bundle is not None:
            bundle_path = case.export(
                bundle,
                comparison=comparison,
                overwrite=overwrite,
            )
            typer.echo(f"Bundle: {bundle_path}", err=True)
    except (KeyError, ValueError, RuntimeError, FileExistsError) as exc:
        typer.echo(f"Could not compare context: {exc}", err=True)
        raise typer.Exit(1) from exc

    if output_format == BriefFormat.json:
        typer.echo(comparison.to_json())
    else:
        typer.echo(comparison.to_markdown())


# ------------------------------------------------------------------ #
# Optional assembly/host-assignment evidence
# ------------------------------------------------------------------ #


@app.command(name="import-assembly-evidence")
def import_assembly_evidence_command(
    input_path: Path = typer.Argument(
        ...,
        exists=True,
        dir_okay=False,
        readable=True,
        help="TSV, CSV, or JSONL with bin_id and contig_id columns.",
    ),
    db: Path = typer.Option(
        Path(DEFAULT_DB),
        "--db",
        "-d",
        help="Core DuckDB used for validation and sidecar discovery.",
    ),
    sidecar: Path | None = typer.Option(
        None,
        "--sidecar",
        help="Output sidecar (default: assembly_evidence.duckdb beside --db).",
    ),
    source: str | None = typer.Option(
        None,
        "--source",
        help="Provenance label for these measurements.",
    ),
    validate: bool = typer.Option(
        True,
        "--validate/--no-validate",
        help="Require every (bin_id, contig_id) to exist in the core dataset.",
    ),
    hash_input: bool = typer.Option(
        True,
        "--hash-input/--skip-input-hash",
        help="SHA-256 the input for provenance (one extra sequential read).",
    ),
):
    """Import optional scalar contig evidence without modifying the core database."""
    from sharur.assembly_evidence import (  # noqa: PLC0415
        default_assembly_evidence_path,
        import_contig_evidence,
    )

    if not db.is_file():
        typer.echo(f"DB not found: {db}", err=True)
        raise typer.Exit(1)
    target = sidecar or default_assembly_evidence_path(db)
    if target.expanduser().resolve() == db.expanduser().resolve():
        raise typer.BadParameter(
            "--sidecar must not be the core DuckDB",
            param_hint="--sidecar",
        )
    try:
        result = import_contig_evidence(
            input_path,
            target,
            source=source,
            validate_db_path=db if validate else None,
            hash_input=hash_input,
        )
    except (OSError, ValueError) as exc:
        typer.echo(f"Could not import assembly evidence: {exc}", err=True)
        raise typer.Exit(1) from exc
    typer.echo(json.dumps(result, indent=2, default=str))


@app.command(name="compute-composition-evidence")
def compute_composition_evidence_command(
    assembly: list[str] = typer.Option(
        [],
        "--assembly",
        help="Repeatable BIN_ID=/path/to/assembly.fna mapping.",
    ),
    db: Path = typer.Option(
        Path(DEFAULT_DB),
        "--db",
        "-d",
        help="Core DuckDB used for validation and sidecar discovery.",
    ),
    sidecar: Path | None = typer.Option(
        None,
        "--sidecar",
        help="Output sidecar (default: assembly_evidence.duckdb beside --db).",
    ),
    validate: bool = typer.Option(
        True,
        "--validate/--no-validate",
        help="Require FASTA contigs to exist in the core dataset.",
    ),
):
    """Explicitly scan FASTAs for scalar GC/4-mer evidence.

    This command is never run by ingestion, preflight, or case inspection.
    It reads every supplied FASTA, keeps 4-mer vectors in memory only, and
    persists scalar leave-one-contig-out distances.
    """
    from sharur.assembly_evidence import (  # noqa: PLC0415
        compute_composition_evidence,
        default_assembly_evidence_path,
    )

    if not assembly:
        typer.echo("At least one --assembly BIN_ID=FASTA is required.", err=True)
        raise typer.Exit(1)
    if not db.is_file():
        typer.echo(f"DB not found: {db}", err=True)
        raise typer.Exit(1)

    assemblies: dict[str, Path] = {}
    for specification in assembly:
        bin_name, separator, fasta = specification.partition("=")
        if not separator or not bin_name.strip() or not fasta.strip():
            raise typer.BadParameter(
                f"Invalid --assembly {specification!r}; use BIN_ID=/path/to/file.fna"
            )
        if bin_name in assemblies:
            raise typer.BadParameter(f"Duplicate --assembly bin ID: {bin_name}")
        assemblies[bin_name] = Path(fasta)

    target = sidecar or default_assembly_evidence_path(db)
    if target.expanduser().resolve() == db.expanduser().resolve():
        raise typer.BadParameter(
            "--sidecar must not be the core DuckDB",
            param_hint="--sidecar",
        )
    typer.echo(
        f"Explicit composition scan: {len(assemblies)} FASTA file(s); "
        "4-mer vectors will not be persisted.",
        err=True,
    )
    try:
        result = compute_composition_evidence(
            assemblies,
            target,
            validate_db_path=db if validate else None,
        )
    except (OSError, ValueError) as exc:
        typer.echo(f"Could not compute composition evidence: {exc}", err=True)
        raise typer.Exit(1) from exc
    typer.echo(json.dumps(result, indent=2, default=str))


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
    output_format: OutputFormat = typer.Option(
        OutputFormat.markdown,
        "--format",
        "-f",
    ),
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

    b = Sharur(db_path, read_only=True)

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

    _emit_result(result, output_format)

    if result.meta.truncated:
        typer.echo(f"\n[Showing {result.meta.rows} of {result.meta.total_rows} results]", err=True)


# ------------------------------------------------------------------ #
# Compute predicates command
# ------------------------------------------------------------------ #


@app.command(name="compute-predicates")
def compute_predicates(
    db: str = typer.Option(DEFAULT_DB, "--db", "-d", help="Path to DuckDB database"),
    protein_id: Optional[str] = typer.Option(None, "--protein", "-p", help="Compute for specific protein only"),
    chunk_size: int = typer.Option(100_000, "--chunk-size", help="V2 generation batch size"),
    workers: Optional[int] = typer.Option(
        None,
        "--workers",
        "-t",
        help="Transform processes (defaults to Slurm CPUs, then host CPUs)",
    ),
    worker_batch_size: Optional[int] = typer.Option(
        None,
        "--worker-batch-size",
        help="Optional proteins per transform task",
    ),
    pipeline_depth: int = typer.Option(
        2,
        "--pipeline-depth",
        min=1,
        help="Bounded full-dataset chunks overlapped across read/transform/write",
    ),
    resume: bool = typer.Option(False, "--resume", help="Continue the latest full V2 checkpoint"),
    review_queue: Optional[Path] = typer.Option(
        None,
        "--review-queue",
        help="Write unresolved-accession review TSV",
    ),
):
    """
    Compute and store V2 predicates for proteins.

    This writes semantic_atoms and semantic_state, then materializes the
    V2-derived protein_predicates compatibility table used by legacy
    search and reports.

    Run this after loading annotations to enable predicate-based search.

    Examples:
        sharur compute-predicates --db data/sharur.duckdb
        sharur compute-predicates --protein PROT_001 --db data/sharur.duckdb
    """
    from sharur.predicates_v2.persistence import generate_and_persist_v2
    from sharur.storage.duckdb_store import DuckDBStore

    db_path = Path(db)
    if not db_path.exists():
        typer.echo(f"DB not found: {db}", err=True)
        raise typer.Exit(1)

    store = DuckDBStore(db_path=db_path)
    if workers is None:
        slurm_workers = os.environ.get("SLURM_CPUS_ON_NODE")
        try:
            workers = int(slurm_workers) if slurm_workers else (os.cpu_count() or 1)
        except ValueError as exc:
            raise typer.BadParameter(
                "SLURM_CPUS_ON_NODE must be a positive integer"
            ) from exc
    if workers <= 0:
        raise typer.BadParameter("--workers must be a positive integer")

    if protein_id:
        if resume:
            typer.echo("--resume applies to full-dataset generation", err=True)
            raise typer.Exit(2)
        typer.echo(f"Computing V2 predicates for {protein_id}...", err=True)
        states = generate_and_persist_v2(
            store,
            protein_ids=[protein_id],
            output_review_queue=str(review_queue) if review_queue else None,
            chunk_size=chunk_size,
            workers=workers,
            worker_batch_size=worker_batch_size,
            pipeline_depth=pipeline_depth,
            update_legacy_predicates=True,
            return_states=True,
            predict_topology=False,
        )
        if protein_id not in states:
            typer.echo(f"No predicates generated for {protein_id}")
            return

        rows = store.execute(
            "SELECT predicates FROM protein_predicates WHERE protein_id = ?",
            [protein_id],
        )
        preds = rows[0][0] if rows else []
        typer.echo(f"\nPredicates for {protein_id}:")
        for pred in preds:
            typer.echo(f"  - {pred}")
        typer.echo(f"\nTotal: {len(preds)} predicates")
    else:
        typer.echo("Computing V2 predicates for all proteins...", err=True)
        generate_and_persist_v2(
            store,
            output_review_queue=str(review_queue) if review_queue else None,
            chunk_size=chunk_size,
            workers=workers,
            worker_batch_size=worker_batch_size,
            pipeline_depth=pipeline_depth,
            resume=resume,
            update_legacy_predicates=True,
            return_states=False,
            predict_topology=False,
        )
        count = store.execute("SELECT COUNT(*) FROM protein_predicates")[0][0]
        atoms = store.execute("SELECT COUNT(*) FROM semantic_atoms")[0][0]
        states = store.execute("SELECT COUNT(*) FROM semantic_state")[0][0]
        terms = store.execute("SELECT COUNT(*) FROM semantic_terms")[0][0]
        system_members = store.execute("SELECT COUNT(*) FROM system_proteins")[0][0]
        unique = store.execute(
            "SELECT COUNT(DISTINCT atom_id) FROM semantic_atoms"
        )[0][0]
        typer.echo(
            f"\nComputed V2 predicates for {count} proteins "
            f"({states} semantic states, {atoms} atoms, {terms} search terms, "
            f"{unique} unique atoms, {system_members} system members)"
        )


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
# Capability preflight and vector-index preparation
# ------------------------------------------------------------------ #


@app.command()
def preflight(
    db: str = typer.Option(DEFAULT_DB, "--db", "-d", help="Path to DuckDB database"),
    assembly_evidence: Path | None = typer.Option(
        None,
        "--assembly-evidence",
        help="Optional non-default assembly-evidence sidecar.",
    ),
    synteny: Path | None = typer.Option(
        None,
        "--synteny",
        help="Optional non-default normalized ELSA sidecar.",
    ),
    output_format: BriefFormat = typer.Option(
        BriefFormat.markdown,
        "--format",
        "-f",
        help="Output format: markdown or json",
    ),
    skip_tools: bool = typer.Option(
        False,
        "--skip-tools",
        help="Skip slower external binary/version probes.",
    ),
    strict: bool = typer.Option(
        False,
        "--strict",
        help="Exit non-zero unless every required dataset capability is available.",
    ),
):
    """Emit one typed dataset/runtime capability brief without modifying data."""
    from sharur.capabilities import CapabilityState, build_capability_brief

    brief = build_capability_brief(
        db,
        include_tools=not skip_tools,
        assembly_evidence_path=assembly_evidence,
        synteny_path=synteny,
    )
    if output_format == BriefFormat.json:
        typer.echo(brief.to_json())
    else:
        typer.echo(brief.to_markdown())
    if strict and brief.overall_state != CapabilityState.available:
        raise typer.Exit(code=1)


@app.command()
def seal(
    db: Path = typer.Option(
        Path(DEFAULT_DB),
        "--db",
        "-d",
        help="Path to the dataset DuckDB.",
    ),
    output: Path | None = typer.Option(
        None,
        "--output",
        "-o",
        help="Seal path (default: DATASET/dataset.seal.json).",
    ),
    full: bool = typer.Option(
        False,
        "--full",
        help="Fully hash every canonical artifact, including large DuckDB/H5/index files.",
    ),
    include_tools: bool = typer.Option(
        False,
        "--include-tools",
        help="Record slower external tool and reference-database probes as provenance.",
    ),
    force: bool = typer.Option(
        False,
        "--force",
        help="Replace an existing seal atomically.",
    ),
):
    """Write a portable integrity seal for a completed dataset."""
    from sharur.dataset_seal import DatasetSealError, seal_dataset  # noqa: PLC0415

    try:
        seal_path, payload = seal_dataset(
            db,
            output_path=output,
            hash_large_files=full,
            include_tools=include_tools,
            force=force,
        )
    except (DatasetSealError, FileExistsError, ValueError) as exc:
        typer.echo(f"Could not seal dataset: {exc}", err=True)
        raise typer.Exit(1) from exc

    typer.echo(f"Dataset ID: {payload['dataset_id']}")
    typer.echo(f"Strength: {payload['seal_strength']}")
    typer.echo(f"Seal: {seal_path}")


@app.command()
def migrate(
    db: Path = typer.Option(
        Path(DEFAULT_DB),
        "--db",
        "-d",
        help="Path to a writable Sharur DuckDB.",
    ),
):
    """Apply pending additive schema/index migrations in a maintenance window."""
    import duckdb  # noqa: PLC0415

    from sharur.storage.migrations import (  # noqa: PLC0415
        get_current_version,
        run_migrations,
    )
    from sharur.storage.schema import SCHEMA_VERSION  # noqa: PLC0415

    database = db.expanduser().resolve()
    if not database.is_file():
        typer.echo(f"DuckDB file does not exist: {database}", err=True)
        raise typer.Exit(1)
    connection = duckdb.connect(str(database))
    try:
        before = get_current_version(connection)
        applied = run_migrations(connection)
        after = get_current_version(connection)
    finally:
        connection.close()
    typer.echo(
        f"Schema: {before} -> {after} "
        f"({applied} migration{'s' if applied != 1 else ''} applied)"
    )
    if after != SCHEMA_VERSION:
        typer.echo(
            f"Current software expects schema {SCHEMA_VERSION}; migration stopped at {after}.",
            err=True,
        )
        raise typer.Exit(1)
    if applied:
        typer.echo(
            "Canonical database state changed. Rebuild dataset.seal.json "
            "with `sharur seal --force` before serving it."
        )


@app.command(name="backfill-contig-context")
def backfill_contig_context_command(
    db: Path = typer.Option(Path(DEFAULT_DB), "--db", "-d", help="Path to a writable Sharur DuckDB."),
    assemblies: Path | None = typer.Option(
        None, "--assemblies",
        help="Directory of assembly FASTAs named BIN_ID.fna|fa|fasta[.gz] "
             "(default: the dataset's stage00_prepared manifest)."),
    proteins: Path | None = typer.Option(
        None, "--proteins",
        help="Prodigal output directory with *.faa files (default: the dataset's stage03_prodigal)."),
):
    """Record assembly contig lengths and Prodigal truncation flags in an existing database.

    Older databases store each contig's length as the end of its last gene and
    omit Prodigal's partial-gene flags; contig-edge context needs both.
    Writes the canonical database: run in a maintenance window, then reseal.
    """
    import duckdb  # noqa: PLC0415

    from sharur.contig_context import (  # noqa: PLC0415
        assemblies_from_dir,
        assemblies_from_stage00,
        backfill_contig_context,
    )

    database = db.expanduser().resolve()
    if not database.is_file():
        typer.echo(f"DuckDB file does not exist: {database}", err=True)
        raise typer.Exit(1)
    root = database.parent
    assembly_paths = (assemblies_from_dir(assemblies) if assemblies
                      else assemblies_from_stage00(root / "stage00_prepared"))
    protein_dir = proteins or root / "stage03_prodigal"
    found = [f for f in protein_dir.glob("**/*.faa") if f.parent.name != "all_protein_symlinks"] \
        if protein_dir.is_dir() else []
    faas = [f for f in found if f.is_file()]
    if len(faas) < len(found):
        typer.echo(f"Skipping {len(found) - len(faas):,} unreadable protein files (broken links?)", err=True)
    if not assembly_paths and not faas:
        typer.echo("No assemblies or Prodigal protein files found; pass --assemblies/--proteins.", err=True)
        raise typer.Exit(1)
    connection = duckdb.connect(str(database))
    try:
        stats = backfill_contig_context(connection, assemblies=assembly_paths, protein_faas=faas)
    except ValueError as exc:
        typer.echo(str(exc), err=True)
        raise typer.Exit(1) from exc
    finally:
        connection.close()
    typer.echo(f"Assemblies read: {len(assembly_paths):,}; protein files read: {len(faas):,}")
    typer.echo(f"Contig lengths from assemblies: {stats['contigs_updated']:,} "
               f"({stats['contigs_unmatched']:,} contigs in those bins without an assembly record)")
    typer.echo(f"Proteins with Prodigal partial flags: {stats['proteins_flagged']:,}")
    typer.echo("Canonical database state changed. Rebuild dataset.seal.json with `sharur seal --force`.")


@app.command(name="verify-seal")
def verify_seal(
    seal_path: Path = typer.Argument(..., help="Path to dataset.seal.json."),
    db: Path | None = typer.Option(
        None,
        "--db",
        "-d",
        help="Current DuckDB path (needed if the seal was moved separately).",
    ),
    output_format: BriefFormat = typer.Option(
        BriefFormat.markdown,
        "--format",
        "-f",
        help="Output format: markdown or json.",
    ),
):
    """Recompute a dataset seal and report canonical identity drift."""
    import json  # noqa: PLC0415

    from sharur.dataset_seal import (  # noqa: PLC0415
        DatasetSealError,
        verify_dataset_seal,
    )

    try:
        verification = verify_dataset_seal(seal_path, db_path=db)
    except (DatasetSealError, OSError, ValueError) as exc:
        typer.echo(f"Could not verify dataset seal: {exc}", err=True)
        raise typer.Exit(1) from exc

    if output_format == BriefFormat.json:
        typer.echo(json.dumps(verification.to_dict(), indent=2))
    else:
        typer.echo(verification.to_markdown())
    if not verification.valid:
        raise typer.Exit(1)


@app.command(name="build-vector-index")
def build_vector_index_command(
    embeddings: Optional[Path] = typer.Option(
        None,
        "--embeddings",
        "-e",
        help="Canonical protein_embeddings.h5 path.",
    ),
    db: Optional[Path] = typer.Option(
        None,
        "--db",
        "-d",
        help="DuckDB path; discovers embeddings beside the dataset.",
    ),
    force: bool = typer.Option(False, "--force", help="Rebuild an already-valid index."),
    chunk_size: int = typer.Option(
        50_000,
        "--chunk-size",
        help="Number of vectors streamed per FAISS add batch.",
    ),
    nprobe: int = typer.Option(32, "--nprobe", help="IVF partitions probed per query."),
    threads: Optional[int] = typer.Option(
        None,
        "--threads",
        "-t",
        help="FAISS CPU build threads (default: FAISS runtime default).",
    ),
):
    """Build mmap-ready FAISS sidecars and a stable disk-backed protein-ID map."""
    import json

    from sharur.storage.vector_store import build_vector_index

    h5_path = embeddings
    if h5_path is None and db is not None:
        for directory in ("embeddings", "stage06_embeddings"):
            candidate = db.expanduser().resolve().parent / directory / "protein_embeddings.h5"
            if candidate.is_file():
                h5_path = candidate
                break
    if h5_path is None:
        raise typer.BadParameter("Provide --embeddings or a --db with discovered embeddings.")
    if not h5_path.expanduser().is_file():
        typer.echo(f"Embeddings not found: {h5_path}", err=True)
        raise typer.Exit(1)

    inspection = build_vector_index(
        h5_path,
        force=force,
        chunk_size=chunk_size,
        nprobe=nprobe,
        threads=threads,
    )
    typer.echo(json.dumps(inspection.to_dict(), indent=2))


# ------------------------------------------------------------------ #
# Doctor command — install verification
# ------------------------------------------------------------------ #


@app.command()
def doctor(
    strict: bool = typer.Option(
        False,
        "--strict",
        help="Exit non-zero if any core tool/database is missing.",
    ),
):
    """Verify external tools, reference databases, and API keys are available.

    Default is informational (always exits 0) so it is safe as a CI smoke check
    and container HEALTHCHECK. Pass --strict to fail when a core component is
    missing.
    """
    from rich.console import Console
    from rich.table import Table

    from sharur import diagnostics

    checks = diagnostics.run_all_checks()

    console = Console()
    table = Table(title=f"sharur doctor  (v{__version__})", show_lines=False)
    table.add_column("", width=6)
    table.add_column("Component", style="bold")
    table.add_column("Version / Detail")
    table.add_column("Purpose", style="dim")

    style_for = {
        diagnostics.OK: ("[green]ok[/green]", ""),
        diagnostics.WARN: ("[yellow]WARN[/yellow]", ""),
        diagnostics.MISSING: ("[red]MISS[/red]", "bold red"),
    }
    n_warn = 0
    n_miss = 0
    for c in checks:
        badge, detail_style = style_for.get(c.status, ("?", ""))
        if c.status == diagnostics.WARN:
            n_warn += 1
        elif c.status == diagnostics.MISSING:
            n_miss += 1
        table.add_row(
            badge,
            c.label if c.core else f"{c.label} [dim](optional)[/dim]",
            f"[{detail_style}]{c.detail}[/{detail_style}]" if detail_style else c.detail,
            c.purpose,
        )

    console.print(table)

    core_failure = diagnostics.has_core_failure(checks)
    summary = f"{n_miss} missing, {n_warn} warnings"
    if core_failure:
        console.print(f"[red]{summary} — core pipeline components missing.[/red]")
    else:
        console.print(f"[green]{summary} — core pipeline OK.[/green]")

    if strict and core_failure:
        raise typer.Exit(code=1)


@app.command(name="health")
def health_command(
    db: Path = typer.Option(Path(DEFAULT_DB), "--db", "-d", help="Dataset DuckDB (opened read-only)."),
    output_format: BriefFormat = typer.Option(BriefFormat.markdown, "--format", "-f", help="markdown or json"),
    strict: bool = typer.Option(False, "--strict", help="Exit non-zero when any check fails."),
):
    """Check a dataset for problems that distort analyses: missing gene calls, unplaced genes,
    annotation gaps, completeness, stale labels, seal drift."""
    from sharur.health import run_checks  # noqa: PLC0415

    if not db.is_file():
        typer.echo(f"DuckDB file does not exist: {db}", err=True)
        raise typer.Exit(1)
    report = run_checks(db)
    if output_format == BriefFormat.json:
        typer.echo(json.dumps(report.to_dict(), indent=2, default=str))
    else:
        typer.echo(report.to_markdown())
    if strict and report.worst == "fail":
        raise typer.Exit(code=1)


@app.command("adopt-astra-hmms")
def adopt_astra_hmms(
    dry_run: bool = typer.Option(False, "--dry-run", help="List what would be adopted; write nothing."),
):
    """Register databases installed through Astra in Aksha's registry.

    Aksha shares Astra's database storage but keeps its own registry, so
    Astra-installed databases stay on disk and unknown to `aksha search`.
    This marks each one installed at its existing directory (no download);
    the previous registry is kept as hmm_databases.json.bak.
    """
    from sharur import hmm_search

    try:
        adopted = hmm_search.adopt_astra_registry(dry_run=dry_run)
    except FileNotFoundError as exc:
        typer.echo(str(exc), err=True)
        raise typer.Exit(code=1)
    verb = "would adopt" if dry_run else "adopted"
    for name, directory in adopted:
        typer.echo(f"{verb} {name}: {directory}")
    typer.echo(f"{len(adopted)} database(s) {verb} into {hmm_search.registry_path('aksha')}")


@app.command(name="architecture")
def architecture_search(
    pattern: str = typer.Argument(..., help='Domain pattern, e.g. "Big_2{5,} . VWA" or "^ SP ( Cadherin | Big_2 )+"'),
    db: str = typer.Option(DEFAULT_DB, "--db", "-d", help="Path to DuckDB database"),
    source: list[str] = typer.Option(["pfam"], "--source", "-s", help="Annotation source(s); repeatable."),
    bins: list[str] = typer.Option([], "--bin", "-b", help="Restrict to genome(s); repeatable."),
    limit: int = typer.Option(50, "--limit", "-n", help="Records to show (the total is always counted)."),
    max_overlap: float = typer.Option(0.5, "--max-overlap", help="Overlap (fraction of the shorter domain) "
                                      "above which the weaker hit is dropped."),
    output_format: BriefFormat = typer.Option(BriefFormat.markdown, "--format", "-f", help="markdown or json"),
):
    """Find proteins whose ordered domains match a pattern.

    Tokens: domain names or accessions, globs (Big_*), '.' for any domain,
    quantifiers (? * + {m} {m,} {m,n}; put a space before ? and *), groups
    with alternatives ( A | B ), and anchors ^ (N-terminus) and $ (C-terminus).
    """
    from sharur.architecture import PatternError, architecture_markdown, search_architecture
    from sharur.storage.duckdb_store import DuckDBStore

    try:
        with DuckDBStore(db, read_only=True) as store:
            result = search_architecture(store, pattern, sources=tuple(source), bins=bins or None,
                                         limit=limit, max_overlap=max_overlap)
    except PatternError as exc:
        typer.echo(str(exc), err=True)
        raise typer.Exit(1) from exc
    typer.echo(json.dumps(result, indent=2, default=str) if output_format == BriefFormat.json
               else architecture_markdown(result))


@app.command(name="import-coverage")
def import_coverage_command(
    paths: list[Path] = typer.Argument(..., exists=True, dir_okay=False, help="CoverM contig table(s) or long TSV/CSV."),
    db: Path = typer.Option(Path(DEFAULT_DB), "--db", "-d", help="Core DuckDB (validation and sidecar location)."),
    fmt: str = typer.Option("coverm", "--format", help="coverm (contig table) or long (sample_id, contig_id, ...)."),
    samples_file: Path | None = typer.Option(None, "--samples", help="TSV/CSV with sample_id and metadata columns."),
    sidecar: Path | None = typer.Option(None, "--sidecar", help="Default: abundance.duckdb beside --db."),
    allow_unknown: bool = typer.Option(False, "--allow-unknown-contigs",
                                       help="Keep matching rows when some contigs are absent from the dataset."),
):
    """Import per-sample contig coverage into the abundance sidecar (the core database is unchanged)."""
    from sharur.abundance import default_path, import_coverage
    from sharur.storage.duckdb_store import DuckDBStore

    if fmt not in ("coverm", "long"):
        typer.echo("--format must be coverm or long", err=True)
        raise typer.Exit(1)
    target = sidecar or default_path(db)
    try:
        with DuckDBStore(str(db), read_only=True) as store:
            result = import_coverage(store, target, list(paths), fmt=fmt, sample_metadata=samples_file,
                                     allow_unknown_contigs=allow_unknown)
    except ValueError as exc:
        typer.echo(str(exc), err=True)
        raise typer.Exit(1) from exc
    typer.echo(f"Imported {result['rows']:,} contig rows for {len(result['samples'])} samples into {target}"
               + (f" ({result['unknown_contigs']:,} unknown contigs skipped)" if result["unknown_contigs"] else ""))


@app.command(name="abundance")
def abundance_command(
    db: Path = typer.Option(Path(DEFAULT_DB), "--db", "-d", help="Core DuckDB."),
    sample: list[str] = typer.Option([], "--sample", help="Sample(s); repeatable (default: all)."),
    bins: list[str] = typer.Option([], "--bin", "-b", help="Genome(s); repeatable."),
    predicate: str | None = typer.Option(None, "--predicate", "-p",
                                         help="Share of reads in genomes carrying this predicate."),
    annotation: str | None = typer.Option(None, "--annotation", "-a",
                                          help="Share of reads in genomes carrying this accession or name."),
    top: int = typer.Option(10, "--top", "-n", help="Genomes shown per sample (markdown)."),
    sidecar: Path | None = typer.Option(None, "--sidecar", help="Default: abundance.duckdb beside --db."),
    output_format: BriefFormat = typer.Option(BriefFormat.markdown, "--format", "-f", help="markdown or json"),
):
    """Genome abundance per sample, or the share of reads in genomes carrying a feature."""
    from sharur.abundance import (
        default_path,
        feature_abundance,
        feature_abundance_markdown,
        genome_abundance,
        genome_abundance_markdown,
    )
    from sharur.storage.duckdb_store import DuckDBStore

    target = sidecar or default_path(db)
    try:
        if predicate or annotation:
            with DuckDBStore(str(db), read_only=True) as store:
                result = feature_abundance(store, target, predicate=predicate, annotation=annotation,
                                           sample_ids=sample or None)
            text = feature_abundance_markdown(result)
        else:
            result = genome_abundance(target, sample_ids=sample or None, bins=bins or None)
            text = genome_abundance_markdown(result, top=top)
    except (FileNotFoundError, ValueError) as exc:
        typer.echo(str(exc), err=True)
        raise typer.Exit(1) from exc
    typer.echo(json.dumps(result, indent=2, default=str) if output_format == BriefFormat.json else text)


@app.command(name="coverage-outliers")
def coverage_outliers_command(
    bin_id: str = typer.Argument(..., help="Genome (bin) ID."),
    db: Path = typer.Option(Path(DEFAULT_DB), "--db", "-d", help="Core DuckDB."),
    min_log2: float = typer.Option(1.0, "--min-log2", help="Minimum |log2(contig depth / genome median)|."),
    sidecar: Path | None = typer.Option(None, "--sidecar", help="Default: abundance.duckdb beside --db."),
    output_format: BriefFormat = typer.Option(BriefFormat.markdown, "--format", "-f", help="markdown or json"),
):
    """Contigs whose depth departs from their genome's median in most samples (possible binning errors)."""
    from sharur.abundance import coverage_outliers, coverage_outliers_markdown, default_path

    try:
        result = coverage_outliers(sidecar or default_path(db), bin_id, min_log2=min_log2)
    except FileNotFoundError as exc:
        typer.echo(str(exc), err=True)
        raise typer.Exit(1) from exc
    typer.echo(json.dumps(result, indent=2, default=str) if output_format == BriefFormat.json
               else coverage_outliers_markdown(result))


@app.command(name="browse")
def browse_command(
    db: Path = typer.Option(Path(DEFAULT_DB), "--db", "-d", help="Dataset DuckDB (opened read-only)."),
    host: str = typer.Option("127.0.0.1", "--host", help="Interface to bind; 127.0.0.1 keeps it on this machine."),
    port: int = typer.Option(8800, "--port", "-p"),
    share: bool = typer.Option(False, "--share", help="Require a random access token (printed in the link)."),
    notes: Path | None = typer.Option(None, "--notes", help="Notes and flags database "
                                      "(default: browser_notes.sqlite beside the dataset)."),
    assemblies: list[Path] = typer.Option([], "--assemblies", help="Directory of genome FASTAs named BIN_ID.fna "
                                          "(for CRISPR repeats and spacers); repeatable. Dataset folders are searched too."),
    summary_cache: bool = typer.Option(True, "--summary-cache/--no-summary-cache",
                                       help="Compute startup summaries once per database version in a worker "
                                            "process and reuse them (keeps the server near 1 GB)."),
):
    """Browse a dataset in a web browser: protein cards, neighborhoods, evidence, genomes, searches."""
    import secrets  # noqa: PLC0415

    import uvicorn  # noqa: PLC0415

    from sharur.browser import create_app  # noqa: PLC0415

    if not db.is_file():
        typer.echo(f"DuckDB file does not exist: {db}", err=True)
        raise typer.Exit(1)
    loopback = host in ("127.0.0.1", "localhost", "::1")
    token = secrets.token_urlsafe(16) if share or not loopback else None
    shown_host = "localhost" if loopback else host
    typer.echo(f"Sharur browser: http://{shown_host}:{port}/" + (f"?token={token}" if token else ""))
    if not loopback and not share:
        typer.echo("Binding beyond this machine: access requires the token in the link above.")
    from sharur.browser.startup_cache import default_dir  # noqa: PLC0415

    app = create_app(db, token=token, notes_path=notes, assemblies=assemblies,
                     cache_dir=default_dir() if summary_cache else None)
    uvicorn.run(app, host=host, port=port,
                log_level="warning")


@app.command(name="describe")
def describe(
    db: str = typer.Option(DEFAULT_DB, "--db", "-d", help="Path to DuckDB database"),
    output_format: BriefFormat = typer.Option(BriefFormat.markdown, "--format", "-f", help="markdown or json"),
):
    """What a dataset holds: annotation sources, curated callers, genome metadata, predicate state."""
    from sharur.operators.introspection import describe_dataset, describe_dataset_markdown
    from sharur.storage.duckdb_store import DuckDBStore

    with DuckDBStore(db, read_only=True) as store:
        result = describe_dataset(store)
    typer.echo(json.dumps(result, indent=2, default=str) if output_format == BriefFormat.json
               else describe_dataset_markdown(result))


@app.command(name="card")
def protein_card(
    protein_id: str = typer.Argument(..., help="Protein ID"),
    db: str = typer.Option(DEFAULT_DB, "--db", "-d", help="Path to DuckDB database"),
    window: int = typer.Option(5, "--window", "-w", help="Neighbors on each side."),
    output_format: BriefFormat = typer.Option(BriefFormat.markdown, "--format", "-f", help="markdown or json"),
):
    """Summarize one protein: context, annotations, evidence-backed predicates (no sequences)."""
    from sharur.operators.cards import card, card_markdown
    from sharur.storage.duckdb_store import DuckDBStore

    with DuckDBStore(db, read_only=True) as store:
        result = card(store, protein_id, window=window)
    typer.echo(json.dumps(result, indent=2, default=str) if output_format == BriefFormat.json
               else card_markdown(result))


@app.command(name="why")
def why_predicate(
    protein_id: str = typer.Argument(..., help="Protein ID"),
    predicate: str = typer.Argument(..., help="Predicate ID, e.g. nad_binding"),
    db: str = typer.Option(DEFAULT_DB, "--db", "-d", help="Path to DuckDB database"),
    output_format: BriefFormat = typer.Option(BriefFormat.markdown, "--format", "-f", help="markdown or json"),
):
    """Explain why a protein carries a predicate: annotation hits, map evidence, expansion chain."""
    from sharur.operators.cards import why, why_markdown
    from sharur.storage.duckdb_store import DuckDBStore

    with DuckDBStore(db, read_only=True) as store:
        result = why(store, protein_id, predicate)
    typer.echo(json.dumps(result, indent=2, default=str) if output_format == BriefFormat.json
               else why_markdown(result))


@app.command(name="modules")
def kegg_modules(
    db: str = typer.Option(DEFAULT_DB, "--db", "-d", help="Path to DuckDB database"),
    bins: Optional[list[str]] = typer.Option(None, "--bin", "-b", help="Genome(s); repeatable."),
    modules: Optional[list[str]] = typer.Option(None, "--module", "-m", help="Module(s), e.g. M00175; repeatable."),
    min_completeness: float = typer.Option(0.0, "--min-completeness", help="Report modules at or above this fraction."),
    around: Optional[str] = typer.Option(None, "--around", help="Protein ID: modules encoded in its neighborhood."),
    window: int = typer.Option(10, "--window", "-w", help="Genes on each side for --around."),
    output_format: BriefFormat = typer.Option(BriefFormat.markdown, "--format", "-f", help="markdown or json"),
):
    """KEGG module completeness per genome, or around one protein (needs `sharur setup-kegg`)."""
    from sharur.modules import genome_modules, locus_modules, modules_markdown
    from sharur.storage.duckdb_store import DuckDBStore

    with DuckDBStore(db, read_only=True) as store:
        if around:
            rows = locus_modules(store, around, window=window)
            title = f"KEGG modules within {window} genes of {around}"
        else:
            rows = genome_modules(store, bins=bins, modules=modules, min_completeness=min_completeness)
            title = "KEGG module completeness"
    typer.echo(json.dumps(rows, indent=2, default=str) if output_format == BriefFormat.json
               else modules_markdown(rows, title))


@app.command(name="setup-kegg")
def setup_kegg(
    directory: Optional[Path] = typer.Option(
        None, "--dir", help="Where to build (default: $SHARUR_KEGG_DIR, else data/reference/kegg)."
    ),
    inputs: Optional[Path] = typer.Option(
        None, "--inputs",
        help="Build from KEGG files already on disk (list/ko as ko_list.tsv, brite/<id>.json, "
             "modules.txt) instead of fetching, e.g. from a licensed KEGG copy.",
    ),
    kofam_ko_list: Path = typer.Option(
        Path("data/reference/ko_list"), "--kofam-ko-list",
        help="KOfam ko_list (downloaded when absent).",
    ),
    report: Optional[Path] = typer.Option(None, "--report", help="TSV of proposed pairs without evidence."),
):
    """Fetch KEGG data and build the KO -> predicate map on this machine.

    KEGG data is subject to KEGG's terms: KEGG REST (used by default, about 100
    requests at under 3 per second) is for academic use; other users need a KEGG
    license and can build from their licensed files with --inputs. The built
    files stay local; Sharur ships only its own rules.
    """
    from sharur.predicates.mappings.kegg_build import setup
    from sharur.predicates.mappings.kegg_map import default_kegg_dir

    out = directory or default_kegg_dir()
    typer.echo("KEGG data is fetched under KEGG's terms of use (https://www.kegg.jp/kegg/legal.html): "
               "KEGG REST is for academic use; non-academic users need a KEGG license.")
    provenance = setup(out, inputs=inputs, kofam_ko_list=kofam_ko_list, report=report, log=typer.echo)
    typer.echo(f"Built {out / 'kegg_predicates.tsv'} (KEGG {provenance['kegg_release']}; "
               f"map sha256 {provenance['map_sha256'][:16]})")
    typer.echo("Regenerate predicates for existing datasets to apply it.")


@app.command(name="import-quality")
def import_quality_command(
    table: Path = typer.Argument(..., help="CheckM2 quality_report.tsv, CheckM --tab_table output, or GTDB "
                                 "ar53/bac120 metadata (.tsv or .tsv.gz)."),
    db: Path = typer.Option(Path(DEFAULT_DB), "--db", "-d", help="Dataset DuckDB."),
    id_column: str | None = typer.Option(None, "--id-column", help="Genome column (default: detected)."),
    completeness_column: str | None = typer.Option(None, "--completeness-column",
                                                   help="Completeness column, percent (default: detected; "
                                                        "CheckM2 preferred)."),
    contamination_column: str | None = typer.Option(None, "--contamination-column"),
    overwrite: bool = typer.Option(False, "--overwrite", help="Replace completeness values already stored."),
    dry_run: bool = typer.Option(False, "--dry-run", help="Report matches and leave the dataset unchanged."),
):
    """Fill genome completeness and contamination from a CheckM, CheckM2 or GTDB table."""
    import duckdb  # noqa: PLC0415

    from sharur.quality import import_quality, read_quality_table  # noqa: PLC0415

    if not db.is_file():
        typer.echo(f"DuckDB file does not exist: {db}", err=True)
        raise typer.Exit(1)
    try:
        values, columns = read_quality_table(table, id_column=id_column, completeness_column=completeness_column,
                                             contamination_column=contamination_column)
    except ValueError as exc:
        typer.echo(str(exc), err=True)
        raise typer.Exit(1) from exc
    conn = duckdb.connect(str(db), read_only=dry_run)
    try:
        result = import_quality(conn, values, overwrite=overwrite, dry_run=dry_run)
    finally:
        conn.close()
    typer.echo(f"{table.name}: {len(values):,} genomes (columns {columns['id']}, {columns['completeness']}"
               + (f", {columns['contamination']}" if columns["contamination"] else "") + ")")
    typer.echo(f"{result['matched']:,} of {result['bins']:,} dataset genomes matched; "
               f"{result['updated']:,} {'would be ' if dry_run else ''}updated"
               + (f" (median completeness {result['median_completeness']:.1f}%)" if result["median_completeness"] else ""))
    if not dry_run and result["updated"]:
        typer.echo("Re-seal the dataset when its writes are done.")


@app.command(name="cas-type")
def cas_type(
    db: Path = typer.Option(Path(DEFAULT_DB), "--db", "-d", help="Dataset DuckDB."),
    cctyper_db: Path | None = typer.Option(None, "--cctyper-db", help="CCTyper database directory "
                                           "(default: $SHARUR_CCTYPER_DB, then Aksha's CCTyper install)."),
    workers: int = typer.Option(0, "--workers", "-w", help="Parallel genome searches (default: all CPUs)."),
    genome: list[str] = typer.Option([], "--genome", "-g", help="Type only these genomes; repeatable."),
    dry_run: bool = typer.Option(False, "--dry-run", help="Report the calls and leave the dataset unchanged."),
):
    """Subtype CRISPR-Cas systems with CRISPRCasTyper's profiles and scoring.

    Searches the CCTyper Cas profiles against every genome's proteins, groups
    hits into operons, scores subtypes, classifies array repeats and links
    operons to arrays within 10 kb. Calls go to ``crispr_cas_systems``,
    ``crispr_array_types`` and ``system_proteins`` (source ``cctyper``),
    replacing earlier calls from this caller.
    """
    from collections import Counter  # noqa: PLC0415

    from sharur.cas_typing import run_dataset, system_rows, write_results  # noqa: PLC0415

    if not db.is_file():
        typer.echo(f"DuckDB file does not exist: {db}", err=True)
        raise typer.Exit(1)

    def progress(done: int, total: int) -> None:
        if done == total or done % 250 == 0:
            typer.echo(f"  searched {done}/{total} genomes", err=True)

    try:
        result = run_dataset(db, cctyper_db=cctyper_db, workers=workers or None, genomes=genome or None,
                             progress=progress)
    except FileNotFoundError as exc:
        typer.echo(str(exc), err=True)
        raise typer.Exit(1) from exc
    systems = result["systems"]
    status = Counter(op["status"] for op in systems)
    typer.echo(f"{result['genomes']} genomes: {len(systems)} Cas operons "
               + ", ".join(f"{k} {v}" for k, v in sorted(status.items())))
    named = Counter(op["joint_prediction"] for op in systems if op["status"] in ("crispr_cas", "cas"))
    if named:
        typer.echo("Subtypes: " + ", ".join(f"{k} {v}" for k, v in named.most_common()))
    arrays = result["arrays"]
    typer.echo(f"{len(arrays)} arrays: {sum(a['prediction'] != 'Unknown' for a in arrays)} with a repeat subtype, "
               f"{sum(a['near_cas'] for a in arrays)} linked to Cas operons")
    if dry_run:
        system_rows(systems)
        return
    counts = write_results(db, result)
    typer.echo(f"Wrote {counts['systems']} systems, {counts['arrays']} array types and "
               f"{counts['members']} member proteins to {db}. Re-seal the dataset when its writes are done.")


# ------------------------------------------------------------------ #
# Main entry point
# ------------------------------------------------------------------ #


def main():
    app(prog_name="sharur")


if __name__ == "__main__":
    main()
