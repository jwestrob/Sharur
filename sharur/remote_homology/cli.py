"""``sharur remote-homology``: CDD per protein and HMM-HMM searches per family ("resembles" rows).

Typical use (``docs/guides/remote-homology.md``)::

    sharur remote-homology check --config rh.json
    # family tier
    sharur remote-homology export --db DB [--families FAMILIES.duckdb --run RUN] --out-dir WORK
    sharur remote-homology search families --work WORK --config rh.json        # or the Slurm shards
    sharur remote-homology load-families --work WORK --config rh.json --staging RH.duckdb
    # CDD
    sharur remote-homology export-cdd --db DB --out-dir WORK
    sharur remote-homology search cdd --work WORK --config rh.json
    sharur remote-homology compact-cdd --work WORK                 # beside the raw tables
    sharur remote-homology load-cdd --work WORK --config rh.json --staging RH.duckdb
    # into the dataset database (staged write)
    sharur remote-homology activate --db DB --staging RH.duckdb [--families-staging FAMILIES.duckdb] --publish
    sharur remote-homology status --db DB
"""

from __future__ import annotations

import json
from pathlib import Path

import typer


remote_homology_app = typer.Typer(
    no_args_is_help=True,
    help="Remote-homology annotation: CDD per protein and HHsearch per protein family (\"resembles\" rows).",
)

_DB_HELP = "Dataset sharur.duckdb."


def _echo_json(payload) -> None:
    typer.echo(json.dumps(payload, indent=2, default=str))


def _fail(label: str, exc: Exception) -> None:
    typer.echo(f"{label}: {exc}", err=True)
    raise typer.Exit(1) from exc


@remote_homology_app.command("check")
def check_command(config: Path = typer.Option(..., "--config", "-c", help="Runner configuration (JSON).")):
    """Check that the configured tools and databases are present; print tool versions."""
    from sharur.remote_homology import runner  # noqa: PLC0415

    try:
        _echo_json(runner.check(runner.load_config(config)))
    except runner.ConfigError as exc:
        _fail("check", exc)


@remote_homology_app.command("export")
def export_command(
    db: Path = typer.Option(..., "--db", "-d", help=_DB_HELP + " Opened read-only."),
    out_dir: Path = typer.Option(..., "--out-dir", "-o", help="Work directory for the family tier."),
    families: Path | None = typer.Option(None, "--families", "-f",
                                            help="Staging file with the clustering run (default: the dataset)."),
    run: str | None = typer.Option(None, "--run", "-r", help="Clustering run (default: the default family run)."),
    min_genomes: int = typer.Option(3, "--min-genomes", help="Families found in at least this many genomes."),
    include_pfam_annotated: bool = typer.Option(False, "--include-pfam-annotated",
                                                help="Also select families whose members carry Pfam hits."),
    max_sequences: int = typer.Option(250, "--max-sequences", help="Unique sequences per family alignment."),
    family_id: list[str] = typer.Option([], "--family-id", help="Export these families only (repeatable)."),
    duckdb_threads: int = typer.Option(4, "--duckdb-threads"),
    memory_limit: str = typer.Option("8GB", "--memory-limit"),
):
    """Select families (default: Pfam-dark, >= 3 genomes) and write one FASTA per family."""
    from sharur.remote_homology import store  # noqa: PLC0415

    try:
        result = store.export_families(db, out_dir, families_db=families, run=run, min_genomes=min_genomes,
                                       require_pfam_dark=not include_pfam_annotated, max_sequences=max_sequences,
                                       family_ids=family_id or None, threads=duckdb_threads,
                                       memory_limit=memory_limit, log=lambda m: typer.echo(m, err=True))
    except store.RemoteHomologyError as exc:
        _fail("export", exc)
    _echo_json(result)


@remote_homology_app.command("export-cdd")
def export_cdd_command(
    db: Path = typer.Option(..., "--db", "-d", help=_DB_HELP + " Opened read-only."),
    out_dir: Path = typer.Option(..., "--out-dir", "-o", help="Work directory (cdd/ is created inside)."),
    chunk_size: int = typer.Option(2000, "--chunk-size", help="Unique sequences per RPS-BLAST chunk."),
    duckdb_threads: int = typer.Option(4, "--duckdb-threads"),
    memory_limit: str = typer.Option("8GB", "--memory-limit"),
):
    """Write every unique protein sequence (named by MD5) in chunks for RPS-BLAST."""
    from sharur.remote_homology import store  # noqa: PLC0415

    _echo_json(store.export_cdd(db, out_dir, chunk_size=chunk_size, threads=duckdb_threads,
                                memory_limit=memory_limit, log=lambda m: typer.echo(m, err=True)))


@remote_homology_app.command("search")
def search_command(
    stage: str = typer.Argument(..., help="families or cdd"),
    work: Path = typer.Option(..., "--work", "-w"),
    config: Path = typer.Option(..., "--config", "-c"),
    threads: int | None = typer.Option(None, "--threads", "-t",
                                          help="Default: $SLURM_CPUS_ON_NODE or all CPUs."),
    threads_per_worker: int = typer.Option(4, "--threads-per-worker"),
    shard: int = typer.Option(0, "--shard"),
    n_shards: int = typer.Option(1, "--n-shards"),
    limit: int | None = typer.Option(None, "--limit"),
):
    """Run the searches on this host (the Slurm wrapper runs the same code per shard)."""
    from sharur.remote_homology import runner  # noqa: PLC0415

    if stage not in ("families", "cdd"):
        raise typer.BadParameter("stage must be 'families' or 'cdd'")
    argv = [stage, "--work", str(work), "--config", str(config), "--shard", str(shard), "--n-shards",
            str(n_shards), "--threads-per-worker", str(threads_per_worker)]
    if threads:
        argv += ["--threads", str(threads)]
    if limit:
        argv += ["--limit", str(limit)]
    raise typer.Exit(runner.main(argv))


@remote_homology_app.command("load-families")
def load_families_command(
    work: Path = typer.Option(..., "--work", "-w"),
    config: Path = typer.Option(..., "--config", "-c"),
    staging: Path = typer.Option(..., "--staging", "-s", help="Staging DuckDB (created when absent)."),
    allow_partial: bool = typer.Option(False, "--allow-partial",
                                       help="Load although some families are pending or failed."),
):
    """Parse the family-tier results into family_annotations (staging file)."""
    from sharur.remote_homology import store  # noqa: PLC0415

    try:
        _echo_json(store.load_families(work, staging, config, allow_partial=allow_partial,
                                       log=lambda m: typer.echo(m, err=True)))
    except store.RemoteHomologyError as exc:
        _fail("load-families", exc)


@remote_homology_app.command("compact-cdd")
def compact_cdd_command(
    work: Path = typer.Option(..., "--work", "-w"),
    max_evalue: float = typer.Option(0.01, "--max-evalue"),
    allow_partial: bool = typer.Option(False, "--allow-partial"),
    threads: int | None = typer.Option(None, "--threads", "-t", help="Default: $SLURM_CPUS_ON_NODE or 4."),
):
    """Reduce the raw RPS-BLAST tables to work/cdd/compact/*.parquet (run beside them, on the cluster)."""
    from sharur.remote_homology import runner, store  # noqa: PLC0415

    try:
        _echo_json(store.compact_cdd(work, max_evalue=max_evalue, allow_partial=allow_partial,
                                     threads=threads or min(runner.allocated_cpus(), 4),
                                     log=lambda m: typer.echo(m, err=True)))
    except store.RemoteHomologyError as exc:
        _fail("compact-cdd", exc)


@remote_homology_app.command("load-cdd")
def load_cdd_command(
    work: Path = typer.Option(..., "--work", "-w"),
    config: Path = typer.Option(..., "--config", "-c"),
    staging: Path = typer.Option(..., "--staging", "-s", help="Staging DuckDB (created when absent)."),
    max_evalue: float = typer.Option(0.01, "--max-evalue"),
    max_overlap: float = typer.Option(0.5, "--max-overlap",
                                      help="Best model per region: hits sharing more than this fraction of the "
                                           "shorter range belong to one region (best bit score kept)."),
    allow_partial: bool = typer.Option(False, "--allow-partial"),
):
    """Load the CDD hits (per unique sequence) and models into a staging file; compacts first if needed."""
    from sharur.remote_homology import store  # noqa: PLC0415

    try:
        _echo_json(store.load_cdd(work, staging, config, max_evalue=max_evalue, max_overlap=max_overlap,
                                  allow_partial=allow_partial, log=lambda m: typer.echo(m, err=True)))
    except store.RemoteHomologyError as exc:
        _fail("load-cdd", exc)


@remote_homology_app.command("compact-staging")
def compact_staging_command(
    staging: Path = typer.Option(..., "--staging", "-s", help="Remote-homology staging file to rewrite."),
    duckdb_threads: int = typer.Option(4, "--duckdb-threads"),
    memory_limit: str = typer.Option("8GB", "--memory-limit"),
):
    """Rewrite a staging file without free blocks or legacy indexes (checked copy, atomic replace)."""
    from sharur.remote_homology import store  # noqa: PLC0415

    try:
        _echo_json(store.compact_staging_file(staging, threads=duckdb_threads, memory_limit=memory_limit))
    except store.RemoteHomologyError as exc:
        _fail("compact-staging", exc)


@remote_homology_app.command("activate")
def activate_command(
    db: Path = typer.Option(..., "--db", "-d", help=_DB_HELP),
    staging: list[Path] = typer.Option(..., "--staging", "-s",
                                       help="Remote-homology staging file(s) (repeatable)."),
    families_staging: list[Path] = typer.Option(
        [], "--families-staging",
        help="Protein-family staging file(s), loaded first in the order given (repeatable; a later file's "
             "run replaces a same-named one)."),
    family_runs: str = typer.Option(
        "searched", "--family-runs",
        help="Clustering runs taken from the families staging files: searched (those the remote-homology "
             "rows refer to) or all (every run, e.g. a dataset's own build and imported clusterings)."),
    label: str = typer.Option("remote-homology", "--label", "-l"),
    work_dir: Path | None = typer.Option(None, "--work-dir",
                                            help="Directory for the staged clone (default: beside DB)."),
    do_publish: bool = typer.Option(False, "--publish", help="Publish the staged copy when the load succeeds."),
    verify_digest: bool = typer.Option(True, "--verify-digest/--no-verify-digest"),
    cdd_rows: str = typer.Option("best_per_region", "--cdd-rows",
                                 help="CDD rows written to annotations: best_per_region (the best model of each "
                                      "region) or all. The full hit table stays in the staging file."),
):
    """Stage DB, load the remote-homology rows into the clone, diff, and (with --publish) swap it in.

    Follows the staged-write protocol (docs/guides/staged-writes.md), as `sharur families
    activate` does: the canonical file changes only at publish, a backup is kept and the
    dataset is resealed.
    """
    from sharur.protein_families.cli import _staged  # noqa: PLC0415
    from sharur.remote_homology import store  # noqa: PLC0415

    if cdd_rows not in store.CDD_ROW_CHOICES:
        raise typer.BadParameter(f"--cdd-rows must be one of {', '.join(store.CDD_ROW_CHOICES)}")

    if family_runs not in (store.FAMILY_RUNS_SEARCHED, store.FAMILY_RUNS_ALL):
        raise typer.BadParameter("--family-runs must be 'searched' or 'all'")
    paths = [p.expanduser().resolve() for p in staging]
    fam = [f.expanduser().resolve() for f in families_staging]

    def writer(path: Path) -> dict:
        con = store.connect(path)
        try:
            return store.load_staging(con, paths, families_staging=fam, family_runs=family_runs,
                                      verify_digest=verify_digest, cdd_rows=cdd_rows,
                                      log=lambda m: typer.echo(m, err=True))
        finally:
            con.close()

    _staged(db.expanduser().resolve(), label, work_dir, writer, do_publish, tables=store.ACTIVATION_TABLES)


@remote_homology_app.command("status")
def status_command(db: Path = typer.Option(..., "--db", "-d", help=_DB_HELP + " Opened read-only.")):
    """Searches loaded, families resembling characterized families, CDD coverage."""
    from sharur.remote_homology import store  # noqa: PLC0415

    con = store.connect(db, read_only=True)
    try:
        _echo_json(store.status(con))
    finally:
        con.close()


__all__ = ["remote_homology_app"]
