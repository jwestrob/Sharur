"""``sharur families``: build, activate, extend and inspect the protein-family layer.

Typical use (``docs/guides/protein-families.md``)::

    sharur families run --db data/DS/sharur.duckdb --out-dir data/DS/stage09_families
    sharur families activate --db data/DS/sharur.duckdb \\
        --staging data/DS/stage09_families/protein_families.duckdb --publish
    sharur families assign --db data/DS/sharur.duckdb --publish     # after adding genomes
    sharur families status --db data/DS/sharur.duckdb
    sharur families import --db data/DS/sharur.duckdb --tsv clu.tsv --name RUN \
        --parameters "mmseqs cluster ..." --staging data/DS/stage09_families/protein_families.duckdb
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

import typer

from sharur.protein_families import mmseqs

families_app = typer.Typer(
    no_args_is_help=True,
    help="Annotation-independent protein families (MMseqs2): build, import, activate, assign, status.",
)

_DB_HELP = "Dataset sharur.duckdb."


def _echo_json(payload) -> None:
    typer.echo(json.dumps(payload, indent=2, default=str))


def _levels(subfamily: float, family: float, coverage: float, sensitivity: float) -> list[mmseqs.LevelParams]:
    try:
        return mmseqs.default_levels(subfamily, family, coverage=coverage, sensitivity=sensitivity)
    except ValueError as exc:
        raise typer.BadParameter(str(exc)) from exc


@families_app.command("run")
def run_command(
    db: Path = typer.Option(..., "--db", "-d", help=_DB_HELP + " Opened read-only."),
    out_dir: Path = typer.Option(..., "--out-dir", "-o", help="Work directory (e.g. DATASET/stage09_families)."),
    subfamily_min_seq_id: float = typer.Option(mmseqs.SUBFAMILY_MIN_SEQ_ID, "--subfamily-min-seq-id"),
    family_min_seq_id: float = typer.Option(mmseqs.FAMILY_MIN_SEQ_ID, "--family-min-seq-id"),
    coverage: float = typer.Option(mmseqs.COVERAGE, "--coverage", help="Coverage of both sequences (-c, cov-mode 0)."),
    sensitivity: float = typer.Option(mmseqs.SENSITIVITY, "--sensitivity", "-s"),
    threads: Optional[int] = typer.Option(None, "--threads", "-t",
                                          help="MMseqs2 threads (default: $SLURM_CPUS_ON_NODE or all CPUs)."),
    duckdb_threads: Optional[int] = typer.Option(None, "--duckdb-threads", help="Default: --threads."),
    memory_limit: Optional[str] = typer.Option(None, "--memory-limit",
                                               help="DuckDB memory limit (default: 25% of the allocation)."),
    split_memory_limit: Optional[str] = typer.Option(None, "--split-memory-limit",
                                                     help="MMseqs2 memory cap (default: 70% of the allocation)."),
    keep_fasta: bool = typer.Option(False, "--keep-fasta"),
    mmseqs_binary: Optional[str] = typer.Option(None, "--mmseqs", help="MMseqs2 executable."),
):
    """Export proteins, cluster subfamilies and families, write a staging DuckDB (the DB is untouched)."""
    from sharur.protein_families import store  # noqa: PLC0415

    try:
        result = store.run_pipeline(
            db, out_dir, levels=_levels(subfamily_min_seq_id, family_min_seq_id, coverage, sensitivity),
            threads=threads, duckdb_threads=duckdb_threads, memory_limit=memory_limit,
            split_memory_limit=split_memory_limit, keep_fasta=keep_fasta, binary=mmseqs_binary,
            log=lambda m: typer.echo(m, err=True))
    except store.FamilyError as exc:
        typer.echo(f"run: {exc}", err=True)
        raise typer.Exit(1) from exc
    _echo_json(result)
    typer.echo(f"Activate with: sharur families activate --db {db} --staging {result['staging']} --publish",
               err=True)


@families_app.command("export")
def export_command(
    db: Path = typer.Option(..., "--db", "-d", help=_DB_HELP + " Opened read-only."),
    out_dir: Path = typer.Option(..., "--out-dir", "-o"),
    duckdb_threads: int = typer.Option(4, "--duckdb-threads"),
    memory_limit: str = typer.Option("8GB", "--memory-limit"),
):
    """Write proteins.faa (protein_id order, trailing '*' removed) and export.json (count, digest)."""
    from sharur.protein_families import store  # noqa: PLC0415

    _echo_json(store.write_export(db.expanduser().resolve(), out_dir, threads=duckdb_threads,
                                  memory_limit=memory_limit))


@families_app.command("build")
def build_command(
    db: Path = typer.Option(..., "--db", "-d", help=_DB_HELP + " Opened read-only."),
    clustering: Path = typer.Option(..., "--clustering", "-c",
                                    help="Directory with clustering.json, level0/1.tsv.gz and export.json."),
    staging: Path = typer.Option(..., "--staging", "-s", help="Staging DuckDB to write (replaced)."),
    duckdb_threads: int = typer.Option(4, "--duckdb-threads"),
    memory_limit: str = typer.Option("8GB", "--memory-limit"),
):
    """Build the staging DuckDB from a finished clustering, such as one run on a cluster node."""
    from sharur.protein_families import store  # noqa: PLC0415

    try:
        result = store.build_staging(db.expanduser().resolve(), clustering, staging, threads=duckdb_threads,
                                     memory_limit=memory_limit, log=lambda m: typer.echo(m, err=True))
    except store.FamilyError as exc:
        typer.echo(f"build: {exc}", err=True)
        raise typer.Exit(1) from exc
    _echo_json(result)


@families_app.command("load")
def load_command(
    db: Path = typer.Option(..., "--db", "-d", help="Writable database (a staged clone of a canonical one)."),
    staging: Path = typer.Option(..., "--staging", "-s"),
    run: list[str] = typer.Option([], "--run", "-r", help="Clustering run to load (repeatable; default: all)."),
    verify_digest: bool = typer.Option(True, "--verify-digest/--no-verify-digest"),
):
    """Writer: copy the clustering runs of a staging file into DB (same-name runs are replaced)."""
    from sharur.protein_families import store  # noqa: PLC0415

    con = store.connect(db)
    try:
        result = store.load_staging(con, staging, runs=run or None, verify_digest=verify_digest,
                                    log=lambda m: typer.echo(m, err=True))
    except store.FamilyError as exc:
        typer.echo(f"load: {exc}", err=True)
        raise typer.Exit(1) from exc
    finally:
        con.close()
    _echo_json(result)


def _staged(db: Path, label: str, work_dir: Optional[Path], writer, do_publish: bool) -> None:
    """Run ``writer(path)`` on a staged clone of ``db``, diff the family tables, optionally publish."""
    from sharur.protein_families import store  # noqa: PLC0415
    from sharur.staging import StagingError, diff_databases, discard, publish, stage  # noqa: PLC0415

    try:
        st = stage(db, work_dir=work_dir, label=label)
    except StagingError as exc:
        typer.echo(f"stage: {exc}", err=True)
        raise typer.Exit(1) from exc
    typer.echo(f"Staged {st.source} -> {st.staged} ({st.method})", err=True)
    try:
        result = writer(Path(st.staged))
    except Exception as exc:
        discard(st)
        typer.echo(f"{label}: {exc}; staged copy discarded, {st.source} unchanged", err=True)
        raise typer.Exit(1) from exc
    diff = diff_databases(st.source, st.staged, tables=list(store.TABLES), examples=3)
    typer.echo(diff.to_markdown())
    _echo_json(result)
    if not do_publish:
        typer.echo(f"Publish with `sharur stage publish {st.staged}`, or discard with "
                   f"`sharur stage discard {st.staged}`.")
        return
    try:
        published = publish(st)
    except StagingError as exc:
        typer.echo(f"publish: {exc}", err=True)
        raise typer.Exit(1) from exc
    typer.echo(f"Published {published.source}; backup {published.backup}; resealed {published.resealed}"
               + (f"; seal error {published.seal_error}" if published.seal_error else ""))
    if published.seal_error:
        raise typer.Exit(1)


@families_app.command("activate")
def activate_command(
    db: Path = typer.Option(..., "--db", "-d", help=_DB_HELP),
    staging: Path = typer.Option(..., "--staging", "-s", help="Staging DuckDB from `sharur families run/build`."),
    label: str = typer.Option("protein-families", "--label", "-l"),
    work_dir: Optional[Path] = typer.Option(None, "--work-dir",
                                            help="Directory for the staged clone (default: beside DB)."),
    run: list[str] = typer.Option([], "--run", "-r", help="Clustering run to load (repeatable; default: all)."),
    do_publish: bool = typer.Option(False, "--publish", help="Publish the staged copy when the load succeeds."),
    verify_digest: bool = typer.Option(True, "--verify-digest/--no-verify-digest"),
):
    """Stage DB, load the staged clustering runs into the clone, diff, and (with --publish) swap it in.

    Uses the staged-write protocol of `sharur stage` (docs/guides/staged-writes.md): the
    canonical file changes only at publish, a backup is kept and the dataset is resealed.
    """
    from sharur.protein_families import store  # noqa: PLC0415

    staging = staging.expanduser().resolve()

    def writer(path: Path) -> dict:
        con = store.connect(path)
        try:
            return store.load_staging(con, staging, runs=run or None, verify_digest=verify_digest,
                                      log=lambda m: typer.echo(m, err=True))
        finally:
            con.close()

    _staged(db.expanduser().resolve(), label, work_dir, writer, do_publish)


@families_app.command("assign")
def assign_command(
    db: Path = typer.Option(..., "--db", "-d", help=_DB_HELP),
    work_dir: Optional[Path] = typer.Option(None, "--work-dir", help="MMseqs2 scratch (default: beside DB)."),
    threads: Optional[int] = typer.Option(None, "--threads", "-t"),
    family_run: Optional[str] = typer.Option(None, "--family-run",
                                             help="Family run of the Sharur build to extend (default: the default run)."),
    in_place: bool = typer.Option(False, "--in-place",
                                  help="Write DB directly (a staged clone or scratch copy); default stages it."),
    label: str = typer.Option("protein-families-assign", "--label", "-l"),
    do_publish: bool = typer.Option(False, "--publish"),
    mmseqs_binary: Optional[str] = typer.Option(None, "--mmseqs"),
):
    """Give families to proteins that lack one (newly added genomes) without re-clustering."""
    from sharur.protein_families import store  # noqa: PLC0415

    db = db.expanduser().resolve()
    scratch = (work_dir or db.parent / f".{db.name}.families-assign").expanduser().resolve()

    def writer(path: Path) -> dict:
        con = store.connect(path)
        try:
            result = store.assign(con, scratch, family_run=family_run, threads=threads, binary=mmseqs_binary,
                                  log=lambda m: typer.echo(m, err=True)).to_dict()
            result["status"] = store.status(con)
            return result
        finally:
            con.close()

    if in_place:
        try:
            _echo_json(writer(db))
        except store.FamilyError as exc:
            typer.echo(f"assign: {exc}", err=True)
            raise typer.Exit(1) from exc
        return
    _staged(db, label, None, writer, do_publish)


@families_app.command("import")
def import_command(
    db: Path = typer.Option(..., "--db", "-d", help=_DB_HELP + " Opened read-only."),
    tsv: Path = typer.Option(..., "--tsv", help="Two-column representative<TAB>member table (mmseqs createtsv)."),
    staging: Path = typer.Option(..., "--staging", "-s",
                                 help="Staging DuckDB to add the run to (created when absent)."),
    name: str = typer.Option(..., "--name", "-n", help="Clustering run identifier, e.g. a tool and parameter tag."),
    level: str = typer.Option("imported", "--level", help="Level label (family, subfamily or any other label)."),
    parameters: str = typer.Option(..., "--parameters", "-p",
                                   help="The clustering command or parameters, stored verbatim."),
    tool: str = typer.Option("mmseqs2", "--tool"),
    tool_version: Optional[str] = typer.Option(None, "--tool-version"),
    set_default: bool = typer.Option(False, "--default",
                                     help="Make this run the default of its level (per-genome views)."),
    duckdb_threads: int = typer.Option(4, "--duckdb-threads"),
    memory_limit: str = typer.Option("8GB", "--memory-limit"),
):
    """Register an externally produced clustering as a run, without re-clustering."""
    from sharur.protein_families import store  # noqa: PLC0415

    try:
        result = store.import_clustering(db, tsv, staging, name=name, level=level, parameters=parameters,
                                         tool=tool, tool_version=tool_version, set_default=set_default,
                                         threads=duckdb_threads, memory_limit=memory_limit,
                                         log=lambda m: typer.echo(m, err=True))
    except store.FamilyError as exc:
        typer.echo(f"import: {exc}", err=True)
        raise typer.Exit(1) from exc
    _echo_json(result)


@families_app.command("status")
def status_command(db: Path = typer.Option(..., "--db", "-d", help=_DB_HELP + " Opened read-only.")):
    """Runs, family counts, proteins without a family, and whether a full re-cluster is due."""
    from sharur.protein_families import store  # noqa: PLC0415

    con = store.connect(db, read_only=True)
    try:
        _echo_json(store.status(con))
    finally:
        con.close()


__all__ = ["families_app"]
