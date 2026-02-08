#!/usr/bin/env python3
"""
Stage 4: Functional Annotation with Astra/PyHMMer
Scan protein sequences against domain databases for functional annotation.
"""

import logging
import json
import subprocess
import time
import shutil
from pathlib import Path
from typing import List, Optional, Dict, Any
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd

import typer
from rich.console import Console
from rich.table import Table
from rich.progress import Progress

console = Console()
logger = logging.getLogger(__name__)


def run_single_astra_scan(database: str, protein_symlink_dir: Path, output_dir: Path, 
                          threads: int, use_cutoffs: bool = True) -> Dict[str, Any]:
    """
    Run astra search for a single database.
    
    Args:
        database: Database name (PFAM, KOFAM, etc.)
        protein_symlink_dir: Directory containing protein file symlinks
        output_dir: Output directory for results
        threads: Number of threads to use
        use_cutoffs: Whether to use gathering cutoffs
        
    Returns:
        Dict containing execution results and statistics
    """
    start_time = time.time()
    
    # Create database-specific output directory
    db_output_dir = output_dir / f"{database.lower()}_results"
    db_output_dir.mkdir(parents=True, exist_ok=True)
    
    result = {
        "database": database,
        "execution_status": "failed",
        "execution_time_seconds": 0.0,
        "error_message": None,
        "output_dir": str(db_output_dir),
        "hits_file": None,
        "total_hits": 0,
        "unique_proteins": 0,
        "unique_domains": 0
    }
    
    try:
        # Build astra search command
        cmd = [
            "astra", "search",
            "--prot_in", str(protein_symlink_dir),
            "--installed_hmms", database,
            "--outdir", str(db_output_dir),
            "--threads", str(threads)
        ]
        
        # Add cutoffs for databases that support them
        # PFAM, KOFAM, and HydDB all have GA (gathering) thresholds
        if use_cutoffs and database.upper() in ["PFAM", "KOFAM", "HYDDB"]:
            cmd.append("--cut_ga")
        
        console.print(f"Running astra search for {database}...")
        console.print(f"Command: {' '.join(cmd)}")
        
        # Execute astra search
        process_result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            # No timeout — KOFAM on large datasets can take hours
        )
        
        if process_result.returncode != 0:
            result["error_message"] = f"Astra search failed: {process_result.stderr}"
            return result
        
        # Find the results file
        hits_file = db_output_dir / f"{database}_hits_df.tsv"
        if not hits_file.exists():
            result["error_message"] = f"Expected results file not found: {hits_file}"
            return result
        
        # Parse results statistics
        try:
            df = pd.read_csv(hits_file, sep='\t')
            result["total_hits"] = len(df)
            result["unique_proteins"] = df['sequence_id'].nunique() if 'sequence_id' in df.columns else 0
            result["unique_domains"] = df['hmm_name'].nunique() if 'hmm_name' in df.columns else 0
            result["hits_file"] = str(hits_file)
        except Exception as e:
            logger.warning(f"Failed to parse results statistics for {database}: {e}")
        
        result["execution_status"] = "success"
        
    except Exception as e:
        result["error_message"] = f"Unexpected error: {str(e)}"
    
    result["execution_time_seconds"] = round(time.time() - start_time, 2)
    return result


def run_astra_scan(
    input_dir: Path = typer.Option(
        Path("data/stage03_prodigal"),
        "--input-dir", "-i",
        help="Directory containing Prodigal protein sequences"
    ),
    output_dir: Path = typer.Option(
        Path("data/stage04_astra"),
        "--output-dir", "-o",
        help="Output directory for Astra scan results"
    ),
    threads: int = typer.Option(
        8,
        "--threads", "-t",
        help="Number of threads to use"
    ),
    databases: List[str] = typer.Option(
        ["PFAM", "KOFAM", "HydDB", "DefenseFinder"],
        "--databases", "-d",
        help="Databases to scan against (PFAM, KOFAM, HydDB, DefenseFinder, VOGdb, etc.)"
    ),
    use_cutoffs: bool = typer.Option(
        True,
        "--cutoffs/--no-cutoffs",
        help="Use gathering cutoffs for databases that support them"
    ),
    force: bool = typer.Option(
        False,
        "--force",
        help="Overwrite existing output directory"
    )
) -> None:
    """
    Run Astra functional annotation scans using PyHMMer.
    
    Scans protein sequences from Prodigal against HMM domain databases
    using the astra tool. Supports PFAM, KOFAM, and other installed databases.
    """
    console.print("[bold blue]Stage 4: Astra Functional Annotation[/bold blue]")
    
    # Validate inputs
    if not input_dir.exists():
        console.print(f"[red]Error: Input directory does not exist: {input_dir}[/red]")
        raise typer.Exit(1)
    
    # Check for symlink directory from prodigal stage
    protein_symlink_dir = input_dir / "genomes" / "all_protein_symlinks"
    if not protein_symlink_dir.exists():
        console.print(f"[red]Error: Protein symlink directory not found: {protein_symlink_dir}[/red]")
        console.print("[yellow]Please run Stage 3 (prodigal) first[/yellow]")
        raise typer.Exit(1)
    
    # Check for prodigal manifest
    manifest_file = input_dir / "processing_manifest.json"
    if not manifest_file.exists():
        console.print(f"[red]Error: Input manifest not found: {manifest_file}[/red]")
        raise typer.Exit(1)
    
    # Load input manifest
    try:
        with open(manifest_file, 'r') as f:
            input_manifest = json.load(f)
    except Exception as e:
        console.print(f"[red]Error: Failed to load input manifest: {e}[/red]")
        raise typer.Exit(1)
    
    # Count protein files
    protein_files = list(protein_symlink_dir.glob("*.faa"))
    if not protein_files:
        console.print(f"[red]Error: No protein files found in {protein_symlink_dir}[/red]")
        raise typer.Exit(1)
    
    console.print(f"Found {len(protein_files)} protein files to process")
    console.print(f"Databases to scan: {', '.join(databases)}")
    console.print(f"Threads: {threads}")
    console.print(f"Use cutoffs: {use_cutoffs}")
    
    # Handle output directory
    if output_dir.exists():
        if not force:
            console.print(f"[red]Error: Output directory already exists: {output_dir}[/red]")
            console.print("[yellow]Use --force to overwrite[/yellow]")
            raise typer.Exit(1)
        else:
            console.print(f"[yellow]Removing existing output directory: {output_dir}[/yellow]")
            shutil.rmtree(output_dir)
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "logs").mkdir(exist_ok=True)
    
    # Process databases sequentially to avoid overwhelming the system
    console.print("\n[bold yellow]Running Astra scans...[/bold yellow]")
    
    start_time = time.time()
    results = []
    
    for database in databases:
        console.print(f"\n[cyan]Processing {database}...[/cyan]")
        
        try:
            result = run_single_astra_scan(
                database=database,
                protein_symlink_dir=protein_symlink_dir,
                output_dir=output_dir,
                threads=threads,
                use_cutoffs=use_cutoffs
            )
            results.append(result)
            
            # Show result
            status = "✓" if result["execution_status"] == "success" else "✗"
            console.print(f"{status} {database}: {result.get('total_hits', 0):,} hits in {result['execution_time_seconds']:.1f}s")
            
            if result["execution_status"] == "failed":
                console.print(f"[red]  Error: {result.get('error_message', 'Unknown error')}[/red]")
                
        except Exception as e:
            error_result = {
                "database": database,
                "execution_status": "failed",
                "error_message": f"Execution error: {str(e)}",
                "execution_time_seconds": 0.0
            }
            results.append(error_result)
            console.print(f"✗ {database}: Execution failed - {str(e)}")
    
    total_time = time.time() - start_time
    
    # Generate summary statistics
    successful_results = [r for r in results if r["execution_status"] == "success"]
    failed_results = [r for r in results if r["execution_status"] == "failed"]
    
    total_hits = sum(r.get("total_hits", 0) for r in successful_results)
    total_unique_proteins = sum(r.get("unique_proteins", 0) for r in successful_results)
    total_unique_domains = sum(r.get("unique_domains", 0) for r in successful_results)
    
    summary_stats = {
        "total_databases": len(databases),
        "successful_databases": len(successful_results),
        "failed_databases": len(failed_results),
        "total_hits": total_hits,
        "total_unique_proteins": total_unique_proteins,
        "total_unique_domains": total_unique_domains,
        "total_execution_time_seconds": round(total_time, 2),
        "protein_files_processed": len(protein_files)
    }
    
    # Create processing manifest
    manifest = {
        "version": "0.1.0",
        "stage": "stage04_astra",
        "timestamp": datetime.now().isoformat(),
        "input_manifest": str(manifest_file.absolute()),
        "execution_parameters": {
            "databases": list(databases),
            "threads": threads,
            "use_cutoffs": bool(use_cutoffs),
            "protein_symlink_dir": str(protein_symlink_dir)
        },
        "summary": summary_stats,
        "database_results": results
    }
    
    # Save manifest
    output_manifest = output_dir / "processing_manifest.json"
    with open(output_manifest, 'w') as f:
        json.dump(manifest, f, indent=2)
    
    # Save summary statistics
    summary_file = output_dir / "summary_stats.json"
    with open(summary_file, 'w') as f:
        json.dump(summary_stats, f, indent=2)
    
    # Display results
    console.print("\n[bold green]Stage 4 Results Summary[/bold green]")
    
    summary_table = Table()
    summary_table.add_column("Metric", style="cyan")
    summary_table.add_column("Value", style="magenta")
    
    summary_table.add_row("Input directory", str(input_dir))
    summary_table.add_row("Output directory", str(output_dir))
    summary_table.add_row("Protein files processed", str(summary_stats["protein_files_processed"]))
    summary_table.add_row("Total databases", str(summary_stats["total_databases"]))
    summary_table.add_row("Successful databases", str(summary_stats["successful_databases"]))
    summary_table.add_row("Failed databases", str(summary_stats["failed_databases"]))
    summary_table.add_row("Total hits", f"{summary_stats['total_hits']:,}")
    summary_table.add_row("Total unique proteins with hits", f"{summary_stats['total_unique_proteins']:,}")
    summary_table.add_row("Total unique domains", f"{summary_stats['total_unique_domains']:,}")
    summary_table.add_row("Execution time", f"{summary_stats['total_execution_time_seconds']:.1f} seconds")
    
    console.print(summary_table)
    
    # Show failed databases if any
    if failed_results:
        console.print("\n[bold red]Failed Databases:[/bold red]")
        for result in failed_results:
            console.print(f"[red]• {result['database']}: {result.get('error_message', 'Unknown error')}[/red]")
    
    # Success message
    if successful_results:
        console.print(f"\n[bold green]✓ Stage 4 completed successfully![/bold green]")
        console.print(f"Annotation results available in: {output_dir}/*/")
        
    logger.info(f"Stage 4 completed: {len(successful_results)} successful, {len(failed_results)} failed")


def main():
    """Entry point for standalone execution."""
    logging.basicConfig(level=logging.INFO)
    logger.info("Stage 4 stub")
    typer.run(run_astra_scan)


if __name__ == "__main__":
    main()
