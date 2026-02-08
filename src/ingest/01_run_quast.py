#!/usr/bin/env python3
"""
Stage 1: Quality Assessment with QUAST
Assess assembly quality metrics for genome assemblies.
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
import os
import multiprocessing

import typer
from rich.console import Console
from rich.table import Table
from rich.progress import Progress, TaskID

console = Console()
logger = logging.getLogger(__name__)


def parse_quast_report(report_file: Path) -> Dict[str, Any]:
    """
    Parse QUAST report.tsv file to extract quality metrics.
    
    Returns:
        Dict containing assembly quality statistics
    """
    stats = {
        "contigs": 0,
        "total_length": 0,
        "largest_contig": 0,
        "n50": 0,
        "n75": 0,
        "gc_content": 0.0,
        "n_count": 0,
        "n_per_100_kbp": 0.0,
        "contigs_1000bp": 0,
        "contigs_5000bp": 0,
        "contigs_10000bp": 0,
        "l50": 0,
        "l75": 0
    }
    
    try:
        if not report_file.exists():
            return stats
            
        with open(report_file, 'r') as f:
            lines = f.readlines()
            
        # QUAST report is tab-separated with metrics in first column, values in second
        for line in lines:
            if '\t' not in line:
                continue
                
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
                
            metric = parts[0].strip()
            value_str = parts[1].strip()
            
            try:
                # Parse different metric types
                if metric == "# contigs":
                    stats["contigs"] = int(value_str)
                elif metric == "Total length":
                    stats["total_length"] = int(value_str.replace(',', ''))
                elif metric == "Largest contig":
                    stats["largest_contig"] = int(value_str.replace(',', ''))
                elif metric == "N50":
                    stats["n50"] = int(value_str.replace(',', ''))
                elif metric == "N75":
                    stats["n75"] = int(value_str.replace(',', ''))
                elif metric == "GC (%)":
                    stats["gc_content"] = float(value_str) / 100.0
                elif metric == "# N's":
                    stats["n_count"] = int(value_str.replace(',', ''))
                elif metric == "# N's per 100 kbp":
                    stats["n_per_100_kbp"] = float(value_str)
                elif metric == "# contigs (>= 1000 bp)":
                    stats["contigs_1000bp"] = int(value_str)
                elif metric == "# contigs (>= 5000 bp)":
                    stats["contigs_5000bp"] = int(value_str)
                elif metric == "# contigs (>= 10000 bp)":
                    stats["contigs_10000bp"] = int(value_str)
                elif metric == "L50":
                    stats["l50"] = int(value_str)
                elif metric == "L75":
                    stats["l75"] = int(value_str)
            except (ValueError, IndexError):
                # Skip metrics that can't be parsed
                continue
                
    except Exception as e:
        logger.warning(f"Failed to parse QUAST report from {report_file}: {e}")
        
    return stats


def validate_quast_outputs(genome_output_dir: Path, genome_id: str) -> Dict[str, Any]:
    """
    Validate QUAST output files and extract basic information.
    
    Returns:
        Dict containing validation results and file information
    """
    validation = {
        "output_files_exist": True,
        "report_valid": True,
        "missing_files": [],
        "report_path": None,
        "has_plots": False
    }
    
    # Check for required QUAST output files and directories
    report_file = genome_output_dir / "report.tsv"
    basic_stats_dir = genome_output_dir / "basic_stats"
    
    # Check report file
    if not report_file.exists():
        validation["output_files_exist"] = False
        validation["missing_files"].append(str(report_file))
        validation["report_valid"] = False
    else:
        validation["report_path"] = str(report_file)
        # Check if report has content
        try:
            with open(report_file, 'r') as f:
                content = f.read()
                if len(content.strip()) == 0:
                    validation["report_valid"] = False
        except Exception:
            validation["report_valid"] = False
    
    # Check if basic stats directory exists
    if basic_stats_dir.exists() and any(basic_stats_dir.iterdir()):
        validation["has_plots"] = True
    
    return validation


def run_quast_single(genome_info: Dict[str, Any], 
                    output_base_dir: Path,
                    min_contig_length: int = 500,
                    threads: int = 1,
                    reference_genome: Optional[Path] = None) -> Dict[str, Any]:
    """
    Run QUAST on a single genome.
    
    Args:
        genome_info: Dictionary containing genome metadata from input manifest
        output_base_dir: Base output directory for all genomes
        min_contig_length: Minimum contig length for analysis
        threads: Number of threads to use per genome
        reference_genome: Optional reference genome for comparison
        
    Returns:
        Dict containing execution results and statistics
    """
    start_time = time.time()
    
    # Extract genome information
    genome_id = genome_info["genome_id"]
    input_file = Path(genome_info["output_path"])  # From stage00 manifest
    
    # Create output directory for this genome
    genome_output_dir = output_base_dir / "genomes" / genome_id
    genome_output_dir.mkdir(parents=True, exist_ok=True)
    
    # Define log file
    log_file = genome_output_dir / "quast.log"
    
    # Initialize result dictionary
    result = {
        "genome_id": genome_id,
        "input_file": str(input_file),
        "output_dir": str(genome_output_dir),
        "execution_status": "failed",
        "execution_time_seconds": 0.0,
        "error_message": None,
        "quast_version": None,
        "statistics": {},
        "validation": {},
        "output_files": {}
    }
    
    try:
        # Build QUAST command
        cmd = [
            "quast.py",
            str(input_file),
            "--output-dir", str(genome_output_dir),
            "--threads", str(threads),
            "--min-contig", str(min_contig_length),
            "--no-html",  # Skip HTML report generation for faster processing
            "--no-icarus",  # Skip Icarus visualization
            "--no-snps",   # Skip SNP analysis
            "--no-plots"   # Skip plot generation for faster processing
        ]
        
        # Add reference if provided
        if reference_genome and reference_genome.exists():
            cmd.extend(["--reference", str(reference_genome)])
            
        # Execute QUAST
        with open(log_file, 'w') as log_f:
            process_result = subprocess.run(
                cmd,
                stdout=log_f,
                stderr=subprocess.STDOUT,
                text=True,
                timeout=600  # 10 minute timeout per genome
            )
            
        if process_result.returncode != 0:
            result["error_message"] = f"QUAST failed with return code {process_result.returncode}"
            return result
            
        # Get QUAST version
        try:
            version_result = subprocess.run(
                ["quast.py", "--version"],
                capture_output=True,
                text=True,
                timeout=10
            )
            if version_result.returncode == 0:
                result["quast_version"] = version_result.stdout.strip()
        except Exception:
            pass  # Version detection is non-critical
            
        # Parse statistics from report
        report_file = genome_output_dir / "report.tsv"
        result["statistics"] = parse_quast_report(report_file)
        
        # Validate outputs
        result["validation"] = validate_quast_outputs(genome_output_dir, genome_id)
        
        # Record output files
        result["output_files"] = {
            "report": str((genome_output_dir / "report.tsv").relative_to(output_base_dir)),
            "log": str(log_file.relative_to(output_base_dir))
        }
        
        # Add additional files if they exist
        if (genome_output_dir / "contigs_reports").exists():
            result["output_files"]["contigs_report"] = str((genome_output_dir / "contigs_reports").relative_to(output_base_dir))
        if (genome_output_dir / "basic_stats").exists():
            result["output_files"]["basic_stats"] = str((genome_output_dir / "basic_stats").relative_to(output_base_dir))
            
        # Check if execution was successful
        if (result["validation"]["output_files_exist"] and 
            result["validation"]["report_valid"] and
            result["statistics"]["contigs"] > 0):
            result["execution_status"] = "success"
        else:
            result["execution_status"] = "failed"
            result["error_message"] = "Output validation failed"
            
    except subprocess.TimeoutExpired:
        result["error_message"] = "QUAST execution timed out (>10 minutes)"
    except Exception as e:
        result["error_message"] = f"Unexpected error: {str(e)}"
        
    result["execution_time_seconds"] = round(time.time() - start_time, 2)
    return result


def process_genomes_parallel(genomes: List[Dict[str, Any]],
                           output_dir: Path,
                           max_workers: int,
                           **quast_kwargs) -> List[Dict[str, Any]]:
    """
    Process multiple genomes in parallel using ProcessPoolExecutor.
    
    Returns:
        List of dictionaries containing results for each genome
    """
    results = []
    
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit all jobs
        future_to_genome = {}
        for genome_info in genomes:
            future = executor.submit(
                run_quast_single,
                genome_info,
                output_dir,
                **quast_kwargs
            )
            future_to_genome[future] = genome_info["genome_id"]
        
        # Collect results as they complete
        with Progress() as progress:
            task = progress.add_task("Processing genomes...", total=len(genomes))
            
            for future in as_completed(future_to_genome):
                genome_id = future_to_genome[future]
                try:
                    result = future.result()
                    results.append(result)
                    
                    # Update progress
                    status = "✓" if result["execution_status"] == "success" else "✗"
                    progress.update(task, 
                                  description=f"Processing genomes... {status} {genome_id}",
                                  advance=1)
                    
                except Exception as e:
                    # Handle unexpected errors
                    error_result = {
                        "genome_id": genome_id,
                        "execution_status": "failed",
                        "error_message": f"Execution error: {str(e)}",
                        "execution_time_seconds": 0.0
                    }
                    results.append(error_result)
                    progress.update(task, 
                                  description=f"Processing genomes... ✗ {genome_id}",
                                  advance=1)
    
    return results


def generate_summary_stats(results: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Generate aggregate statistics across all genomes."""
    successful = [r for r in results if r["execution_status"] == "success"]
    failed = [r for r in results if r["execution_status"] == "failed"]
    
    # Calculate aggregate statistics
    total_contigs = sum(r.get("statistics", {}).get("contigs", 0) for r in successful)
    total_length = sum(r.get("statistics", {}).get("total_length", 0) for r in successful)
    
    mean_execution_time = 0.0
    if successful:
        mean_execution_time = sum(r.get("execution_time_seconds", 0) for r in successful) / len(successful)
        
    # Quality metrics aggregation
    n50_values = [r.get("statistics", {}).get("n50", 0) for r in successful if r.get("statistics", {}).get("n50", 0) > 0]
    gc_values = [r.get("statistics", {}).get("gc_content", 0) for r in successful if r.get("statistics", {}).get("gc_content", 0) > 0]
    
    mean_n50 = sum(n50_values) / len(n50_values) if n50_values else 0
    mean_gc = sum(gc_values) / len(gc_values) if gc_values else 0
    
    summary = {
        "total_genomes": len(results),
        "successful": len(successful),
        "failed": len(failed),
        "success_rate": len(successful) / len(results) if results else 0.0,
        "total_contigs": total_contigs,
        "total_assembly_length": total_length,
        "mean_n50": round(mean_n50, 0),
        "mean_gc_content": round(mean_gc, 3),
        "mean_execution_time_seconds": round(mean_execution_time, 2),
        "total_execution_time_seconds": round(sum(r.get("execution_time_seconds", 0) for r in results), 2)
    }
    
    return summary


def run_quast(
    input_dir: Path = typer.Option(
        Path("data/stage00_prepared"),
        "--input-dir", "-i",
        help="Directory containing input genome assemblies"
    ),
    output_dir: Path = typer.Option(
        Path("data/stage01_quast"),
        "--output-dir", "-o",
        help="Output directory for QUAST results"
    ),
    min_contig_length: int = typer.Option(
        500,
        "--min-contig", "-m",
        help="Minimum contig length for analysis"
    ),
    threads_per_genome: int = typer.Option(
        1,
        "--threads-per-genome", "-t",
        help="Number of threads per genome (total threads = max_workers * threads_per_genome)"
    ),
    max_workers: Optional[int] = typer.Option(
        None,
        "--max-workers", "-w",
        help="Maximum number of parallel workers (default: CPU count - 1)"
    ),
    reference_genome: Optional[Path] = typer.Option(
        None,
        "--reference", "-r",
        help="Reference genome for comparison (optional)"
    ),
    force: bool = typer.Option(
        False,
        "--force",
        help="Overwrite existing output directory"
    )
) -> None:
    """
    Run QUAST quality assessment on genome assemblies.
    
    Analyzes assembly quality metrics including N50, contig counts, total length,
    GC content, and other statistics. Processes genomes in parallel for efficiency.
    """
    console.print("[bold blue]Stage 1: QUAST Quality Assessment[/bold blue]")
    
    # Validate inputs
    if not input_dir.exists():
        console.print(f"[red]Error: Input directory does not exist: {input_dir}[/red]")
        raise typer.Exit(1)
        
    manifest_file = input_dir / "processing_manifest.json"
    if not manifest_file.exists():
        console.print(f"[red]Error: Input manifest not found: {manifest_file}[/red]")
        console.print("[yellow]Please run Stage 0 (prepare_inputs) first[/yellow]")
        raise typer.Exit(1)
    
    # Load input manifest
    try:
        with open(manifest_file, 'r') as f:
            input_manifest = json.load(f)
    except Exception as e:
        console.print(f"[red]Error: Failed to load input manifest: {e}[/red]")
        raise typer.Exit(1)
    
    # Filter for valid genomes only
    valid_genomes = [g for g in input_manifest["genomes"] if g.get("format_valid", True)]
    if not valid_genomes:
        console.print("[red]Error: No valid genomes found in input manifest[/red]")
        raise typer.Exit(1)
    
    console.print(f"Found {len(valid_genomes)} valid genomes to process")
    
    # Validate reference genome if provided
    if reference_genome and not reference_genome.exists():
        console.print(f"[red]Error: Reference genome not found: {reference_genome}[/red]")
        raise typer.Exit(1)
    
    # Handle output directory
    if output_dir.exists():
        if not force:
            console.print(f"[red]Error: Output directory already exists: {output_dir}[/red]")
            console.print("[yellow]Use --force to overwrite[/yellow]")
            raise typer.Exit(1)
        else:
            console.print(f"[yellow]Removing existing output directory: {output_dir}[/yellow]")
            shutil.rmtree(output_dir)
    
    # Create output directories
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "genomes").mkdir(exist_ok=True)
    (output_dir / "logs").mkdir(exist_ok=True)
    
    # Determine number of workers
    if max_workers is None:
        max_workers = max(1, multiprocessing.cpu_count() - 1)
    
    console.print(f"Using {max_workers} parallel workers")
    console.print(f"Threads per genome: {threads_per_genome}")
    console.print(f"Min contig length: {min_contig_length}")
    if reference_genome:
        console.print(f"Reference genome: {reference_genome}")
    
    # Process genomes in parallel
    console.print("\n[bold yellow]Processing genomes...[/bold yellow]")
    
    start_time = time.time()
    results = process_genomes_parallel(
        valid_genomes,
        output_dir,
        max_workers,
        min_contig_length=min_contig_length,
        threads=threads_per_genome,
        reference_genome=reference_genome
    )
    total_time = time.time() - start_time
    
    # Generate summary statistics
    summary_stats = generate_summary_stats(results)
    
    # Create processing manifest
    manifest = {
        "version": "0.1.0",
        "stage": "stage01_quast",
        "timestamp": datetime.now().isoformat(),
        "input_manifest": str(manifest_file.absolute()),
        "execution_parameters": {
            "min_contig_length": min_contig_length,
            "threads_per_genome": threads_per_genome,
            "max_workers": max_workers,
            "reference_genome": str(reference_genome) if reference_genome else None
        },
        "summary": summary_stats,
        "genomes": results
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
    console.print("\n[bold green]Stage 1 Results Summary[/bold green]")
    
    summary_table = Table()
    summary_table.add_column("Metric", style="cyan")
    summary_table.add_column("Value", style="magenta")
    
    summary_table.add_row("Input directory", str(input_dir))
    summary_table.add_row("Output directory", str(output_dir))
    summary_table.add_row("Total genomes", str(summary_stats["total_genomes"]))
    summary_table.add_row("Successful", str(summary_stats["successful"]))
    summary_table.add_row("Failed", str(summary_stats["failed"]))
    summary_table.add_row("Success rate", f"{summary_stats['success_rate']:.1%}")
    summary_table.add_row("Total contigs", f"{summary_stats['total_contigs']:,}")  
    summary_table.add_row("Total assembly length", f"{summary_stats['total_assembly_length']:,} bp")
    summary_table.add_row("Mean N50", f"{summary_stats['mean_n50']:,.0f} bp")
    summary_table.add_row("Mean GC content", f"{summary_stats['mean_gc_content']:.1%}")
    summary_table.add_row("Execution time", f"{total_time:.1f} seconds")
    summary_table.add_row("Mean time per genome", f"{summary_stats['mean_execution_time_seconds']:.1f} seconds")
    
    console.print(summary_table)
    
    # Show failed genomes if any
    failed_results = [r for r in results if r["execution_status"] == "failed"]
    if failed_results:
        console.print("\n[bold red]Failed Genomes:[/bold red]")
        for result in failed_results:
            console.print(f"[red]• {result['genome_id']}: {result.get('error_message', 'Unknown error')}[/red]")
    
    # Success message
    if summary_stats["successful"] > 0:
        console.print(f"\n[bold green]✓ Stage 1 completed successfully![/bold green]")
        console.print(f"QUAST reports available in: {output_dir}/genomes/*/")
        
    logger.info(f"Stage 1 completed: {summary_stats['successful']} successful, {summary_stats['failed']} failed")


def main():
    """Entry point for standalone execution."""
    logging.basicConfig(level=logging.INFO)
    typer.run(run_quast)


if __name__ == "__main__":
    main()
