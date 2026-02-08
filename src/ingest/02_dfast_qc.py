#!/usr/bin/env python3
"""
Stage 2: Taxonomic Classification with DFAST_QC
Evaluate genome completeness and assign taxonomic classifications using DFAST_QC.
"""

import json
import logging
import subprocess
import shutil
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import typer
from rich.console import Console
from rich.progress import Progress
from rich.table import Table

console = Console()
logger = logging.getLogger(__name__)


def parse_dqc_json(dqc_result_path: Path) -> Dict[str, Any]:
    """
    Parse DFAST_QC JSON result and extract key taxonomy/quality information.
    
    Args:
        dqc_result_path: Path to dqc_result.json file
        
    Returns:
        Dict containing parsed taxonomy and quality data
    """
    try:
        with open(dqc_result_path, 'r') as f:
            dqc_data = json.load(f)
        
        # Extract taxonomy information
        taxonomy = dqc_data.get("taxonomy", {})
        
        # Extract quality information
        quality = dqc_data.get("quality", {})
        checkm = quality.get("checkm", {})
        
        # Build summary dict with standardized keys
        summary = {
            "rank": taxonomy.get("rank", "unknown"),
            "name": taxonomy.get("species", taxonomy.get("genus", "unknown")),
            "ani": taxonomy.get("ani", 0.0),
            "status": "complete" if checkm.get("completeness", 0) > 90 else "partial",
            "completeness": checkm.get("completeness", 0.0),
            "contamination": checkm.get("contamination", 0.0),
            "tool": "dfast_qc",
            "version": dqc_data.get("version", "unknown"),
            "confidence": taxonomy.get("confidence", 0.0)
        }
        
        return summary
        
    except Exception as e:
        logger.warning(f"Failed to parse DFAST_QC result from {dqc_result_path}: {e}")
        return {
            "rank": "unknown",
            "name": "unknown", 
            "ani": 0.0,
            "status": "failed",
            "completeness": 0.0,
            "contamination": 0.0,
            "tool": "dfast_qc",
            "version": "unknown",
            "confidence": 0.0
        }


def run_dfast(
    fasta: Path,
    out: Path, 
    threads: int = 8,
    enable_cc: bool = True
) -> Dict[str, Any]:
    """
    Run DFAST_QC on a single genome.
    
    Args:
        fasta: Input FASTA file path
        out: Output directory path
        threads: Number of threads to use
        enable_cc: Enable completeness/contamination analysis
        
    Returns:
        Dict containing execution results and parsed taxonomy data
    """
    start_time = time.time()
    
    # Create output directory
    out.mkdir(parents=True, exist_ok=True)
    
    # Initialize result dictionary
    result = {
        "genome_id": fasta.stem,
        "input_file": str(fasta),
        "output_dir": str(out),
        "execution_status": "failed",
        "execution_time_seconds": 0.0,
        "error_message": None,
        "dfast_qc_version": None,
        "taxonomy": {},
        "output_files": {}
    }
    
    # Light guardrail: skip very large inputs to avoid runaway runtimes
    max_bytes = 25 * 1024 * 1024  # 25 MB
    if fasta.stat().st_size > max_bytes:
        result["error_message"] = f"Input file {fasta.name} exceeds 25 MB; skipping."
        logger.warning(result["error_message"])
        return result

    try:
        exe = shutil.which("dfast_qc")
        if exe is None:
            raise FileNotFoundError("dfast_qc executable not found on PATH")
        # Build DFAST_QC command
        cmd = [
            exe,
            "-i", str(fasta),
            "-o", str(out),
            "--num_threads", str(threads)
        ]
        
        # Add completeness/contamination flag
        if not enable_cc:
            cmd.append("--disable_cc")
        
        # Create log file
        log_file = out / "dfast_qc.log"
        
        # Execute DFAST_QC with tempfile for scratch space
        with tempfile.TemporaryDirectory() as temp_dir:
            # Run DFAST_QC
            with open(log_file, 'w') as log_f:
                process_result = subprocess.run(
                    cmd,
                    stdout=log_f,
                    stderr=subprocess.STDOUT,
                    text=True,
                    timeout=1800,  # 30 minute timeout per genome
                    env={"TMPDIR": temp_dir}
                )
            
            if process_result.returncode != 0:
                # Pull tail of log for easier debugging
                tail = ""
                try:
                    tail = "".join(log_file.read_text().splitlines()[-10:])
                except Exception:
                    pass
                result["error_message"] = f"DFAST_QC failed with return code {process_result.returncode}. Log tail: {tail}"
                return result
        
        # Get DFAST_QC version
        try:
            version_result = subprocess.run(
                [exe, "--version"],
                capture_output=True,
                text=True,
                timeout=10
            )
            if version_result.returncode == 0:
                result["dfast_qc_version"] = version_result.stdout.strip()
        except Exception:
            pass  # Version detection is non-critical
        
        # Check for output files
        dqc_result_file = out / "dqc_result.json"
        if not dqc_result_file.exists():
            result["error_message"] = "DFAST_QC result file not found"
            return result
        
        # Parse DFAST_QC results
        taxonomy_summary = parse_dqc_json(dqc_result_file)
        result["taxonomy"] = taxonomy_summary
        
        # Write trimmed summary
        tax_summary_file = out / "tax_summary.json"
        with open(tax_summary_file, 'w') as f:
            json.dump(taxonomy_summary, f, indent=2)
        
        # Record output files
        result["output_files"] = {
            "dqc_result": str(dqc_result_file.relative_to(out.parent)),
            "tax_summary": str(tax_summary_file.relative_to(out.parent)),
            "log": str(log_file.relative_to(out.parent))
        }
        
        # Check if execution was successful
        if (dqc_result_file.exists() and 
            tax_summary_file.exists() and
            taxonomy_summary.get("rank") != "unknown"):
            result["execution_status"] = "success"
        else:
            result["execution_status"] = "failed"
            result["error_message"] = "Output validation failed"
            
    except subprocess.TimeoutExpired:
        result["error_message"] = "DFAST_QC execution timed out (>30 minutes)"
    except Exception as e:
        result["error_message"] = f"Unexpected error: {str(e)}"
        
    result["execution_time_seconds"] = round(time.time() - start_time, 2)
    return result


def run_dfast_single(
    genome_info: Dict[str, Any],
    output_base_dir: Path,
    threads: int = 8,
    enable_cc: bool = True
) -> Dict[str, Any]:
    """
    Run DFAST_QC on a single genome (wrapper for parallel processing).
    
    Args:
        genome_info: Dictionary containing genome metadata
        output_base_dir: Base output directory
        threads: Number of threads per genome
        enable_cc: Enable completeness/contamination analysis
        
    Returns:
        Dict containing execution results
    """
    genome_id = genome_info["genome_id"]
    input_file = Path(genome_info["output_path"])
    
    # Create genome-specific output directory
    genome_output_dir = output_base_dir / "genomes" / genome_id
    
    return run_dfast(
        fasta=input_file,
        out=genome_output_dir,
        threads=threads,
        enable_cc=enable_cc
    )


def process_genomes_parallel(
    genomes: List[Dict[str, Any]],
    output_dir: Path,
    max_workers: int,
    **dfast_kwargs
) -> List[Dict[str, Any]]:
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
                run_dfast_single,
                genome_info,
                output_dir,
                **dfast_kwargs
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
    
    # Calculate completion and contamination stats
    completeness_values = [
        r.get("taxonomy", {}).get("completeness", 0) 
        for r in successful
    ]
    contamination_values = [
        r.get("taxonomy", {}).get("contamination", 0) 
        for r in successful
    ]
    
    # Calculate means
    mean_completeness = sum(completeness_values) / len(completeness_values) if completeness_values else 0
    mean_contamination = sum(contamination_values) / len(contamination_values) if contamination_values else 0
    mean_execution_time = sum(r.get("execution_time_seconds", 0) for r in successful) / len(successful) if successful else 0
    
    # Count by taxonomic rank
    rank_counts = {}
    for result in successful:
        rank = result.get("taxonomy", {}).get("rank", "unknown")
        rank_counts[rank] = rank_counts.get(rank, 0) + 1
    
    summary = {
        "total_genomes": len(results),
        "successful": len(successful),
        "failed": len(failed),
        "success_rate": len(successful) / len(results) if results else 0.0,
        "mean_completeness": round(mean_completeness, 1),
        "mean_contamination": round(mean_contamination, 1),
        "mean_execution_time_seconds": round(mean_execution_time, 2),
        "total_execution_time_seconds": round(sum(r.get("execution_time_seconds", 0) for r in results), 2),
        "taxonomic_ranks": rank_counts
    }
    
    return summary


def call(
    input_dir: Path = typer.Option(
        Path("data/stage00_prepared"),
        "--input-dir", "-i",
        help="Directory containing input genome assemblies"
    ),
    output_dir: Path = typer.Option(
        Path("data/stage02_dfast_qc"),
        "--output-dir", "-o",
        help="Output directory for DFAST_QC results"
    ),
    threads: int = typer.Option(
        8,
        "--threads", "-t",
        help="Number of threads per genome"
    ),
    max_workers: Optional[int] = typer.Option(
        None,
        "--max-workers", "-w",
        help="Maximum number of parallel workers (default: CPU count - 1)"
    ),
    enable_cc: bool = typer.Option(
        False,
        "--enable-cc/--disable-cc",
        help="Enable/disable completeness and contamination analysis"
    ),
    force: bool = typer.Option(
        False,
        "--force",
        help="Overwrite existing output directory"
    )
) -> None:
    """
    Run DFAST_QC taxonomic classification on genome assemblies.
    
    Performs taxonomic classification using Average Nucleotide Identity (ANI)
    and evaluates genome completeness/contamination using CheckM-style analysis.
    """
    console.print("[bold blue]Stage 2: DFAST_QC Taxonomic Classification[/bold blue]")
    
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
    
    # Handle output directory
    if output_dir.exists():
        if not force:
            console.print(f"[red]Error: Output directory already exists: {output_dir}[/red]")
            console.print("[yellow]Use --force to overwrite[/yellow]")
            raise typer.Exit(1)
        else:
            console.print(f"[yellow]Removing existing output directory: {output_dir}[/yellow]")
            import shutil
            shutil.rmtree(output_dir)
    
    # Create output directories
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "genomes").mkdir(exist_ok=True)
    
    # Determine number of workers
    if max_workers is None:
        import multiprocessing
        max_workers = max(1, multiprocessing.cpu_count() - 1)
    
    console.print(f"Using {max_workers} parallel workers")
    console.print(f"Threads per genome: {threads}")
    console.print(f"Completeness/contamination analysis: {'enabled' if enable_cc else 'disabled'}")
    
    # Process genomes in parallel
    console.print("\n[bold yellow]Processing genomes...[/bold yellow]")
    
    start_time = time.time()
    results = process_genomes_parallel(
        valid_genomes,
        output_dir,
        max_workers,
        threads=threads,
        enable_cc=enable_cc
    )
    total_time = time.time() - start_time
    
    # Generate summary statistics
    summary_stats = generate_summary_stats(results)
    
    # Create processing manifest
    manifest = {
        "version": "0.1.0",
        "stage": "stage02_dfast_qc",
        "timestamp": datetime.now().isoformat(),
        "input_manifest": str(manifest_file.absolute()),
        "execution_parameters": {
            "threads": threads,
            "max_workers": max_workers,
            "enable_cc": enable_cc
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
    console.print("\n[bold green]Stage 2 Results Summary[/bold green]")
    
    summary_table = Table()
    summary_table.add_column("Metric", style="cyan")
    summary_table.add_column("Value", style="magenta")
    
    summary_table.add_row("Input directory", str(input_dir))
    summary_table.add_row("Output directory", str(output_dir))
    summary_table.add_row("Total genomes", str(summary_stats["total_genomes"]))
    summary_table.add_row("Successful", str(summary_stats["successful"]))
    summary_table.add_row("Failed", str(summary_stats["failed"]))
    summary_table.add_row("Success rate", f"{summary_stats['success_rate']:.1%}")
    summary_table.add_row("Mean completeness", f"{summary_stats['mean_completeness']:.1f}%")
    summary_table.add_row("Mean contamination", f"{summary_stats['mean_contamination']:.1f}%")
    summary_table.add_row("Execution time", f"{total_time:.1f} seconds")
    summary_table.add_row("Mean time per genome", f"{summary_stats['mean_execution_time_seconds']:.1f} seconds")
    
    console.print(summary_table)
    
    # Show taxonomic rank distribution
    if summary_stats["taxonomic_ranks"]:
        console.print("\n[bold yellow]Taxonomic Distribution:[/bold yellow]")
        for rank, count in summary_stats["taxonomic_ranks"].items():
            console.print(f"[cyan]• {rank}:[/cyan] {count} genomes")
    
    # Show failed genomes if any
    failed_results = [r for r in results if r["execution_status"] == "failed"]
    if failed_results:
        console.print("\n[bold red]Failed Genomes:[/bold red]")
        for result in failed_results:
            console.print(f"[red]• {result['genome_id']}: {result.get('error_message', 'Unknown error')}[/red]")
    
    # Success message  
    if summary_stats["successful"] > 0:
        console.print(f"\n[bold green]✓ Stage 2 completed successfully![/bold green]")
        console.print(f"Taxonomy summaries available in: {output_dir}/genomes/*/tax_summary.json")
        
    logger.info(f"Stage 2 completed: {summary_stats['successful']} successful, {summary_stats['failed']} failed")


def main():
    """Entry point for standalone execution."""
    logging.basicConfig(level=logging.INFO)
    typer.run(call)


if __name__ == "__main__":
    main()
