#!/usr/bin/env python3
"""
Stage 3: Gene Prediction with Prodigal
Predict protein-coding sequences and extract genomic features.
"""

import logging
import json
import subprocess
import time
from pathlib import Path
from typing import List, Optional, Dict, Any
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
import multiprocessing
import shutil

import typer
from rich.console import Console
from rich.table import Table
from rich.progress import Progress, TaskID

console = Console()
logger = logging.getLogger(__name__)


def parse_prodigal_stats(log_file: Path) -> Dict[str, Any]:
    """
    Parse statistics from prodigal log file.
    
    Returns:
        Dict containing gene prediction statistics
    """
    stats = {
        "genes_predicted": 0,
        "mean_gene_length": 0,
        "coding_density": 0.0,
        "gc_content": 0.0,
        "complete_genes": 0,
        "partial_genes": 0
    }
    
    try:
        if not log_file.exists():
            return stats
            
        with open(log_file, 'r') as f:
            content = f.read()
            
        # Parse gene counts and statistics from prodigal output
        lines = content.split('\n')
        for line in lines:
            line = line.strip()
            if 'genes predicted' in line.lower():
                try:
                    # Extract number from line like "5432 genes predicted"
                    parts = line.split()
                    for i, part in enumerate(parts):
                        if part.isdigit():
                            stats["genes_predicted"] = int(part)
                            break
                except (ValueError, IndexError):
                    pass
            elif 'avg gene length' in line.lower():
                try:
                    # Extract from line like "Avg gene length: 924 bp"
                    parts = line.split(':')
                    if len(parts) > 1:
                        length_str = parts[1].split()[0]
                        stats["mean_gene_length"] = int(float(length_str))
                except (ValueError, IndexError):
                    pass
            elif 'coding density' in line.lower():
                try:
                    # Extract from line like "Coding density: 91.2%"
                    parts = line.split(':')
                    if len(parts) > 1:
                        density_str = parts[1].strip().rstrip('%')
                        stats["coding_density"] = float(density_str) / 100.0
                except (ValueError, IndexError):
                    pass
            elif 'gc content' in line.lower():
                try:
                    # Extract from line like "GC content: 67.3%"
                    parts = line.split(':')
                    if len(parts) > 1:
                        gc_str = parts[1].strip().rstrip('%')
                        stats["gc_content"] = float(gc_str) / 100.0
                except (ValueError, IndexError):
                    pass
                    
    except Exception as e:
        logger.warning(f"Failed to parse prodigal stats from {log_file}: {e}")
        
    return stats


def count_sequences_in_fasta(file_path: Path) -> int:
    """Count number of sequences in a FASTA file."""
    if not file_path.exists():
        return 0
        
    count = 0
    try:
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    count += 1
    except Exception as e:
        logger.warning(f"Failed to count sequences in {file_path}: {e}")
        
    return count


def validate_prodigal_outputs(genome_dir: Path, genome_id: str) -> Dict[str, Any]:
    """
    Validate prodigal output files and extract basic statistics.
    
    Returns:
        Dict containing validation results and file statistics
    """
    validation = {
        "output_files_exist": True,
        "files_non_empty": True,
        "protein_count": 0,
        "nucleotide_count": 0,
        "missing_files": [],
        "empty_files": []
    }
    
    # Check required output files
    protein_file = genome_dir / f"{genome_id}.faa"
    nucleotide_file = genome_dir / f"{genome_id}.genes.fna"
    
    for file_path, file_type in [(protein_file, "protein"), (nucleotide_file, "nucleotide")]:
        if not file_path.exists():
            validation["output_files_exist"] = False
            validation["missing_files"].append(str(file_path))
        elif file_path.stat().st_size == 0:
            validation["files_non_empty"] = False
            validation["empty_files"].append(str(file_path))
        else:
            # Count sequences
            seq_count = count_sequences_in_fasta(file_path)
            if file_type == "protein":
                validation["protein_count"] = seq_count
            else:
                validation["nucleotide_count"] = seq_count
    
    return validation


def run_prodigal_single(genome_info: Dict[str, Any], 
                       output_base_dir: Path,
                       mode: str = "single",
                       genetic_code: int = 11,
                       min_gene_length: int = 90,
                       include_nucleotides: bool = True) -> Dict[str, Any]:
    """
    Run prodigal on a single genome.
    
    Args:
        genome_info: Dictionary containing genome metadata from input manifest
        output_base_dir: Base output directory for all genomes
        mode: Prodigal mode (single, meta, train)
        genetic_code: Genetic code table number
        min_gene_length: Minimum gene length in bp
        include_nucleotides: Whether to output nucleotide gene sequences
        
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
    
    # Define output files
    protein_file = genome_output_dir / f"{genome_id}.faa"
    nucleotide_file = genome_output_dir / f"{genome_id}.genes.fna"
    log_file = genome_output_dir / "prodigal.log"
    
    # Initialize result dictionary
    result = {
        "genome_id": genome_id,
        "input_file": str(input_file),
        "output_dir": str(genome_output_dir),
        "execution_status": "failed",
        "execution_time_seconds": 0.0,
        "error_message": None,
        "prodigal_version": None,
        "statistics": {},
        "validation": {},
        "output_files": {}
    }
    
    try:
        # Build prodigal command
        cmd = [
            "prodigal",
            "-i", str(input_file),
            "-a", str(protein_file),  # protein output
            "-p", mode,               # procedure (single/meta/train)
            "-g", str(genetic_code),  # genetic code
            "-m"                      # output metadata to stderr
        ]
        
        # Add nucleotide output if requested
        if include_nucleotides:
            cmd.extend(["-d", str(nucleotide_file)])
            
        # Add minimum gene length if different from default
        if min_gene_length != 90:
            cmd.extend(["-t", str(min_gene_length)])
            
        # Execute prodigal
        with open(log_file, 'w') as log_f:
            process_result = subprocess.run(
                cmd,
                stdout=log_f,
                stderr=subprocess.STDOUT,
                text=True,
                timeout=300  # 5 minute timeout per genome
            )
            
        if process_result.returncode != 0:
            result["error_message"] = f"Prodigal failed with return code {process_result.returncode}"
            return result
            
        # Get prodigal version
        try:
            version_result = subprocess.run(
                ["prodigal", "-v"],
                capture_output=True,
                text=True,
                timeout=10
            )
            if version_result.returncode == 0:
                # Parse version from output
                version_lines = version_result.stderr.split('\n')
                for line in version_lines:
                    if 'Prodigal V' in line:
                        result["prodigal_version"] = line.strip()
                        break
        except Exception:
            pass  # Version detection is non-critical
            
        # Parse statistics from log
        result["statistics"] = parse_prodigal_stats(log_file)
        
        # Validate outputs
        result["validation"] = validate_prodigal_outputs(genome_output_dir, genome_id)
        
        # Record output files
        result["output_files"] = {
            "proteins": str(protein_file.relative_to(output_base_dir)),
        }
        if include_nucleotides:
            result["output_files"]["nucleotides"] = str(nucleotide_file.relative_to(output_base_dir))
            
        # Check if execution was successful
        if (result["validation"]["output_files_exist"] and 
            result["validation"]["files_non_empty"] and
            result["validation"]["protein_count"] > 0):
            result["execution_status"] = "success"
        else:
            result["execution_status"] = "failed"
            result["error_message"] = "Output validation failed"
            
    except subprocess.TimeoutExpired:
        result["error_message"] = "Prodigal execution timed out (>5 minutes)"
    except Exception as e:
        result["error_message"] = f"Unexpected error: {str(e)}"
        
    result["execution_time_seconds"] = round(time.time() - start_time, 2)
    return result


def process_genomes_parallel(genomes: List[Dict[str, Any]],
                           output_dir: Path,
                           max_workers: int,
                           **prodigal_kwargs) -> List[Dict[str, Any]]:
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
                run_prodigal_single,
                genome_info,
                output_dir,
                **prodigal_kwargs
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


def create_protein_symlinks(output_dir: Path, results: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Create symlinks to all .faa files in a centralized directory for easy access.
    
    Args:
        output_dir: Base output directory 
        results: List of processing results from prodigal runs
        
    Returns:
        Dict containing symlink creation statistics
    """
    symlink_stats = {
        "symlinks_created": 0,
        "symlinks_failed": 0,
        "symlink_directory": "",
        "failed_genomes": []
    }
    
    # Create the symlink directory
    symlink_dir = output_dir / "genomes" / "all_protein_symlinks"
    symlink_dir.mkdir(parents=True, exist_ok=True)
    symlink_stats["symlink_directory"] = str(symlink_dir)
    
    # Process only successful results
    successful_results = [r for r in results if r["execution_status"] == "success"]
    
    for result in successful_results:
        genome_id = result["genome_id"]
        
        try:
            # Construct path to the .faa file
            genome_dir = output_dir / "genomes" / genome_id
            faa_file = genome_dir / f"{genome_id}.faa"
            
            # Create symlink target path
            symlink_target = symlink_dir / f"{genome_id}.faa"
            
            # Remove existing symlink if it exists
            if symlink_target.exists() or symlink_target.is_symlink():
                symlink_target.unlink()
            
            # Create symlink (relative path for portability)
            relative_source = os.path.relpath(faa_file, symlink_dir)
            symlink_target.symlink_to(relative_source)
            
            symlink_stats["symlinks_created"] += 1
            
        except Exception as e:
            logger.warning(f"Failed to create symlink for {genome_id}: {e}")
            symlink_stats["symlinks_failed"] += 1
            symlink_stats["failed_genomes"].append({
                "genome_id": genome_id,
                "error": str(e)
            })
    
    return symlink_stats


def generate_summary_stats(results: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Generate aggregate statistics across all genomes."""
    successful = [r for r in results if r["execution_status"] == "success"]
    failed = [r for r in results if r["execution_status"] == "failed"]
    
    total_genes = sum(r.get("statistics", {}).get("genes_predicted", 0) for r in successful)
    total_proteins = sum(r.get("validation", {}).get("protein_count", 0) for r in successful)
    
    mean_execution_time = 0.0
    if successful:
        mean_execution_time = sum(r.get("execution_time_seconds", 0) for r in successful) / len(successful)
    
    summary = {
        "total_genomes": len(results),
        "successful": len(successful),
        "failed": len(failed),
        "success_rate": len(successful) / len(results) if results else 0.0,
        "total_genes_predicted": total_genes,
        "total_proteins": total_proteins,
        "mean_execution_time_seconds": round(mean_execution_time, 2),
        "total_execution_time_seconds": round(sum(r.get("execution_time_seconds", 0) for r in results), 2)
    }
    
    return summary


def run_prodigal(
    input_dir: Path = typer.Option(
        Path("data/stage00_prepared"),
        "--input-dir", "-i",
        help="Directory containing input genome assemblies"
    ),
    output_dir: Path = typer.Option(
        Path("data/stage03_prodigal"),
        "--output-dir", "-o",
        help="Output directory for Prodigal results"
    ),
    mode: str = typer.Option(
        "single",
        "--mode", "-m",
        help="Prodigal mode: single, meta, or train"
    ),
    genetic_code: int = typer.Option(
        11,
        "--genetic-code", "-g",
        help="Genetic code table (11 for bacteria/archaea)"
    ),
    min_gene_length: int = typer.Option(
        90,
        "--min-gene-length",
        help="Minimum gene length in nucleotides"
    ),
    max_workers: Optional[int] = typer.Option(
        None,
        "--max-workers", "-w",
        help="Maximum number of parallel workers (default: CPU count - 1)"
    ),
    include_nucleotides: bool = typer.Option(
        True,
        "--nucleotides/--no-nucleotides",
        help="Include nucleotide gene sequences in output"
    ),
    force: bool = typer.Option(
        False,
        "--force",
        help="Overwrite existing output directory"
    )
) -> None:
    """
    Run Prodigal gene prediction on genome assemblies.
    
    Predicts protein-coding sequences from prepared genome assemblies,
    generating protein sequences (.faa) and optionally nucleotide gene 
    sequences (.genes.fna). Processes genomes in parallel for efficiency.
    """
    console.print("[bold blue]Stage 3: Prodigal Gene Prediction[/bold blue]")
    
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
    (output_dir / "logs").mkdir(exist_ok=True)
    
    # Determine number of workers
    if max_workers is None:
        max_workers = max(1, multiprocessing.cpu_count() - 1)
    
    console.print(f"Using {max_workers} parallel workers")
    console.print(f"Mode: {mode}")
    console.print(f"Genetic code: {genetic_code}")
    console.print(f"Include nucleotides: {include_nucleotides}")
    
    # Process genomes in parallel
    console.print("\n[bold yellow]Processing genomes...[/bold yellow]")
    
    start_time = time.time()
    results = process_genomes_parallel(
        valid_genomes,
        output_dir,
        max_workers,
        mode=mode,
        genetic_code=genetic_code,
        min_gene_length=min_gene_length,
        include_nucleotides=include_nucleotides
    )
    total_time = time.time() - start_time
    
    # Create symlinks to all .faa files for easy access
    console.print("\n[bold yellow]Creating protein file symlinks...[/bold yellow]")
    symlink_stats = create_protein_symlinks(output_dir, results)
    
    # Generate summary statistics
    summary_stats = generate_summary_stats(results)
    
    # Create processing manifest
    manifest = {
        "version": "0.1.0",
        "stage": "stage03_prodigal",
        "timestamp": datetime.now().isoformat(),
        "input_manifest": str(manifest_file.absolute()),
        "execution_parameters": {
            "mode": mode,
            "genetic_code": genetic_code,
            "min_gene_length": min_gene_length,
            "include_nucleotides": include_nucleotides,
            "max_workers": max_workers
        },
        "summary": summary_stats,
        "symlink_stats": symlink_stats,
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
    console.print("\n[bold green]Stage 3 Results Summary[/bold green]")
    
    summary_table = Table()
    summary_table.add_column("Metric", style="cyan")
    summary_table.add_column("Value", style="magenta")
    
    summary_table.add_row("Input directory", str(input_dir))
    summary_table.add_row("Output directory", str(output_dir))
    summary_table.add_row("Total genomes", str(summary_stats["total_genomes"]))
    summary_table.add_row("Successful", str(summary_stats["successful"]))
    summary_table.add_row("Failed", str(summary_stats["failed"]))
    summary_table.add_row("Success rate", f"{summary_stats['success_rate']:.1%}")
    summary_table.add_row("Total genes predicted", f"{summary_stats['total_genes_predicted']:,}")  
    summary_table.add_row("Total proteins", f"{summary_stats['total_proteins']:,}")
    summary_table.add_row("Protein symlinks created", str(symlink_stats["symlinks_created"]))
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
        console.print(f"\n[bold green]✓ Stage 3 completed successfully![/bold green]")
        console.print(f"Protein sequences available in: {output_dir}/genomes/*/")
        if symlink_stats["symlinks_created"] > 0:
            console.print(f"All protein files (.faa) symlinked in: {symlink_stats['symlink_directory']}")
        
    logger.info(f"Stage 3 completed: {summary_stats['successful']} successful, {summary_stats['failed']} failed")


def main():
    """Entry point for standalone execution."""
    logging.basicConfig(level=logging.INFO)
    typer.run(run_prodigal)


if __name__ == "__main__":
    main()
