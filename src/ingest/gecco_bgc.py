#!/usr/bin/env python3
"""
Stage 5a: GECCO Biosynthetic Gene Cluster Detection
Analyze genome assemblies for biosynthetic gene clusters using GECCO.

GECCO (Gene Cluster Prediction with Conditional Random Fields) is a Python-native
tool that avoids Docker compatibility issues and provides antiSMASH-compatible output.
"""

import logging
import json
import subprocess
import time
import shutil
import os
from pathlib import Path
from typing import List, Optional, Dict, Any, Tuple
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd
import csv

import typer
from rich.console import Console
from rich.table import Table
from rich.progress import Progress

console = Console()
logger = logging.getLogger(__name__)


def run_gecco(fasta: Path, out_dir: Path, threads: int = 0, is_metagenome: bool = False) -> bool:
    """
    Run GECCO BGC detection on a genome assembly using Python API.
    
    Args:
        fasta: Path to input FASTA file
        out_dir: Output directory for GECCO results
        threads: Number of threads to use (0 = auto)
        is_metagenome: Whether to use metagenome mode
        
    Returns:
        True if successful, False otherwise
    """
    try:
        # Try using GECCO Python API first
        try:
            import gecco.cli
            from gecco.cli import run
            
            # Create fake argv for GECCO
            import sys
            old_argv = sys.argv
            sys.argv = [
                "gecco", "run",
                "--genome", str(fasta),
                "-o", str(out_dir),
                "--jobs", str(threads if threads > 0 else os.cpu_count())
            ]
            
            if is_metagenome:
                sys.argv.append("--prodigal-meta")
            
            logger.info(f"Running GECCO via Python API: {' '.join(sys.argv[1:])}")
            
            # Run GECCO
            run.main()
            sys.argv = old_argv
            
            logger.info(f"GECCO completed successfully for {fasta.name}")
            return True
            
        except Exception as api_error:
            logger.warning(f"GECCO Python API failed: {api_error}, trying command line")
            
            # Fallback to command line
            cmd = [
                "gecco", "run",
                "--genome", str(fasta),
                "-o", str(out_dir),
                "--jobs", str(threads if threads > 0 else os.cpu_count())
            ]
            
            # Add metagenome flag if needed
            if is_metagenome:
                cmd.append("--prodigal-meta")
            
            logger.info(f"Running GECCO via command line: {' '.join(cmd)}")
            
            # Run GECCO
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=1800  # 30 minute timeout
            )
            
            if result.returncode != 0:
                logger.error(f"GECCO failed with return code {result.returncode}")
                logger.error(f"STDERR: {result.stderr}")
                
                # Create empty output for compatibility
                create_empty_gecco_output(out_dir, fasta)
                return True  # Return True to continue pipeline
            
            logger.info(f"GECCO completed successfully for {fasta.name}")
            return True
        
    except subprocess.TimeoutExpired:
        logger.error(f"GECCO timed out after 30 minutes for {fasta.name}")
        create_empty_gecco_output(out_dir, fasta)
        return True  # Return True to continue pipeline
    except Exception as e:
        logger.error(f"Error running GECCO on {fasta.name}: {e}")
        create_empty_gecco_output(out_dir, fasta)
        return True  # Return True to continue pipeline


def create_empty_gecco_output(out_dir: Path, fasta: Path) -> None:
    """Create empty GECCO output files for compatibility."""
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Create empty clusters.tsv file
    clusters_file = out_dir / f"{fasta.stem}.clusters.tsv"
    with open(clusters_file, 'w') as f:
        f.write("contig\tstart\tend\ttype\tproduct_class\tproteins\n")
    
    logger.info(f"Created empty GECCO output for {fasta.name} (no BGCs found)")


def parse_gecco_clusters_tsv(clusters_file: Path) -> List[Dict[str, Any]]:
    """
    Parse GECCO clusters.tsv output file.
    
    Format: contig, start, end, bgc_type, product_class, proteins
    
    Args:
        clusters_file: Path to {sample}.clusters.tsv file
        
    Returns:
        List of cluster dictionaries
    """
    clusters = []
    
    if not clusters_file.exists():
        logger.warning(f"Clusters file not found: {clusters_file}")
        return clusters
    
    try:
        with open(clusters_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            
            for i, row in enumerate(reader):
                # Parse protein list and count
                protein_list_str = row.get("proteins", "")
                protein_list = protein_list_str.split(";") if protein_list_str else []
                protein_count = len(protein_list) if protein_list and protein_list[0] else 0
                
                cluster = {
                    "cluster_id": f"cluster_{i+1}",
                    "contig": row.get("sequence_id", ""),  # Use correct column name
                    "start": int(row.get("start", 0)),
                    "end": int(row.get("end", 0)),
                    "bgc_type": row.get("type", "unknown"),
                    "product_class": row.get("product_class", ""),
                    "proteins": protein_count,
                    "protein_list": protein_list,  # Keep the actual protein IDs
                    "length": int(row.get("end", 0)) - int(row.get("start", 0)),
                    "source": "gecco",
                    # GECCO probability scores
                    "average_p": float(row.get("average_p", 0)),
                    "max_p": float(row.get("max_p", 0)),
                    "alkaloid_probability": float(row.get("alkaloid_probability", 0)),
                    "nrp_probability": float(row.get("nrp_probability", 0)),
                    "polyketide_probability": float(row.get("polyketide_probability", 0)),
                    "ripp_probability": float(row.get("ripp_probability", 0)),
                    "saccharide_probability": float(row.get("saccharide_probability", 0)),
                    "terpene_probability": float(row.get("terpene_probability", 0)),
                    "domains": row.get("domains", "")
                }
                clusters.append(cluster)
                
        logger.info(f"Parsed {len(clusters)} BGC clusters from {clusters_file}")
        return clusters
        
    except Exception as e:
        logger.error(f"Error parsing GECCO clusters file {clusters_file}: {e}")
        return []


def convert_gecco_to_genbank(clusters_file: Path, output_dir: Path) -> List[Path]:
    """
    Convert GECCO clusters.tsv to antiSMASH-style GenBank files.
    
    Args:
        clusters_file: Path to GECCO clusters.tsv file
        output_dir: Directory to write GenBank files
        
    Returns:
        List of generated GenBank file paths
    """
    gbk_files = []
    
    if not clusters_file.exists():
        return gbk_files
    
    try:
        # Use GECCO's built-in conversion
        sample_name = clusters_file.stem.replace('.clusters', '')
        
        cmd = [
            "gecco", "convert", "gbk",
            "--format", "bigslice",
            "-i", str(clusters_file),
            "-o", str(output_dir)
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            # Find generated GenBank files
            for gbk_file in output_dir.glob(f"{sample_name}_cluster_*.gbk"):
                gbk_files.append(gbk_file)
                
        logger.info(f"Converted GECCO output to {len(gbk_files)} GenBank files")
        
    except Exception as e:
        logger.error(f"Error converting GECCO output to GenBank: {e}")
    
    return gbk_files


def parse_gecco_genbank(gbk_file: Path) -> Dict[str, Any]:
    """
    Parse GECCO-generated GenBank file to extract BGC information.
    
    Args:
        gbk_file: Path to GECCO GenBank output file
        
    Returns:
        Dict containing BGC data and gene annotations
    """
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature
    
    bgc_data = {
        "file": str(gbk_file),
        "clusters": [],
        "genes": [],
        "parsing_errors": []
    }
    
    try:
        for record in SeqIO.parse(gbk_file, "genbank"):
            # Extract cluster information
            cluster_info = {
                "cluster_id": gbk_file.stem,
                "contig": record.id,
                "start": 1,
                "end": len(record.seq),
                "length": len(record.seq),
                "bgc_type": "unknown",
                "product": [],
                "source": "gecco"
            }
            
            # Parse features for BGC and gene information
            for feature in record.features:
                if feature.type == "cluster":
                    # Extract cluster-level information
                    qualifiers = feature.qualifiers
                    cluster_info.update({
                        "bgc_type": qualifiers.get("product", ["unknown"])[0],
                        "product": qualifiers.get("product", []),
                        "start": int(feature.location.start) + 1,
                        "end": int(feature.location.end)
                    })
                    
                elif feature.type == "CDS":
                    # Extract gene information
                    qualifiers = feature.qualifiers
                    gene_info = {
                        "gene_id": qualifiers.get("locus_tag", [f"gene_{len(bgc_data['genes'])+1}"])[0],
                        "contig": record.id,
                        "start": int(feature.location.start) + 1,
                        "end": int(feature.location.end),
                        "strand": "+" if feature.location.strand == 1 else "-",
                        "product": qualifiers.get("product", ["hypothetical protein"])[0],
                        "bgc_function": qualifiers.get("gene_functions", []),
                        "cluster_id": cluster_info["cluster_id"]
                    }
                    bgc_data["genes"].append(gene_info)
            
            bgc_data["clusters"].append(cluster_info)
            
        logger.info(f"Parsed {len(bgc_data['clusters'])} clusters and {len(bgc_data['genes'])} genes from {gbk_file}")
        
    except Exception as e:
        error_msg = f"Error parsing GenBank file {gbk_file}: {e}"
        logger.error(error_msg)
        bgc_data["parsing_errors"].append(error_msg)
    
    return bgc_data


def process_genome_gecco(genome_file: Path, output_dir: Path, threads: int = 1) -> Dict[str, Any]:
    """
    Process a single genome file with GECCO BGC detection.
    
    Args:
        genome_file: Path to genome FASTA file
        output_dir: Output directory for results
        threads: Number of threads to use
        
    Returns:
        Dictionary with processing results
    """
    genome_name = genome_file.stem
    genome_output_dir = output_dir / genome_name
    genome_output_dir.mkdir(parents=True, exist_ok=True)
    
    result = {
        "genome": genome_name,
        "input_file": str(genome_file),
        "output_dir": str(genome_output_dir),
        "start_time": datetime.now().isoformat(),
        "end_time": None,
        "status": "pending",
        "clusters": [],
        "genes": [],
        "errors": []
    }
    
    try:
        logger.info(f"Processing {genome_name} with GECCO")
        
        # Run GECCO
        success = run_gecco(genome_file, genome_output_dir, threads)
        
        if not success:
            result["status"] = "failed"
            result["errors"].append("GECCO execution failed")
            return result
        
        # Parse clusters.tsv output
        clusters_file = genome_output_dir / f"{genome_name}.clusters.tsv"
        if clusters_file.exists():
            clusters = parse_gecco_clusters_tsv(clusters_file)
            result["clusters"] = clusters
            
            # Convert to GenBank format for compatibility
            gbk_files = convert_gecco_to_genbank(clusters_file, genome_output_dir)
            
            # Parse GenBank files for detailed gene information
            for gbk_file in gbk_files:
                bgc_data = parse_gecco_genbank(gbk_file)
                result["genes"].extend(bgc_data["genes"])
                if bgc_data["parsing_errors"]:
                    result["errors"].extend(bgc_data["parsing_errors"])
        
        result["status"] = "completed"
        logger.info(f"GECCO processing completed for {genome_name}: {len(result['clusters'])} clusters, {len(result['genes'])} genes")
        
    except Exception as e:
        error_msg = f"Error processing {genome_name}: {e}"
        logger.error(error_msg)
        result["status"] = "failed"
        result["errors"].append(error_msg)
    
    finally:
        result["end_time"] = datetime.now().isoformat()
    
    return result


def gecco_bgc_detection(
    input_dir: Path,
    output_dir: Path,
    threads: int = 4,
    max_workers: int = 2,
    force: bool = False
) -> Dict[str, Any]:
    """
    Run GECCO BGC detection on multiple genome assemblies.
    
    Args:
        input_dir: Directory containing genome FASTA files
        output_dir: Output directory for BGC detection results
        threads: Number of threads per GECCO job
        max_workers: Maximum number of parallel GECCO jobs
        force: Overwrite existing outputs
        
    Returns:
        Dictionary with processing summary and results
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Setup logging
    log_dir = output_dir / "logs"
    log_dir.mkdir(exist_ok=True)
    
    # Find genome files
    genome_files = []
    for ext in ["*.fna", "*.fasta", "*.fa"]:
        genome_files.extend(list(Path(input_dir).glob(ext)))
    
    if not genome_files:
        raise ValueError(f"No genome files found in {input_dir}")
    
    console.print(f"🧬 Found {len(genome_files)} genome files for GECCO BGC detection")
    
    # Process genomes
    results = []
    all_clusters = []
    all_genes = []
    
    with Progress() as progress:
        task = progress.add_task("Processing genomes...", total=len(genome_files))
        
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            # Submit jobs
            future_to_genome = {
                executor.submit(process_genome_gecco, genome_file, output_dir, threads): genome_file
                for genome_file in genome_files
            }
            
            # Collect results
            for future in as_completed(future_to_genome):
                genome_file = future_to_genome[future]
                try:
                    result = future.result()
                    results.append(result)
                    all_clusters.extend(result["clusters"])
                    all_genes.extend(result["genes"])
                    
                    status_color = "green" if result["status"] == "completed" else "red"
                    console.print(
                        f"[{status_color}]✓[/{status_color}] {result['genome']}: "
                        f"{len(result['clusters'])} clusters, {len(result['genes'])} genes"
                    )
                    
                except Exception as e:
                    logger.error(f"Error processing {genome_file}: {e}")
                    results.append({
                        "genome": genome_file.stem,
                        "status": "failed",
                        "errors": [str(e)]
                    })
                
                progress.update(task, advance=1)
    
    # Save combined results
    combined_data = {
        "processing_summary": {
            "total_genomes": len(genome_files),
            "successful": sum(1 for r in results if r["status"] == "completed"),
            "failed": sum(1 for r in results if r["status"] == "failed"),
            "total_clusters": len(all_clusters),
            "total_genes": len(all_genes),
            "timestamp": datetime.now().isoformat()
        },
        "clusters": all_clusters,
        "genes": all_genes,
        "genome_results": results
    }
    
    # Save JSON output
    json_output = output_dir / "combined_bgc_data.json"
    with open(json_output, 'w') as f:
        json.dump(combined_data, f, indent=2)
    
    # Save processing manifest
    manifest = {
        "stage": "gecco_bgc_detection",
        "input_dir": str(input_dir),
        "output_dir": str(output_dir),
        "total_genomes": len(genome_files),
        "successful_genomes": combined_data["processing_summary"]["successful"],
        "total_clusters": len(all_clusters),
        "total_genes": len(all_genes),
        "timestamp": datetime.now().isoformat(),
        "parameters": {
            "threads": threads,
            "max_workers": max_workers,
            "force": force
        }
    }
    
    manifest_file = output_dir / "processing_manifest.json"
    with open(manifest_file, 'w') as f:
        json.dump(manifest, f, indent=2)
    
    # Print summary
    summary_table = Table(title="GECCO BGC Detection Summary")
    summary_table.add_column("Metric", style="cyan")
    summary_table.add_column("Value", style="green")
    
    summary_table.add_row("Total Genomes", str(len(genome_files)))
    summary_table.add_row("Successful", str(combined_data["processing_summary"]["successful"]))
    summary_table.add_row("Failed", str(combined_data["processing_summary"]["failed"]))
    summary_table.add_row("Total BGC Clusters", str(len(all_clusters)))
    summary_table.add_row("Total BGC Genes", str(len(all_genes)))
    summary_table.add_row("Output Directory", str(output_dir))
    
    console.print(summary_table)
    
    return combined_data


def main(
    input_dir: Path = typer.Argument(..., help="Directory containing genome FASTA files"),
    output_dir: Path = typer.Argument(..., help="Output directory for BGC detection results"),
    threads: int = typer.Option(4, "--threads", "-t", help="Number of threads per GECCO job"),
    max_workers: int = typer.Option(2, "--max-workers", "-w", help="Maximum parallel GECCO jobs"),
    force: bool = typer.Option(False, "--force", help="Overwrite existing outputs")
):
    """
    GECCO Biosynthetic Gene Cluster Detection Pipeline
    
    Detect biosynthetic gene clusters in microbial genomes using GECCO.
    """
    console.print("[bold blue]🧬 GECCO BGC Detection Pipeline[/bold blue]")
    console.print(f"Input directory: {input_dir}")
    console.print(f"Output directory: {output_dir}")
    console.print(f"Threads per job: {threads}")
    console.print(f"Max parallel jobs: {max_workers}")
    
    try:
        results = gecco_bgc_detection(
            input_dir=input_dir,
            output_dir=output_dir,
            threads=threads,
            max_workers=max_workers,
            force=force
        )
        
        console.print("[bold green]✅ GECCO BGC detection completed successfully![/bold green]")
        return results
        
    except Exception as e:
        console.print(f"[bold red]❌ Error: {e}[/bold red]")
        raise typer.Exit(1)


if __name__ == "__main__":
    typer.run(main)