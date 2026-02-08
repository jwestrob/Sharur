#!/usr/bin/env python3
"""
Stage 0: Input Preparation
Validate and organize input genome assemblies for processing.
"""

import logging
import hashlib
import json
from pathlib import Path
from typing import List, Optional, Dict, Any
from datetime import datetime
import shutil
import re

import typer
from rich.console import Console
from rich.table import Table
from rich.progress import Progress, TaskID

console = Console()
logger = logging.getLogger(__name__)


def validate_fasta_format(file_path: Path) -> Dict[str, Any]:
    """
    Validate FASTA file format and extract basic statistics.
    
    Returns:
        Dict containing validation results and statistics
    """
    stats = {
        "is_valid": True,
        "error_message": None,
        "sequence_count": 0,
        "total_length": 0,
        "sequence_ids": [],
        "duplicate_ids": [],
        "invalid_characters": []
    }
    
    try:
        with open(file_path, 'r') as f:
            current_seq_id = None
            current_seq_length = 0
            seen_ids = set()
            line_number = 0
            
            for line in f:
                line_number += 1
                line = line.strip()
                
                if not line:  # Skip empty lines
                    continue
                
                if line.startswith('>'):
                    # Header line
                    if current_seq_id is not None:
                        # Finish previous sequence
                        stats["total_length"] += current_seq_length
                    
                    # Extract sequence ID
                    seq_id = line[1:].split()[0]  # Take first part after >
                    
                    if seq_id in seen_ids:
                        stats["duplicate_ids"].append(seq_id)
                        stats["is_valid"] = False
                    
                    seen_ids.add(seq_id)
                    stats["sequence_ids"].append(seq_id)
                    stats["sequence_count"] += 1
                    current_seq_id = seq_id
                    current_seq_length = 0
                    
                else:
                    # Sequence line
                    if current_seq_id is None:
                        stats["is_valid"] = False
                        stats["error_message"] = f"Sequence data before header at line {line_number}"
                        break
                    
                    # Check for invalid characters
                    valid_chars = set('ATCGNRYSWKMBVDHX-')  # IUPAC nucleotide codes
                    invalid_chars = set(line.upper()) - valid_chars
                    if invalid_chars:
                        stats["invalid_characters"].extend(list(invalid_chars))
                        stats["is_valid"] = False
                    
                    current_seq_length += len(line)
            
            # Finish last sequence
            if current_seq_id is not None:
                stats["total_length"] += current_seq_length
            
            # Check if file has any sequences (only if no previous error)
            if stats["sequence_count"] == 0 and not stats["error_message"]:
                stats["is_valid"] = False
                stats["error_message"] = "No sequences found in file"
    
    except Exception as e:
        stats["is_valid"] = False
        stats["error_message"] = f"Error reading file: {str(e)}"
    
    return stats


def calculate_file_checksum(file_path: Path) -> str:
    """Calculate MD5 checksum of file."""
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def find_genome_files(input_dir: Path, extensions: List[str]) -> List[Path]:
    """Find all genome files with specified extensions."""
    genome_files = []
    
    for ext in extensions:
        # Handle both .ext and ext formats
        if not ext.startswith('.'):
            ext = '.' + ext
        
        pattern = f"*{ext}"
        genome_files.extend(input_dir.glob(pattern))
    
    return sorted(genome_files)


def generate_genome_id(file_path: Path) -> str:
    """Generate genome ID from filename."""
    # Remove all common extensions
    name = file_path.name
    extensions = ['.fasta', '.fa', '.fna', '.fas', '.faa']
    
    for ext in extensions:
        if name.lower().endswith(ext.lower()):
            name = name[:-len(ext)]
            break
    
    # Clean up name (remove special characters, replace with underscores)
    genome_id = re.sub(r'[^a-zA-Z0-9_]', '_', name)
    genome_id = re.sub(r'_+', '_', genome_id)  # Replace multiple underscores
    genome_id = genome_id.strip('_')  # Remove leading/trailing underscores
    
    return genome_id


def prepare_inputs(
    input_dir: Path = typer.Option(
        Path("data/raw"),
        "--input-dir", "-i",
        help="Directory containing input genome assemblies"
    ),
    output_dir: Path = typer.Option(
        Path("data/stage00_prepared"),
        "--output-dir", "-o",
        help="Output directory for validated assemblies"
    ),
    file_extensions: List[str] = typer.Option(
        [".fasta", ".fa", ".fna"],
        "--extensions", "-e",
        help="File extensions to search for"
    ),
    validate_format: bool = typer.Option(
        True,
        "--validate/--no-validate",
        help="Validate FASTA format"
    ),
    copy_files: bool = typer.Option(
        False,
        "--copy/--symlink",
        help="Copy files instead of creating symlinks"
    ),
    force: bool = typer.Option(
        False,
        "--force",
        help="Overwrite existing output directory"
    )
) -> None:
    """
    Prepare input genome assemblies for pipeline processing.
    
    Validates FASTA file formats, checks for duplicates, generates processing
    manifest, and creates organized output directory structure.
    """
    console.print("[bold blue]Stage 0: Input Preparation[/bold blue]")
    
    # Validate input directory
    if not input_dir.exists():
        console.print(f"[red]Error: Input directory does not exist: {input_dir}[/red]")
        raise typer.Exit(1)
    
    if not input_dir.is_dir():
        console.print(f"[red]Error: Input path is not a directory: {input_dir}[/red]")
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
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Find genome files
    console.print(f"Searching for genome files in: {input_dir}")
    genome_files = find_genome_files(input_dir, file_extensions)
    
    if not genome_files:
        console.print(f"[red]Error: No genome files found with extensions: {file_extensions}[/red]")
        raise typer.Exit(1)
    
    console.print(f"Found {len(genome_files)} genome files")
    
    # Initialize manifest
    manifest = {
        "version": "0.1.0",
        "timestamp": datetime.now().isoformat(),
        "input_dir": str(input_dir.absolute()),
        "output_dir": str(output_dir.absolute()),
        "file_extensions": file_extensions,
        "validate_format": validate_format,
        "copy_files": copy_files,
        "genomes": []
    }
    
    # Process each genome file
    valid_genomes = 0
    invalid_genomes = 0
    
    with Progress() as progress:
        task = progress.add_task("Processing genomes...", total=len(genome_files))
        
        for genome_file in genome_files:
            progress.update(task, description=f"Processing {genome_file.name}")
            
            # Generate genome ID
            genome_id = generate_genome_id(genome_file)
            
            # Initialize genome entry
            genome_entry = {
                "filename": genome_file.name,
                "genome_id": genome_id,
                "file_path": str(genome_file.absolute()),
                "file_size": genome_file.stat().st_size,
                "checksum": calculate_file_checksum(genome_file),
                "format_valid": True,
                "validation_errors": []
            }
            
            # Validate FASTA format if requested
            if validate_format:
                validation_results = validate_fasta_format(genome_file)
                genome_entry["format_valid"] = validation_results["is_valid"]
                genome_entry["sequence_count"] = validation_results["sequence_count"]
                genome_entry["total_length"] = validation_results["total_length"]
                genome_entry["sequence_ids"] = validation_results["sequence_ids"]
                
                if not validation_results["is_valid"]:
                    errors = []
                    if validation_results["error_message"]:
                        errors.append(validation_results["error_message"])
                    if validation_results["duplicate_ids"]:
                        errors.append(f"Duplicate sequence IDs: {validation_results['duplicate_ids']}")
                    if validation_results["invalid_characters"]:
                        errors.append(f"Invalid characters: {set(validation_results['invalid_characters'])}")
                    
                    genome_entry["validation_errors"] = errors
                    invalid_genomes += 1
                else:
                    valid_genomes += 1
            
            # Create output file (symlink or copy)
            output_file = output_dir / genome_file.name
            
            try:
                if copy_files:
                    shutil.copy2(genome_file, output_file)
                else:
                    output_file.symlink_to(genome_file.absolute())
                
                genome_entry["output_path"] = str(output_file.absolute())
                genome_entry["linked_successfully"] = True
                
            except Exception as e:
                genome_entry["linked_successfully"] = False
                genome_entry["link_error"] = str(e)
                console.print(f"[red]Warning: Failed to link {genome_file.name}: {e}[/red]")
            
            manifest["genomes"].append(genome_entry)
            progress.advance(task)
    
    # Save manifest
    manifest_file = output_dir / "processing_manifest.json"
    with open(manifest_file, 'w') as f:
        json.dump(manifest, f, indent=2)
    
    # Display summary
    console.print("\n[bold green]Input Preparation Summary[/bold green]")
    
    summary_table = Table()
    summary_table.add_column("Metric", style="cyan")
    summary_table.add_column("Value", style="magenta")
    
    summary_table.add_row("Input directory", str(input_dir))
    summary_table.add_row("Output directory", str(output_dir))
    summary_table.add_row("Total files found", str(len(genome_files)))
    
    if validate_format:
        summary_table.add_row("Valid genomes", str(valid_genomes))
        summary_table.add_row("Invalid genomes", str(invalid_genomes))
    
    summary_table.add_row("File operation", "Copy" if copy_files else "Symlink")
    summary_table.add_row("Manifest file", str(manifest_file))
    
    console.print(summary_table)
    
    # Show validation errors if any
    if validate_format and invalid_genomes > 0:
        console.print("\n[bold red]Validation Errors:[/bold red]")
        for genome in manifest["genomes"]:
            if not genome["format_valid"]:
                console.print(f"[red]• {genome['filename']}:[/red]")
                for error in genome["validation_errors"]:
                    console.print(f"  - {error}")
    
    console.print(f"\n[bold green]✓ Stage 0 completed successfully![/bold green]")
    logger.info(f"Stage 0 completed: {valid_genomes} valid, {invalid_genomes} invalid genomes")


def main():
    """Entry point for standalone execution."""
    logging.basicConfig(level=logging.INFO)
    logger.info("Stage 0 stub")
    typer.run(prepare_inputs)


if __name__ == "__main__":
    main()
