#!/usr/bin/env python3
"""
dbCAN CAZyme annotation module for carbohydrate-active enzyme detection.

This module processes protein sequences through dbCAN to identify CAZyme families
(GH, GT, PL, CE, AA) and integrates results into the genomic knowledge graph.

Usage:
    python -m src.ingest.dbcan_cazyme --input-dir data/stage03_prodigal --output-dir data/stage05b_dbcan
"""

import argparse
import json
import logging
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
from pydantic import BaseModel
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
console = Console()


class CAZymeAnnotation(BaseModel):
    """Represents a single CAZyme annotation result."""
    protein_id: str
    cazyme_family: str
    family_type: str  # GH, GT, PL, CE, AA, CBM
    evalue: float
    coverage: float
    start_pos: int
    end_pos: int
    hmm_length: int
    substrate_prediction: Optional[str] = None
    ec_number: Optional[str] = None


class CAZymeResult(BaseModel):
    """Complete CAZyme analysis results for a genome."""
    genome_id: str
    total_proteins: int
    cazyme_proteins: int
    annotations: List[CAZymeAnnotation]
    family_counts: Dict[str, int]
    processing_time: float


def run_dbcan_analysis(
    protein_file: Path,
    output_dir: Path,
    threads: int = 4,
    evalue_threshold: float = 1e-15,
    coverage_threshold: float = 0.35
) -> Optional[CAZymeResult]:
    """
    Run dbCAN analysis on a protein FASTA file.
    
    Args:
        protein_file: Path to protein FASTA file
        output_dir: Directory for dbCAN output
        threads: Number of threads for dbCAN
        evalue_threshold: E-value cutoff for HMM hits
        coverage_threshold: Coverage threshold for domain hits
        
    Returns:
        CAZymeResult object or None if analysis failed
    """
    import time
    start_time = time.time()
    
    genome_id = protein_file.stem.replace('.faa', '')
    genome_output_dir = output_dir / genome_id
    genome_output_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        # Enhanced: Use absolute paths for better compatibility
        abs_protein_file = protein_file.resolve()
        abs_output_dir = genome_output_dir.resolve()
        abs_db_dir = Path("data/dbcan_db").resolve()
        
        # Run dbCAN using CAZyme_annotation command with DIAMOND only
        cmd = [
            'run_dbcan',
            'CAZyme_annotation',
            '--input_raw_data', str(abs_protein_file),
            '--mode', 'protein',
            '--output_dir', str(abs_output_dir),
            '--db_dir', str(abs_db_dir),  # Use absolute path for database
            '--methods', 'diamond',  # Use only DIAMOND method (CAZy.dmnd works)
            '--threads', str(threads),
            '--e_value_threshold', str(evalue_threshold)  # DIAMOND parameter
        ]
        
        logger.info(f"Running dbCAN for {genome_id}: {' '.join(cmd)}")
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=3600  # 1 hour timeout
        )
        
        if result.returncode != 0:
            logger.error(f"dbCAN failed for {genome_id}: {result.stderr}")
            return None
            
        # Parse dbCAN overview.tsv output
        overview_file = genome_output_dir / "overview.tsv"
        if not overview_file.exists():
            logger.error(f"dbCAN overview file not found: {overview_file}")
            return None
            
        annotations = parse_dbcan_overview(overview_file)
        
        # Count total proteins
        total_proteins = count_proteins_in_fasta(protein_file)
        cazyme_proteins = len(set(ann.protein_id for ann in annotations))
        
        # Count families
        family_counts = {}
        for ann in annotations:
            family_counts[ann.cazyme_family] = family_counts.get(ann.cazyme_family, 0) + 1
            
        processing_time = time.time() - start_time
        
        return CAZymeResult(
            genome_id=genome_id,
            total_proteins=total_proteins,
            cazyme_proteins=cazyme_proteins,
            annotations=annotations,
            family_counts=family_counts,
            processing_time=processing_time
        )
        
    except subprocess.TimeoutExpired:
        logger.error(f"dbCAN analysis timed out for {genome_id}")
        return None
    except Exception as e:
        logger.error(f"Error running dbCAN for {genome_id}: {e}")
        return None


def parse_dbcan_overview(overview_file: Path) -> List[CAZymeAnnotation]:
    """
    Parse dbCAN overview.tsv file to extract CAZyme annotations.
    
    The overview file format:
    Gene ID	EC#	dbCAN_hmm	dbCAN_sub	DIAMOND	#ofTools	Recommend Results
    protein_1	-	GH13(1-295)	GH13	-	2	GH13
    """
    annotations = []
    
    try:
        df = pd.read_csv(overview_file, sep='\t')
        
        # Load substrate mapping for enhanced annotations
        substrate_mapping = load_cazyme_substrate_mapping()
        
        for _, row in df.iterrows():
            protein_id = row['Gene ID']
            
            # Try different columns for CAZyme family information
            # DIAMOND column should be prioritized for DIAMOND-only runs
            cazyme_result = None
            for col in ['DIAMOND', 'dbCAN_hmm', 'HMMER', 'Recommend Results']:
                if col in row and row[col] not in ['-', ''] and not pd.isna(row[col]):
                    cazyme_result = row[col]
                    break
            
            if not cazyme_result:
                continue
                
            # Handle multiple families in DIAMOND results (e.g., "CBM48+GH13_8")
            families = []
            if '+' in cazyme_result:
                # Split on + and process each family
                family_parts = cazyme_result.split('+')
                for part in family_parts:
                    families.append(part.strip())
            else:
                families = [cazyme_result]
            
            # Process each family found
            for family_result in families:
                # Parse format like "GH13(1-295)" or just "GH13"
                if '(' in family_result and ')' in family_result:
                    family = family_result.split('(')[0]
                    coords = family_result.split('(')[1].rstrip(')')
                    try:
                        start_pos, end_pos = map(int, coords.split('-'))
                    except ValueError:
                        start_pos, end_pos = 0, 0
                else:
                    family = family_result.strip()
                    start_pos, end_pos = 0, 0
                    
                # Determine family type
                family_type = get_cazyme_family_type(family)
                
                # Get EC number if available
                ec_number = row.get('EC#', '-')
                if ec_number == '-' or pd.isna(ec_number):
                    ec_number = None
                    
                # Enhanced: Add substrate prediction from mapping
                substrate_prediction = substrate_mapping.get(family, None)
                    
                annotation = CAZymeAnnotation(
                    protein_id=protein_id,
                    cazyme_family=family,
                    family_type=family_type,
                    evalue=1e-16,  # Default for overview format
                    coverage=0.8,  # Default for overview format
                    start_pos=start_pos,
                    end_pos=end_pos,
                    hmm_length=max(end_pos - start_pos + 1, 0) if end_pos > start_pos else 0,
                    substrate_prediction=substrate_prediction,
                    ec_number=ec_number
                )
                
                annotations.append(annotation)
            
    except Exception as e:
        logger.error(f"Error parsing dbCAN overview file {overview_file}: {e}")
        
    return annotations


def load_cazyme_substrate_mapping() -> Dict[str, str]:
    """
    Load CAZyme family to substrate mapping from fam-substrate-mapping.tsv.
    
    Returns:
        Dictionary mapping CAZyme family IDs to substrate descriptions
    """
    mapping = {}
    
    try:
        # Look for substrate mapping file
        mapping_file = Path("data/dbcan_db/fam-substrate-mapping.tsv")
        if not mapping_file.exists():
            logger.warning(f"CAZyme substrate mapping file not found: {mapping_file}")
            return {}
            
        if mapping_file.stat().st_size == 0:
            logger.warning("CAZyme substrate mapping file is empty")
            return {}
            
        df = pd.read_csv(mapping_file, sep='\t')
        
        # New format: Family, Substrate_high_level, Substrate_curated  
        if 'Family' in df.columns and 'Substrate_high_level' in df.columns:
            for _, row in df.iterrows():
                family = row['Family']
                substrate_high = row['Substrate_high_level']
                substrate_curated = row.get(' Substrate_curated', '') # Note the space
                
                if pd.notna(family) and pd.notna(substrate_high):
                    # Use high-level substrate, with curated details if available
                    substrate = str(substrate_high)
                    if pd.notna(substrate_curated) and str(substrate_curated).strip():
                        substrate += f" ({substrate_curated.strip()})"
                    mapping[str(family)] = substrate
        else:
            logger.warning(f"Unexpected format in substrate mapping file: {list(df.columns)}")
            
    except Exception as e:
        logger.warning(f"Could not load CAZyme substrate mapping: {e}")
        
    logger.info(f"Loaded substrate predictions for {len(mapping)} CAZyme families")
    return mapping


def get_cazyme_family_type(family: str) -> str:
    """Determine CAZyme family type from family name."""
    if family.startswith('GH'):
        return 'GH'  # Glycoside Hydrolase
    elif family.startswith('GT'):
        return 'GT'  # Glycosyltransferase
    elif family.startswith('PL'):
        return 'PL'  # Polysaccharide Lyase
    elif family.startswith('CE'):
        return 'CE'  # Carbohydrate Esterase
    elif family.startswith('AA'):
        return 'AA'  # Auxiliary Activity
    elif family.startswith('CBM'):
        return 'CBM'  # Carbohydrate-Binding Module
    else:
        return 'Unknown'


def count_proteins_in_fasta(fasta_file: Path) -> int:
    """Count number of protein sequences in FASTA file."""
    count = 0
    try:
        with open(fasta_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    count += 1
    except Exception as e:
        logger.error(f"Error counting proteins in {fasta_file}: {e}")
    return count


def run_single_dbcan_analysis(args: Tuple[Path, Path, int, float, float]) -> Optional[CAZymeResult]:
    """Wrapper function for parallel processing."""
    protein_file, output_dir, threads, evalue_threshold, coverage_threshold = args
    return run_dbcan_analysis(protein_file, output_dir, threads, evalue_threshold, coverage_threshold)


def run_dbcan_batch_analysis(
    input_dir: Path,
    output_dir: Path,
    max_workers: int = 2,
    threads_per_job: int = 4,
    evalue_threshold: float = 1e-15,
    coverage_threshold: float = 0.35
) -> Dict[str, CAZymeResult]:
    """
    Run dbCAN analysis on multiple genomes in parallel.
    
    Args:
        input_dir: Directory containing protein FASTA files
        output_dir: Output directory for results
        max_workers: Number of parallel dbCAN processes
        threads_per_job: Threads per dbCAN job
        evalue_threshold: E-value cutoff
        coverage_threshold: Coverage threshold
        
    Returns:
        Dictionary mapping genome IDs to CAZymeResult objects
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Find all protein FASTA files
    protein_files = list(input_dir.glob("*.faa"))
    if not protein_files:
        logger.warning(f"No .faa files found in {input_dir}")
        return {}
        
    logger.info(f"Found {len(protein_files)} protein files for dbCAN analysis")
    
    # Prepare arguments for parallel processing
    args_list = [
        (pfile, output_dir, threads_per_job, evalue_threshold, coverage_threshold)
        for pfile in protein_files
    ]
    
    results = {}
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        console=console
    ) as progress:
        
        task = progress.add_task("Running dbCAN analysis...", total=len(protein_files))
        
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            # Submit all jobs
            future_to_file = {
                executor.submit(run_single_dbcan_analysis, args): args[0]
                for args in args_list
            }
            
            # Collect results
            for future in as_completed(future_to_file):
                protein_file = future_to_file[future]
                try:
                    result = future.result()
                    if result:
                        results[result.genome_id] = result
                        logger.info(f"Completed dbCAN for {result.genome_id}: "
                                  f"{result.cazyme_proteins}/{result.total_proteins} CAZyme proteins")
                    else:
                        logger.warning(f"dbCAN analysis failed for {protein_file.name}")
                except Exception as e:
                    logger.error(f"Error processing {protein_file.name}: {e}")
                finally:
                    progress.advance(task)
    
    return results


def save_results(results: Dict[str, CAZymeResult], output_dir: Path) -> None:
    """Save dbCAN analysis results to JSON files."""
    
    # Save individual results
    for genome_id, result in results.items():
        result_file = output_dir / f"{genome_id}_cazyme_results.json"
        with open(result_file, 'w') as f:
            json.dump(result.dict(), f, indent=2)
    
    # Save summary
    summary = {
        'total_genomes': len(results),
        'total_proteins': sum(r.total_proteins for r in results.values()),
        'total_cazyme_proteins': sum(r.cazyme_proteins for r in results.values()),
        'genome_results': {
            genome_id: {
                'total_proteins': result.total_proteins,
                'cazyme_proteins': result.cazyme_proteins,
                'cazyme_families': len(result.family_counts),
                'processing_time': result.processing_time
            }
            for genome_id, result in results.items()
        }
    }
    
    summary_file = output_dir / "dbcan_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    logger.info(f"Saved dbCAN results for {len(results)} genomes to {output_dir}")


def create_processing_manifest(results: Dict[str, CAZymeResult], output_dir: Path) -> None:
    """Create processing manifest for pipeline integration."""
    manifest = {
        'stage': 'dbcan_cazyme',
        'version': '1.0.0',
        'timestamp': pd.Timestamp.now().isoformat(),
        'input_genomes': len(results),
        'successful_analyses': len(results),
        'total_proteins_analyzed': sum(r.total_proteins for r in results.values()),
        'total_cazyme_proteins': sum(r.cazyme_proteins for r in results.values()),
        'output_files': [
            f"{genome_id}_cazyme_results.json" for genome_id in results.keys()
        ] + ['dbcan_summary.json']
    }
    
    manifest_file = output_dir / "processing_manifest.json"
    with open(manifest_file, 'w') as f:
        json.dump(manifest, f, indent=2)


def main():
    """Main CLI function for dbCAN CAZyme analysis."""
    parser = argparse.ArgumentParser(
        description="Run dbCAN CAZyme analysis on protein sequences"
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Directory containing protein FASTA files (.faa)"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Output directory for dbCAN results"
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=2,
        help="Number of parallel dbCAN processes (default: 2)"
    )
    parser.add_argument(
        "--threads-per-job",
        type=int,
        default=4,
        help="Threads per dbCAN job (default: 4)"
    )
    parser.add_argument(
        "--evalue-threshold",
        type=float,
        default=1e-15,
        help="E-value threshold for HMM hits (default: 1e-15)"
    )
    parser.add_argument(
        "--coverage-threshold",
        type=float,
        default=0.35,
        help="Coverage threshold for domain hits (default: 0.35)"
    )
    
    args = parser.parse_args()
    
    console.print(f"[bold green]Starting dbCAN CAZyme Analysis[/bold green]")
    console.print(f"Input directory: {args.input_dir}")
    console.print(f"Output directory: {args.output_dir}")
    console.print(f"Max workers: {args.max_workers}")
    console.print(f"Threads per job: {args.threads_per_job}")
    
    # Run analysis
    results = run_dbcan_batch_analysis(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        max_workers=args.max_workers,
        threads_per_job=args.threads_per_job,
        evalue_threshold=args.evalue_threshold,
        coverage_threshold=args.coverage_threshold
    )
    
    if results:
        # Save results
        save_results(results, args.output_dir)
        create_processing_manifest(results, args.output_dir)
        
        # Print summary
        total_proteins = sum(r.total_proteins for r in results.values())
        total_cazymes = sum(r.cazyme_proteins for r in results.values())
        
        console.print(f"\n[bold green]dbCAN Analysis Complete![/bold green]")
        console.print(f"Genomes analyzed: {len(results)}")
        console.print(f"Total proteins: {total_proteins:,}")
        console.print(f"CAZyme proteins: {total_cazymes:,} ({total_cazymes/total_proteins*100:.1f}%)")
        
        # Show family distribution
        all_families = {}
        for result in results.values():
            for family, count in result.family_counts.items():
                all_families[family] = all_families.get(family, 0) + count
                
        console.print(f"\nTop CAZyme families:")
        for family, count in sorted(all_families.items(), key=lambda x: x[1], reverse=True)[:10]:
            console.print(f"  {family}: {count}")
    else:
        console.print("[bold red]No successful dbCAN analyses completed[/bold red]")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())