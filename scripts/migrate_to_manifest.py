#!/usr/bin/env python3
"""
Migration utility to create manifest.json from existing analysis outputs.

Scans an analysis directory for:
- findings.jsonl
- structures/*.pdb
- figures/*.png
- Database annotation counts

And builds a comprehensive manifest.json.

Usage:
    python scripts/migrate_to_manifest.py /path/to/data/dataset/
    python scripts/migrate_to_manifest.py /path/to/data/dataset/bennu.duckdb
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional


def find_database(dataset_dir: Path) -> Optional[Path]:
    """Find the DuckDB database in a dataset directory."""
    candidates = [
        dataset_dir / "bennu.duckdb",
        dataset_dir / "database.duckdb",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate

    # Search for any .duckdb file
    duckdb_files = list(dataset_dir.glob("*.duckdb"))
    if duckdb_files:
        return duckdb_files[0]

    return None


def scan_findings(dataset_dir: Path) -> dict:
    """Scan findings.jsonl and build findings summary."""
    findings_info = {
        "count": 0,
        "by_category": {},
        "high_priority": [],
    }
    proteins_with_findings: set[str] = set()

    # Check both exploration and survey directories
    for subdir in ["exploration", "survey"]:
        findings_path = dataset_dir / subdir / "findings.jsonl"
        if not findings_path.exists():
            continue

        print(f"  Found findings: {findings_path.relative_to(dataset_dir)}")

        try:
            with open(findings_path) as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue

                    finding = json.loads(line)
                    findings_info["count"] += 1

                    # Category tracking
                    cat = finding.get("category", "Uncategorized")
                    findings_info["by_category"][cat] = findings_info["by_category"].get(cat, 0) + 1

                    # High priority
                    if finding.get("priority") == "HIGH":
                        findings_info["high_priority"].append({
                            "id": finding.get("id"),
                            "title": finding.get("title"),
                        })

                    # Track proteins
                    for gene in finding.get("genes", []):
                        proteins_with_findings.add(str(gene))

        except Exception as e:
            print(f"    Warning: Error reading {findings_path}: {e}")

    return findings_info, len(proteins_with_findings)


def scan_structures(dataset_dir: Path) -> list[dict]:
    """Scan structures directory for PDB files."""
    structures = []

    structures_dir = dataset_dir / "structures"
    if not structures_dir.exists():
        return structures

    print(f"  Found structures directory: {structures_dir.relative_to(dataset_dir)}")

    for pdb_file in sorted(structures_dir.glob("*.pdb")):
        # Try to extract protein ID from filename
        # Common patterns: gene_123.pdb, protein_id.pdb, gene_123_nterm1024.pdb
        stem = pdb_file.stem
        protein_id = stem

        # Try to extract metrics from file if it has a sidecar JSON
        metrics_file = pdb_file.with_suffix(".json")
        metrics = {}
        if metrics_file.exists():
            try:
                metrics = json.loads(metrics_file.read_text())
            except Exception:
                pass

        structure_entry = {
            "protein_id": protein_id,
            "pdb_file": str(pdb_file.relative_to(dataset_dir)),
            "predicted_at": datetime.fromtimestamp(
                pdb_file.stat().st_mtime, tz=timezone.utc
            ).isoformat(),
        }

        # Add metrics if available
        if "plddt_mean" in metrics:
            structure_entry["plddt_mean"] = metrics["plddt_mean"]
        if "ptm" in metrics:
            structure_entry["ptm"] = metrics["ptm"]
        if "length" in metrics:
            structure_entry["length"] = metrics["length"]

        structures.append(structure_entry)

    print(f"    Found {len(structures)} PDB files")
    return structures


def scan_figures(dataset_dir: Path) -> list[dict]:
    """Scan for figure files."""
    figures = []

    # Check both exploration and survey figures directories
    for subdir in ["exploration", "survey"]:
        figures_dir = dataset_dir / subdir / "figures"
        if not figures_dir.exists():
            continue

        print(f"  Found figures directory: {figures_dir.relative_to(dataset_dir)}")

        for fig_file in sorted(figures_dir.glob("*.png")):
            # Try to infer figure type from filename
            stem = fig_file.stem.lower()

            figure_type = "unknown"
            if "neighborhood" in stem:
                figure_type = "neighborhood"
            elif "domain" in stem:
                figure_type = "domain"
            elif "heatmap" in stem:
                figure_type = "heatmap"
            elif "dendrogram" in stem:
                figure_type = "dendrogram"
            elif "ecotype" in stem:
                figure_type = "ecotype"

            # Try to extract center protein from neighborhood figures
            center_protein = None
            if figure_type == "neighborhood":
                # Pattern: something_gene123_neighborhood.png
                match = re.search(r"gene[_-]?(\d+)", stem)
                if match:
                    center_protein = f"gene_{match.group(1)}"

            figure_entry = {
                "path": str(fig_file.relative_to(dataset_dir)),
                "type": figure_type,
                "created_at": datetime.fromtimestamp(
                    fig_file.stat().st_mtime, tz=timezone.utc
                ).isoformat(),
            }

            # Use filename as title (cleaned up)
            title = fig_file.stem.replace("_", " ").replace("-", " ").title()
            figure_entry["title"] = title

            if center_protein:
                figure_entry["center_protein"] = center_protein

            figures.append(figure_entry)

        print(f"    Found {len([f for f in figures if subdir in f['path']])} figures")

    return figures


def scan_reports(dataset_dir: Path) -> list[dict]:
    """Scan for report files (PDF, MD)."""
    reports = []

    # Check reports/ directory (new standard) and root (legacy)
    search_dirs = [dataset_dir / "reports", dataset_dir]

    for search_dir in search_dirs:
        if not search_dir.exists():
            continue

        for pattern in ["*.pdf", "*.PDF"]:
            for report_file in sorted(search_dir.glob(pattern)):
                report_entry = {
                    "path": str(report_file.relative_to(dataset_dir)),
                    "generated": datetime.fromtimestamp(
                        report_file.stat().st_mtime, tz=timezone.utc
                    ).isoformat(),
                }
                reports.append(report_entry)

    # Deduplicate by path
    seen_paths = set()
    unique_reports = []
    for r in reports:
        if r["path"] not in seen_paths:
            seen_paths.add(r["path"])
            unique_reports.append(r)

    if unique_reports:
        print(f"  Found {len(unique_reports)} report files")

    return unique_reports


def get_annotation_counts(db_path: Path) -> dict:
    """Query database for annotation counts."""
    annotations = {}

    try:
        import duckdb

        conn = duckdb.connect(str(db_path), read_only=True)

        # Get total protein count
        result = conn.execute("SELECT COUNT(*) FROM proteins").fetchone()
        total_proteins = result[0] if result else 0

        # Get annotation counts by source
        result = conn.execute("""
            SELECT source, COUNT(DISTINCT protein_id), COUNT(*)
            FROM annotations
            GROUP BY source
        """).fetchall()

        for source, unique_proteins, total_annotations in result:
            coverage = unique_proteins / total_proteins if total_proteins > 0 else 0
            annotations[source.lower()] = {
                "count": total_annotations,
                "proteins_annotated": unique_proteins,
                "coverage": round(coverage, 3),
            }

        # Get predicate count
        try:
            result = conn.execute("SELECT COUNT(*) FROM protein_predicates").fetchone()
            if result and result[0] > 0:
                annotations["predicates"] = {"count": result[0]}
        except Exception:
            pass

        conn.close()
        print(f"  Annotation sources: {', '.join(sorted(annotations.keys()))}")

    except ImportError:
        print("  Warning: duckdb not installed, skipping annotation counts")
    except Exception as e:
        print(f"  Warning: Error reading database: {e}")

    return annotations


def get_dataset_info(db_path: Path) -> dict:
    """Get basic dataset info from database."""
    info = {
        "protein_count": 0,
        "genome_count": 0,
    }

    try:
        import duckdb

        conn = duckdb.connect(str(db_path), read_only=True)

        result = conn.execute("SELECT COUNT(*) FROM proteins").fetchone()
        info["protein_count"] = result[0] if result else 0

        try:
            result = conn.execute("SELECT COUNT(DISTINCT bin_id) FROM proteins WHERE bin_id IS NOT NULL").fetchone()
            info["genome_count"] = result[0] if result else 0
        except Exception:
            pass

        conn.close()

    except Exception as e:
        print(f"  Warning: Error getting dataset info: {e}")

    return info


def infer_exploration_status(
    findings_count: int,
    structures_count: int,
    figures_count: int,
) -> str:
    """Infer exploration status from outputs."""
    if findings_count == 0 and structures_count == 0 and figures_count == 0:
        return "not_started"
    elif findings_count > 0 or structures_count > 0:
        return "in_progress"
    return "in_progress"


def migrate_to_manifest(dataset_path: str | Path, dry_run: bool = False) -> dict:
    """
    Create manifest.json from existing analysis outputs.

    Args:
        dataset_path: Path to dataset directory or database file
        dry_run: If True, print manifest but don't write it

    Returns:
        The generated manifest dict
    """
    dataset_path = Path(dataset_path)

    # Determine dataset directory and database path
    if dataset_path.suffix == ".duckdb":
        db_path = dataset_path
        dataset_dir = dataset_path.parent
    else:
        dataset_dir = dataset_path
        db_path = find_database(dataset_dir)

    print(f"Dataset directory: {dataset_dir}")
    if db_path:
        print(f"Database: {db_path.name}")
    else:
        print("Warning: No database found")

    manifest_path = dataset_dir / "manifest.json"

    # Check for existing manifest
    if manifest_path.exists():
        print(f"Warning: manifest.json already exists at {manifest_path}")
        if not dry_run:
            backup_path = manifest_path.with_suffix(".json.bak")
            manifest_path.rename(backup_path)
            print(f"  Backed up to {backup_path.name}")

    # Build manifest
    print("\nScanning for analysis outputs...")

    # Dataset info
    dataset_info = {"name": dataset_dir.name, "database": db_path.name if db_path else None}
    if db_path:
        db_info = get_dataset_info(db_path)
        dataset_info["source_files"] = [{"type": "proteins", "count": db_info["protein_count"]}]
        if db_info["genome_count"] > 0:
            dataset_info["genome_count"] = db_info["genome_count"]

    # Annotations
    annotations = {}
    if db_path:
        annotations = get_annotation_counts(db_path)

    # Findings
    findings_info, proteins_with_findings = scan_findings(dataset_dir)

    # Structures
    structures = scan_structures(dataset_dir)

    # Figures
    figures = scan_figures(dataset_dir)

    # Reports
    reports = scan_reports(dataset_dir)

    # Infer status
    status = infer_exploration_status(
        findings_info["count"],
        len(structures),
        len(figures),
    )

    # Build manifest
    manifest = {
        "version": "1.0",
        "dataset": dataset_info,
        "annotations": annotations,
        "exploration": {
            "status": status,
            "phases_completed": [],
            "coverage": {
                "proteins_examined": 0,
                "proteins_with_findings": proteins_with_findings,
                "predicates_queried": [],
                "contigs_browsed": [],
            },
        },
        "findings": findings_info,
        "structures": {
            "predicted": len(structures),
            "proteins": structures,
        },
        "figures": figures,
        "reports": reports,
        "unexplored": {
            "proteins_without_findings": (
                dataset_info.get("source_files", [{}])[0].get("count", 0) - proteins_with_findings
                if dataset_info.get("source_files")
                else 0
            ),
            "predicates_not_queried": [],
            "suggested_next_steps": [],
        },
        "session_log": [
            {
                "timestamp": datetime.now(timezone.utc).isoformat(),
                "action": "migration",
                "summary": f"Manifest created from existing outputs. Found {findings_info['count']} findings, {len(structures)} structures, {len(figures)} figures.",
            }
        ],
    }

    # Summary
    print("\n" + "=" * 50)
    print("Migration Summary")
    print("=" * 50)
    print(f"  Dataset: {dataset_info['name']}")
    print(f"  Status: {status}")
    print(f"  Findings: {findings_info['count']}")
    print(f"  Structures: {len(structures)}")
    print(f"  Figures: {len(figures)}")
    print(f"  Reports: {len(reports)}")
    print(f"  Annotation sources: {len(annotations)}")

    if dry_run:
        print("\n[DRY RUN] Manifest would be written to:", manifest_path)
        print("\nManifest preview:")
        print(json.dumps(manifest, indent=2)[:2000] + "...")
    else:
        # Write manifest
        manifest_path.write_text(json.dumps(manifest, indent=2))
        print(f"\nManifest written to: {manifest_path}")

    return manifest


def main():
    parser = argparse.ArgumentParser(
        description="Create manifest.json from existing analysis outputs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python scripts/migrate_to_manifest.py data/my_dataset/
    python scripts/migrate_to_manifest.py data/my_dataset/bennu.duckdb
    python scripts/migrate_to_manifest.py data/my_dataset/ --dry-run
        """,
    )
    parser.add_argument(
        "dataset_path",
        help="Path to dataset directory or bennu.duckdb file",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print manifest without writing to disk",
    )

    args = parser.parse_args()

    if not Path(args.dataset_path).exists():
        print(f"Error: Path does not exist: {args.dataset_path}")
        sys.exit(1)

    migrate_to_manifest(args.dataset_path, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
