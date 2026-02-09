"""
Analysis Manifest system for Bennu.

Provides a single-file manifest that captures the complete state of an analysis,
enabling session continuity, report generation, and reproducibility.

The manifest tracks:
- Dataset metadata (source files, protein counts)
- Annotation coverage (PFAM, KEGG, CAZy, etc.)
- Exploration status and coverage
- Findings (references to findings.jsonl)
- Structure predictions with metrics
- Figures with legends
- Session history (append-only log)
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from bennu.storage.duckdb_store import DuckDBStore


def _now_iso() -> str:
    """Return current UTC timestamp in ISO format."""
    return datetime.now(timezone.utc).isoformat()


class AnalysisManifest:
    """
    Manages the analysis manifest for a Bennu dataset.

    The manifest is stored as a JSON file alongside the database and provides
    a comprehensive view of analysis state, enabling:
    - Session continuity (resume after break)
    - Report generation (figure legends, findings)
    - Reproducibility (session log with timestamps)

    Location: `{db_directory}/manifest.json`
    """

    VERSION = "1.0"

    def __init__(self, db_path: str | Path, store: Optional["DuckDBStore"] = None):
        """
        Initialize manifest manager.

        Args:
            db_path: Path to the DuckDB database file
            store: Optional DuckDBStore instance for annotation queries
        """
        self.db_path = Path(db_path) if db_path else None
        self.path = self.db_path.parent / "manifest.json" if self.db_path else None
        self._store = store
        self.data = self._load_or_create()

    def _load_or_create(self) -> dict:
        """Load existing manifest or create new one."""
        if self.path and self.path.exists():
            try:
                return json.loads(self.path.read_text())
            except (json.JSONDecodeError, IOError):
                # Corrupted manifest - create fresh
                pass

        # Create new manifest with minimal structure
        return self._create_empty_manifest()

    def _create_empty_manifest(self) -> dict:
        """Create a new empty manifest with required fields."""
        dataset_name = self.db_path.parent.name if self.db_path else "unknown"

        return {
            "version": self.VERSION,
            "dataset": {
                "name": dataset_name,
                "created": _now_iso(),
                "database": self.db_path.name if self.db_path else None,
                "source_files": [],
            },
            "annotations": {},
            "exploration": {
                "status": "not_started",
                "phases_completed": [],
                "coverage": {
                    "proteins_examined": 0,
                    "proteins_with_findings": 0,
                    "predicates_queried": [],
                    "contigs_browsed": [],
                },
            },
            "findings": {
                "count": 0,
                "by_category": {},
                "high_priority": [],
            },
            "structures": {
                "predicted": 0,
                "proteins": [],
            },
            "figures": [],
            "hypotheses": {
                "count": 0,
                "active": 0,
                "registry_path": "exploration/hypotheses.json",
            },
            "reports": [],
            "unexplored": {
                "proteins_without_findings": 0,
                "predicates_not_queried": [],
                "suggested_next_steps": [],
            },
            "session_log": [],
        }

    # ------------------------------------------------------------------ #
    # Annotation tracking
    # ------------------------------------------------------------------ #

    def update_annotations(self, store: Optional["DuckDBStore"] = None) -> None:
        """
        Update annotation counts from database.

        Queries the database for annotation counts by source and updates
        the manifest with coverage statistics.

        Args:
            store: DuckDBStore instance (uses cached store if not provided)
        """
        store = store or self._store
        if not store:
            return

        try:
            # Get total protein count
            result = store.execute("SELECT COUNT(*) FROM proteins")
            total_proteins = result[0][0] if result else 0

            # Get annotation counts by source
            result = store.execute("""
                SELECT source, COUNT(DISTINCT protein_id), COUNT(*)
                FROM annotations
                GROUP BY source
            """)

            for source, unique_proteins, total_annotations in result:
                coverage = unique_proteins / total_proteins if total_proteins > 0 else 0
                self.data["annotations"][source.lower()] = {
                    "count": total_annotations,
                    "proteins_annotated": unique_proteins,
                    "coverage": round(coverage, 3),
                    "updated": _now_iso(),
                }

            # Get predicate count if available
            try:
                result = store.execute("SELECT COUNT(*) FROM protein_predicates")
                predicate_count = result[0][0] if result else 0
                if predicate_count > 0:
                    self.data["annotations"]["predicates"] = {
                        "count": predicate_count,
                        "generated": _now_iso(),
                    }
            except Exception:
                pass  # Table may not exist

            # Update source files from proteins table if not set
            if not self.data["dataset"]["source_files"]:
                result = store.execute("SELECT COUNT(*) FROM proteins")
                protein_count = result[0][0] if result else 0
                if protein_count > 0:
                    self.data["dataset"]["source_files"] = [
                        {"type": "proteins", "count": protein_count}
                    ]

        except Exception as e:
            # Log error but don't fail
            self.log_session("error", f"Failed to update annotations: {e}")

    # ------------------------------------------------------------------ #
    # Session logging
    # ------------------------------------------------------------------ #

    def log_session(self, action: str, summary: str) -> None:
        """
        Append entry to session log.

        Args:
            action: Type of action (e.g., "survey", "exploration", "structure_prediction")
            summary: Brief description of what was done
        """
        entry = {
            "timestamp": _now_iso(),
            "action": action,
            "summary": summary,
        }
        self.data["session_log"].append(entry)

    # ------------------------------------------------------------------ #
    # Structure tracking
    # ------------------------------------------------------------------ #

    def add_structure(
        self,
        protein_id: str,
        pdb_path: str,
        metrics: Optional[dict] = None,
        foldseek_hits: Optional[int] = None,
        best_hit: Optional[dict] = None,
        interpretation: Optional[str] = None,
    ) -> None:
        """
        Record a structure prediction.

        Args:
            protein_id: Protein identifier
            pdb_path: Path to PDB file (relative to dataset directory)
            metrics: Dict with plddt_mean, ptm, etc.
            foldseek_hits: Number of Foldseek hits (if searched)
            best_hit: Best Foldseek hit details
            interpretation: Human interpretation of structure
        """
        metrics = metrics or {}

        # Make path relative to dataset directory if possible
        if self.db_path:
            try:
                pdb_path = str(Path(pdb_path).relative_to(self.db_path.parent))
            except ValueError:
                pass  # Keep absolute path if not relative

        structure_entry = {
            "protein_id": protein_id,
            "pdb_file": pdb_path,
            "predicted_at": _now_iso(),
        }

        # Add optional fields
        if "length" in metrics:
            structure_entry["length"] = metrics["length"]
        if "plddt_mean" in metrics:
            structure_entry["plddt_mean"] = metrics["plddt_mean"]
        if "ptm" in metrics:
            structure_entry["ptm"] = metrics["ptm"]
        if "plddt_by_region" in metrics:
            structure_entry["plddt_by_region"] = metrics["plddt_by_region"]
        if foldseek_hits is not None:
            structure_entry["foldseek_hits"] = foldseek_hits
        if best_hit:
            structure_entry["best_hit"] = best_hit
        if interpretation:
            structure_entry["interpretation"] = interpretation

        # Check if protein already has a structure entry
        existing_idx = None
        for i, s in enumerate(self.data["structures"]["proteins"]):
            if s.get("protein_id") == protein_id:
                existing_idx = i
                break

        if existing_idx is not None:
            # Update existing entry
            self.data["structures"]["proteins"][existing_idx] = structure_entry
        else:
            # Add new entry
            self.data["structures"]["proteins"].append(structure_entry)
            self.data["structures"]["predicted"] = len(self.data["structures"]["proteins"])

    def get_structure(self, protein_id: str) -> Optional[dict]:
        """Get structure entry for a protein if it exists."""
        for s in self.data["structures"]["proteins"]:
            if s.get("protein_id") == protein_id:
                return s
        return None

    # ------------------------------------------------------------------ #
    # Figure tracking
    # ------------------------------------------------------------------ #

    def add_figure(
        self,
        path: str,
        figure_type: str,
        title: Optional[str] = None,
        legend: Optional[str] = None,
        center_protein: Optional[str] = None,
        finding_id: Optional[int] = None,
        **metadata: Any,
    ) -> None:
        """
        Record a generated figure with metadata.

        Args:
            path: Path to figure file (relative to dataset directory)
            figure_type: Type of figure ("neighborhood", "domain", "heatmap", etc.)
            title: Figure title
            legend: Figure legend/caption
            center_protein: Center protein ID (for neighborhood figures)
            finding_id: Associated finding ID
            **metadata: Additional metadata
        """
        # Make path relative to dataset directory if possible
        if self.db_path:
            try:
                path = str(Path(path).relative_to(self.db_path.parent))
            except ValueError:
                pass  # Keep absolute path if not relative

        figure_entry = {
            "path": path,
            "type": figure_type,
            "created_at": _now_iso(),
        }

        if title:
            figure_entry["title"] = title
        if legend:
            figure_entry["legend"] = legend
        if center_protein:
            figure_entry["center_protein"] = center_protein
        if finding_id is not None:
            figure_entry["finding_id"] = finding_id

        # Add any extra metadata
        for key, value in metadata.items():
            if value is not None:
                figure_entry[key] = value

        # Check for existing figure with same path
        existing_idx = None
        for i, f in enumerate(self.data["figures"]):
            if f.get("path") == path:
                existing_idx = i
                break

        if existing_idx is not None:
            self.data["figures"][existing_idx] = figure_entry
        else:
            self.data["figures"].append(figure_entry)

    def get_figure(self, path: str) -> Optional[dict]:
        """Get figure entry by path."""
        # Normalize path for comparison
        if self.db_path:
            try:
                path = str(Path(path).relative_to(self.db_path.parent))
            except ValueError:
                pass

        for f in self.data["figures"]:
            if f.get("path") == path:
                return f
        return None

    # ------------------------------------------------------------------ #
    # Findings tracking
    # ------------------------------------------------------------------ #

    def update_findings(self, findings_path: Optional[str] = None) -> None:
        """
        Update findings summary from findings.jsonl.

        Args:
            findings_path: Path to findings.jsonl (auto-detected if not provided)
        """
        if not self.db_path:
            return

        # Auto-detect findings file
        if not findings_path:
            for subdir in ["exploration", "survey"]:
                candidate = self.db_path.parent / subdir / "findings.jsonl"
                if candidate.exists():
                    findings_path = str(candidate)
                    break

        if not findings_path or not Path(findings_path).exists():
            return

        try:
            findings = []
            with open(findings_path) as f:
                for line in f:
                    line = line.strip()
                    if line:
                        findings.append(json.loads(line))

            # Update counts
            self.data["findings"]["count"] = len(findings)

            # Group by category
            by_category: dict[str, int] = {}
            high_priority: list[dict] = []
            proteins_with_findings: set[str] = set()

            for finding in findings:
                cat = finding.get("category", "Uncategorized")
                by_category[cat] = by_category.get(cat, 0) + 1

                if finding.get("priority") == "HIGH":
                    high_priority.append({
                        "id": finding.get("id"),
                        "title": finding.get("title"),
                    })

                # Track proteins mentioned in findings
                for gene in finding.get("genes", []):
                    proteins_with_findings.add(str(gene))

            self.data["findings"]["by_category"] = by_category
            self.data["findings"]["high_priority"] = high_priority
            self.data["exploration"]["coverage"]["proteins_with_findings"] = len(proteins_with_findings)

        except Exception as e:
            self.log_session("error", f"Failed to update findings: {e}")

    # ------------------------------------------------------------------ #
    # Exploration status
    # ------------------------------------------------------------------ #

    def set_exploration_status(self, status: str) -> None:
        """
        Set exploration status.

        Args:
            status: One of "not_started", "in_progress", "complete"
        """
        valid_statuses = {"not_started", "in_progress", "complete"}
        if status not in valid_statuses:
            status = "in_progress"
        self.data["exploration"]["status"] = status

    def add_phase_completed(self, phase: str) -> None:
        """Record a completed exploration phase."""
        if phase not in self.data["exploration"]["phases_completed"]:
            self.data["exploration"]["phases_completed"].append(phase)

    def add_predicate_queried(self, predicate: str) -> None:
        """Record a predicate that was queried."""
        if predicate not in self.data["exploration"]["coverage"]["predicates_queried"]:
            self.data["exploration"]["coverage"]["predicates_queried"].append(predicate)

    def add_contig_browsed(self, contig_id: str) -> None:
        """Record a contig that was browsed."""
        if contig_id not in self.data["exploration"]["coverage"]["contigs_browsed"]:
            self.data["exploration"]["coverage"]["contigs_browsed"].append(contig_id)

    # ------------------------------------------------------------------ #
    # Unexplored tracking
    # ------------------------------------------------------------------ #

    def get_unexplored(self, store: Optional["DuckDBStore"] = None) -> dict:
        """
        Calculate what hasn't been examined yet.

        Args:
            store: DuckDBStore instance for queries

        Returns:
            Dict with unexplored proteins, predicates, and suggested next steps
        """
        store = store or self._store
        unexplored = self.data.get("unexplored", {})

        if store:
            try:
                # Get total protein count
                result = store.execute("SELECT COUNT(*) FROM proteins")
                total_proteins = result[0][0] if result else 0

                proteins_with_findings = self.data["exploration"]["coverage"].get("proteins_with_findings", 0)
                unexplored["proteins_without_findings"] = total_proteins - proteins_with_findings

            except Exception:
                pass

        # Update predicates not queried (would need vocabulary list)
        # This is set manually or by exploration skill

        return unexplored

    def add_suggested_step(self, step: str) -> None:
        """Add a suggested next step."""
        if step not in self.data["unexplored"]["suggested_next_steps"]:
            self.data["unexplored"]["suggested_next_steps"].append(step)

    # ------------------------------------------------------------------ #
    # Reports
    # ------------------------------------------------------------------ #

    def add_report(
        self,
        path: str,
        pages: Optional[int] = None,
        figures_included: Optional[int] = None,
    ) -> None:
        """
        Record a generated report.

        Args:
            path: Path to report file
            pages: Number of pages
            figures_included: Number of figures included
        """
        if self.db_path:
            try:
                path = str(Path(path).relative_to(self.db_path.parent))
            except ValueError:
                pass

        report_entry = {
            "path": path,
            "generated": _now_iso(),
        }
        if pages is not None:
            report_entry["pages"] = pages
        if figures_included is not None:
            report_entry["figures_included"] = figures_included

        # Check for existing report with same path
        existing_idx = None
        for i, r in enumerate(self.data["reports"]):
            if r.get("path") == path:
                existing_idx = i
                break

        if existing_idx is not None:
            self.data["reports"][existing_idx] = report_entry
        else:
            self.data["reports"].append(report_entry)

    # ------------------------------------------------------------------ #
    # Status summary for session continuity
    # ------------------------------------------------------------------ #

    def get_status_summary(self) -> str:
        """
        Generate a summary of current analysis state for a new session.

        This provides context for an agent resuming work on this dataset.

        Returns:
            Formatted markdown summary of analysis state
        """
        d = self.data
        lines = [
            f"## Analysis State: {d['dataset']['name']}",
            "",
        ]

        # Status overview
        status = d.get("exploration", {}).get("status", "unknown")
        lines.append(f"**Status:** {status}")

        # Annotation summary
        annotations = d.get("annotations", {})
        if annotations:
            ann_summary = []
            for source, info in sorted(annotations.items()):
                if isinstance(info, dict) and "count" in info:
                    ann_summary.append(f"{source}: {info['count']}")
            if ann_summary:
                lines.append(f"**Annotations:** {', '.join(ann_summary)}")

        # Findings
        findings = d.get("findings", {})
        if findings.get("count", 0) > 0:
            lines.append(f"**Findings:** {findings['count']}")
            if findings.get("by_category"):
                cats = [f"{k}: {v}" for k, v in sorted(findings["by_category"].items())]
                lines.append(f"  Categories: {', '.join(cats[:5])}")

        # Structures
        structures = d.get("structures", {})
        if structures.get("predicted", 0) > 0:
            lines.append(f"**Structures predicted:** {structures['predicted']}")

        # Phases completed
        phases = d.get("exploration", {}).get("phases_completed", [])
        if phases:
            lines.append(f"**Phases completed:** {', '.join(phases)}")

        # Hypotheses from persistent registry
        if self.db_path:
            registry_path = self.db_path.parent / "exploration" / "hypotheses.json"
            if registry_path.exists():
                try:
                    from bennu.core.hypothesis_registry import HypothesisRegistry

                    registry = HypothesisRegistry(registry_path)
                    hypotheses = registry.list_all()
                    if hypotheses:
                        active = registry.list_active()
                        lines.append(f"**Hypotheses:** {len(hypotheses)} ({len(active)} active)")
                        lines.extend(["", "### Hypotheses"])
                        for h in hypotheses:
                            n_for = sum(1 for e in h.evidence if e.supports)
                            n_against = sum(1 for e in h.evidence if not e.supports)
                            lines.append(
                                f"- [{h.status.value}] {h.statement} "
                                f"(+{n_for}/-{n_against}, conf={h.confidence:.2f})"
                            )
                        # Update manifest counts
                        if "hypotheses" not in self.data:
                            self.data["hypotheses"] = {}
                        self.data["hypotheses"]["count"] = len(hypotheses)
                        self.data["hypotheses"]["active"] = len(active)
                except Exception:
                    pass  # Don't fail resume if registry is corrupted

        # Recent session activity
        session_log = d.get("session_log", [])
        if session_log:
            lines.extend(["", "### Recent Activity"])
            for entry in session_log[-5:]:
                ts = entry.get("timestamp", "")[:10]  # Just date
                action = entry.get("action", "")
                summary = entry.get("summary", "")[:80]
                lines.append(f"- [{ts}] **{action}**: {summary}")

        # Suggested next steps
        next_steps = d.get("unexplored", {}).get("suggested_next_steps", [])
        if next_steps:
            lines.extend(["", "### Suggested Next Steps"])
            for step in next_steps[:5]:
                lines.append(f"- {step}")

        # High priority findings
        high_priority = d.get("findings", {}).get("high_priority", [])
        if high_priority:
            lines.extend(["", "### High Priority Findings"])
            for hp in high_priority[:5]:
                lines.append(f"- [{hp.get('id')}] {hp.get('title', '')[:60]}")

        return "\n".join(lines)

    # ------------------------------------------------------------------ #
    # Persistence
    # ------------------------------------------------------------------ #

    def save(self) -> None:
        """Write manifest to disk with pretty-printing for git-friendliness."""
        if not self.path:
            return

        # Ensure directory exists
        self.path.parent.mkdir(parents=True, exist_ok=True)

        # Write with indent for readability
        self.path.write_text(json.dumps(self.data, indent=2, default=str))

    def __enter__(self) -> "AnalysisManifest":
        """Context manager entry."""
        return self

    def __exit__(self, *args) -> None:
        """Context manager exit - auto-save."""
        self.save()


__all__ = ["AnalysisManifest"]
