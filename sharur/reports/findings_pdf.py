"""
PDF report generation from normalized findings records.
"""

from __future__ import annotations

import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from sharur.core.analysis_record_compat import (
    FindingNormalizationResult,
    build_findings_summary,
    load_findings_file,
)
from sharur.reports.template import SharurReport


def _infer_theme(dataset_name: str) -> str:
    name_lower = dataset_name.lower()
    if any(x in name_lower for x in ["archae", "thermo", "methano"]):
        return "archaeal"
    if any(x in name_lower for x in ["virus", "phage", "ncldv"]):
        return "viral"
    if any(x in name_lower for x in ["bacter", "pseudo", "strepto"]):
        return "bacterial"
    return "default"


def _phase_title(phase: str) -> str:
    return phase.replace("_", " ").title()


def _load_manifest(dataset_dir: Path) -> dict[str, Any]:
    manifest_path = dataset_dir / "manifest.json"
    if not manifest_path.exists():
        return {}
    with manifest_path.open() as handle:
        return json.load(handle)


def _format_compact_evidence(evidence: Any) -> str | None:
    if not evidence:
        return None
    if isinstance(evidence, dict):
        parts = []
        for key, value in evidence.items():
            if isinstance(value, (dict, list)):
                continue
            parts.append(f"{key}: {value}")
            if len(parts) == 6:
                break
        if parts:
            return ", ".join(parts)
        return json.dumps(evidence, default=str)[:300]
    return str(evidence)[:300]


def _issue_frequency(records: list[FindingNormalizationResult]) -> list[tuple[str, int]]:
    counts: Counter[str] = Counter()
    for record in records:
        counts.update(record.issues)
    return counts.most_common()


def _group_by_category(records: list[FindingNormalizationResult]) -> dict[str, list[FindingNormalizationResult]]:
    grouped: dict[str, list[FindingNormalizationResult]] = defaultdict(list)
    for record in records:
        grouped[record.finding.get("category", "general")].append(record)
    return dict(grouped)


def generate_phase_report(
    dataset_dir: str | Path,
    phase: str,
    *,
    output_path: str | Path | None = None,
    theme: str | None = None,
) -> str:
    """Generate a PDF report for one findings phase."""
    dataset_dir = Path(dataset_dir)
    phase_path = dataset_dir / phase / "findings.jsonl"
    records = load_findings_file(phase_path)
    if not records:
        raise FileNotFoundError(f"No findings found at {phase_path}")

    manifest = _load_manifest(dataset_dir)
    dataset_name = manifest.get("dataset", {}).get("name", dataset_dir.name)
    findings_summary, proteins_with_findings = build_findings_summary(records)

    if output_path is None:
        output_path = dataset_dir / f"{dataset_name}_{phase}_report.pdf"
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    report_theme = theme or _infer_theme(dataset_name)
    pdf = SharurReport(
        theme=report_theme,
        title=f"{dataset_name} {_phase_title(phase)} Report",
    )
    pdf.create_title_page(
        title=dataset_name.replace("_", " ").title(),
        subtitle=f"{_phase_title(phase)} Findings",
        stats={
            "Findings": findings_summary["count"],
            "Key Findings": len(findings_summary["key_findings"]),
            "Proteins": proteins_with_findings,
            "Validation Issues": len(findings_summary["validation_issues"]),
        },
    )
    toc_page = pdf.create_toc()

    # Executive Summary
    pdf.add_page()
    pdf.chapter_title(1, "Executive Summary")
    pdf.body_text(
        f"This report summarizes {findings_summary['count']} normalized {_phase_title(phase).lower()} findings "
        f"for {dataset_name}. The current archive references {proteins_with_findings} canonical protein IDs."
    )
    if findings_summary["validation_issues"]:
        pdf.add_evidence_box(
            f"{len(findings_summary['validation_issues'])} findings still need manual cleanup. "
            f"Most common issues: {', '.join(f'{issue} ({count})' for issue, count in _issue_frequency(records)[:3])}"
        )
    if findings_summary["key_findings"]:
        pdf.section_title("Key Findings")
        for finding in findings_summary["key_findings"][:5]:
            pdf.add_bullet(f"[{finding['id']}] {finding['title']}")

    # Category Summary
    pdf.add_page()
    pdf.chapter_title(2, "Category Summary")
    table_data = [
        [category, str(count)]
        for category, count in sorted(
            findings_summary["by_category"].items(),
            key=lambda item: (-item[1], item[0]),
        )
    ]
    if table_data:
        pdf.add_table(["Category", "Count"], table_data, [110, 40])

    # Validation summary
    if findings_summary["validation_issues"]:
        pdf.add_page()
        pdf.chapter_title(3, "Validation Summary")
        for issue, count in _issue_frequency(records)[:8]:
            pdf.add_bullet(f"{issue}: {count}")

    grouped = _group_by_category(records)
    chapter_num = 4 if findings_summary["validation_issues"] else 3
    for category in sorted(grouped):
        pdf.add_page()
        pdf.chapter_title(chapter_num, category.replace("_", " ").title())
        chapter_num += 1
        for record in grouped[category]:
            finding = record.finding
            pdf.subsection_title(f"[{finding.get('id')}] {finding.get('title', 'Untitled')}")
            meta_bits = [f"Phase: {finding.get('phase', 'unknown')}"]
            if finding.get("n_genomes") is not None:
                meta_bits.append(f"Genomes: {finding['n_genomes']}")
            if finding.get("protein_ids"):
                meta_bits.append(f"Protein IDs: {len(finding['protein_ids'])}")
            pdf.body_text(" | ".join(meta_bits))
            pdf.body_text(finding.get("description", ""))

            evidence_text = _format_compact_evidence(finding.get("evidence"))
            if evidence_text:
                pdf.add_evidence_box(f"Evidence: {evidence_text}")

            verification = finding.get("verification", [])
            if verification:
                pdf.section_title("Verification")
                for entry in verification[:3]:
                    message = entry.get("claim", "Verification claim")
                    if entry.get("expected") is not None:
                        message += f" (expected {entry['expected']})"
                    if entry.get("query"):
                        message += f": {str(entry['query'])[:140]}"
                    pdf.add_bullet(message)

            if record.issues:
                pdf.add_evidence_box(
                    "Validation issues: " + "; ".join(record.issues)
                )

            if finding.get("related_findings"):
                pdf.body_text(
                    "Related findings: " + ", ".join(finding["related_findings"])
                )
            pdf.add_separator()

    figures_dir = dataset_dir / phase / "figures"
    figure_files = sorted(figures_dir.glob("*.png")) if figures_dir.exists() else []
    manifest_figures = {
        figure.get("path"): figure
        for figure in manifest.get("figures", [])
        if isinstance(figure, dict)
    }
    if figure_files:
        pdf.add_page()
        pdf.chapter_title(chapter_num, "Figures")
        chapter_num += 1
        for figure_path in figure_files[:12]:
            rel_path = str(figure_path.relative_to(dataset_dir))
            figure_meta = manifest_figures.get(rel_path, {})
            pdf.add_image_with_legend(
                str(figure_path),
                figure_meta.get("title", figure_path.stem.replace("_", " ").title()),
                figure_meta.get("legend", "Generated figure associated with this findings phase."),
            )

    pdf.add_page()
    pdf.chapter_title(chapter_num, "Methods")
    pdf.body_text(
        "Report generated from normalized findings records using the canonical findings compatibility layer. "
        "Legacy IDs are preserved on read, verification blocks are rendered directly when present, and "
        "remaining schema issues are surfaced without inventing missing scientific content."
    )

    pdf.fill_toc(toc_page)
    pdf.output(str(output_path))
    return str(output_path)


__all__ = ["generate_phase_report"]
