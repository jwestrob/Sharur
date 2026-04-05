#!/usr/bin/env python3
"""Compile structured findings + report manifest into a comprehensive Markdown report.

Reads findings.jsonl from survey/ and exploration/ directories, plus an optional
reports/report_manifest.json that organizes findings into thematic sections.

Usage:
    python scripts/compile_comprehensive_report.py --dataset data/omni_production/
    python scripts/compile_comprehensive_report.py --dataset data/omni_production/ --output COMPREHENSIVE_REPORT.md
"""

import argparse
import json
from collections import Counter
from datetime import datetime
from pathlib import Path

from sharur.core.analysis_record_compat import (
    FindingNormalizationResult,
    load_dataset_findings,
)


def load_findings(dataset_dir: Path) -> list[FindingNormalizationResult]:
    """Load and normalize canonical findings across survey and exploration."""
    return load_dataset_findings(dataset_dir)


def load_manifest(dataset_dir: Path) -> dict | None:
    """Load report manifest if it exists."""
    manifest_path = dataset_dir / "reports" / "report_manifest.json"
    if not manifest_path.exists():
        return None
    with open(manifest_path) as f:
        return json.load(f)


def _truncate(text: str, max_length: int = 200) -> str:
    if len(text) <= max_length:
        return text
    return f"{text[: max_length - 3]}..."


def render_finding(record: FindingNormalizationResult) -> list[str]:
    """Render a single finding as Markdown lines."""
    finding = record.finding
    lines = []
    fid = finding.get("id", "???")
    title = finding.get("title", "Untitled")
    lines.append(f"### [{fid}] {title}")

    category = finding.get("category", "uncategorized")
    phase = finding.get("phase", "unknown")
    n_genomes = finding.get("n_genomes")

    meta_parts = [f"**Category:** {category}", f"**Phase:** {phase}"]
    if n_genomes is not None:
        meta_parts.append(f"**Genomes:** {n_genomes}")
    lines.append(" | ".join(meta_parts))
    lines.append("")

    if finding.get("description"):
        lines.append(finding["description"])
        lines.append("")

    if finding.get("evidence"):
        evidence = finding["evidence"]
        if isinstance(evidence, dict):
            lines.append("**Evidence:**")
            lines.append("```json")
            lines.append(json.dumps(evidence, indent=2, sort_keys=True, default=str))
            lines.append("```")
        else:
            lines.append(f"**Evidence:** {evidence}")
        lines.append("")

    prov = finding.get("provenance", {})
    if prov.get("query_used"):
        lines.append(f"**Provenance:** `{_truncate(str(prov['query_used']))}`")
        lines.append("")

    verification = finding.get("verification", [])
    if verification:
        lines.append("**Verification:**")
        for entry in verification:
            claim = entry.get("claim", "Unlabeled verification claim")
            expected = entry.get("expected")
            query = entry.get("query")
            detail = f"- {claim}"
            if expected is not None:
                detail += f" (expected `{expected}`)"
            if query:
                detail += f": `{_truncate(str(query))}`"
            lines.append(detail)
        lines.append("")

    if finding.get("figures"):
        for fig in finding["figures"]:
            lines.append(f"![{fid}]({fig})")
        lines.append("")

    if finding.get("related_findings"):
        related = ", ".join(finding["related_findings"])
        lines.append(f"**Related:** {related}")
        lines.append("")

    if record.issues:
        lines.append(f"**Validation issues:** {'; '.join(record.issues)}")
        lines.append("")

    return lines


def compile_report(dataset_dir: Path) -> str:
    """Compile the full comprehensive report."""
    findings = load_findings(dataset_dir)
    manifest = load_manifest(dataset_dir)

    if not findings:
        return "# Comprehensive Report\n\nNo findings found.\n"

    # Index findings by ID
    by_id = {}
    for record in findings:
        fid = record.finding.get("id")
        if fid:
            by_id[fid] = record

    # Count stats
    phases = Counter(record.finding.get("phase", "unknown") for record in findings)
    categories = Counter(
        record.finding.get("category", "uncategorized") for record in findings
    )
    n_figures = sum(len(record.finding.get("figures", [])) for record in findings)

    dataset_name = dataset_dir.name
    lines = [
        f"# Comprehensive Report: {dataset_name}",
        "",
        f"*Auto-generated from {len(findings)} findings across {len(phases)} phases. "
        f"{n_figures} figures referenced.*",
        f"*Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}*",
        "",
        "---",
        "",
    ]

    # Track which findings have been rendered (to find orphans)
    rendered_ids = set()

    if manifest and manifest.get("sections"):
        # Manifest-driven rendering
        for section in manifest["sections"]:
            section_title = section.get("title", "Untitled Section")
            lines.append(f"## {section_title}")
            lines.append("")

            if section.get("narrative"):
                lines.append(f"*Narrative: {section['narrative']}*")
                lines.append("")

            finding_ids = section.get("finding_ids", [])
            if not finding_ids:
                lines.append("*No findings assigned to this section.*")
                lines.append("")
                continue

            for fid in finding_ids:
                if fid in by_id:
                    lines.extend(render_finding(by_id[fid]))
                    rendered_ids.add(fid)
                else:
                    lines.append(f"### [{fid}] *(finding not found)*")
                    lines.append("")

            if section.get("figures"):
                lines.append("**Section figures:**")
                for fig in section["figures"]:
                    lines.append(f"- {fig}")
                lines.append("")

        # Excluded findings
        excluded = manifest.get("excluded_findings", [])
        if excluded:
            lines.append("## Excluded Findings")
            lines.append("")
            lines.append("| ID | Reason |")
            lines.append("|----|--------|")
            for ex in excluded:
                lines.append(f"| {ex.get('id', '?')} | {ex.get('reason', 'No reason given')} |")
                rendered_ids.add(ex.get("id"))
            lines.append("")

        # Orphan findings (in findings.jsonl but not in manifest)
        orphan_ids = [
            record.finding.get("id") for record in findings
            if record.finding.get("id") and record.finding["id"] not in rendered_ids
        ]
        if orphan_ids:
            lines.append("## Unassigned Findings")
            lines.append("")
            lines.append(
                f"*{len(orphan_ids)} findings not referenced by any manifest section.*"
            )
            lines.append("")
            for fid in orphan_ids:
                if fid in by_id:
                    lines.extend(render_finding(by_id[fid]))

    else:
        # No manifest — group by category
        lines.append("*No report manifest found. Findings grouped by category.*")
        lines.append("")

        findings_by_cat: dict[str, list[FindingNormalizationResult]] = {}
        for record in findings:
            cat = record.finding.get("category", "uncategorized")
            findings_by_cat.setdefault(cat, []).append(record)

        for cat in sorted(findings_by_cat.keys()):
            cat_findings = findings_by_cat[cat]
            lines.append(
                f"## {cat.replace('_', ' ').title()} ({len(cat_findings)} findings)"
            )
            lines.append("")
            for record in cat_findings:
                lines.extend(render_finding(record))

    # Summary statistics
    lines.append("---")
    lines.append("")
    lines.append("## Summary Statistics")
    lines.append("")
    lines.append("### Findings by Phase")
    lines.append("")
    lines.append("| Phase | Count |")
    lines.append("|-------|-------|")
    for phase, count in phases.most_common():
        lines.append(f"| {phase} | {count} |")
    lines.append("")

    lines.append("### Findings by Category")
    lines.append("")
    lines.append("| Category | Count |")
    lines.append("|----------|-------|")
    for cat, count in categories.most_common():
        lines.append(f"| {cat} | {count} |")
    lines.append("")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Compile findings into a comprehensive Markdown report."
    )
    parser.add_argument(
        "--dataset",
        required=True,
        help="Path to the dataset directory (e.g., data/omni_production/)",
    )
    parser.add_argument(
        "--output",
        default="COMPREHENSIVE_REPORT.md",
        help="Output filename (default: COMPREHENSIVE_REPORT.md)",
    )
    args = parser.parse_args()

    dataset_dir = Path(args.dataset)
    if not dataset_dir.exists():
        print(f"Error: dataset directory not found: {dataset_dir}")
        return 1

    report = compile_report(dataset_dir)
    output_path = dataset_dir / args.output
    output_path.write_text(report)

    # Count for summary
    findings = load_findings(dataset_dir)
    manifest = load_manifest(dataset_dir)
    n_sections = len(manifest["sections"]) if manifest and manifest.get("sections") else 0

    print(f"Compiled {len(findings)} findings into {output_path}")
    if manifest:
        print(f"  Manifest: {n_sections} sections")
    else:
        print("  No manifest found; findings grouped by category")
    print(f"  Output: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
