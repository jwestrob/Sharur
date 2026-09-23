#!/usr/bin/env python3
"""Validate the provenance chain from findings → manuscript claims → citations.

Checks:
1. Every source_findings ID in MANUSCRIPT_CLAIMS.jsonl exists in some findings.jsonl
2. Every source_citations ID exists in literature_citations.jsonl
3. Every figure path referenced in claims/manifest exists on disk
4. Quantitative claims in MANUSCRIPT.md without corresponding claim entries
5. Orphan findings (in findings.jsonl but unreferenced by any claim or manifest)

Usage:
    python scripts/validate_provenance.py --dataset data/DATASET/
"""

import argparse
import json
import re
from pathlib import Path


def load_jsonl(path: Path) -> list[dict]:
    """Load a JSONL file, returning empty list if missing."""
    if not path.exists():
        return []
    entries = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                entries.append(json.loads(line))
    return entries


def load_json(path: Path) -> dict | None:
    if not path.exists():
        return None
    with open(path) as f:
        return json.load(f)


def find_quantitative_claims(manuscript_path: Path) -> list[dict]:
    """Find lines in MANUSCRIPT.md with numbers/percentages that might be claims."""
    if not manuscript_path.exists():
        return []
    claims = []
    with open(manuscript_path) as f:
        for i, line in enumerate(f, 1):
            line = line.strip()
            # Skip headings, figure captions, references, empty lines
            if not line or line.startswith("#") or line.startswith("!") or line.startswith("**Figure"):
                continue
            # Match patterns like "74.6%", "1,831 genomes", "2,034 proteins"
            if re.search(r'\d{2,}[,.]?\d*\s*(%|genomes|proteins|aa\b|amino acid)', line):
                claims.append({"line": i, "text": line[:120]})
    return claims


def validate_claims_schema(claims: list[dict]) -> list[dict]:
    """Check MANUSCRIPT_CLAIMS.jsonl entries for schema drift.

    Canonical fields:
      Required: claim_id, section, claim_text, source_findings (list), claim_type
      Optional: source_citations (list), figures (list), review_status
    """
    # Fields required for provenance chain integrity
    REQUIRED_FIELDS = {"claim_id", "claim_text", "source_findings", "claim_type"}
    # Fields recommended for full traceability but not chain-breaking
    RECOMMENDED_FIELDS = {"section"}
    # Common misspellings / drift from the canonical schema
    FIELD_CORRECTIONS = {
        "id": "claim_id",
        "claim": "claim_text",
        "text": "claim_text",
        "type": "claim_type",
        "source_finding": "source_findings",
        "finding": "source_findings",
        "findings": "source_findings",
        "citation": "source_citations",
        "source_citation": "source_citations",
    }
    issues = []
    for i, claim in enumerate(claims):
        # Check for drifted field names
        for wrong, correct in FIELD_CORRECTIONS.items():
            if wrong in claim and correct not in claim:
                issues.append({
                    "type": "schema_drift",
                    "entry": i + 1,
                    "detail": f"Entry {i+1}: uses '{wrong}' instead of '{correct}'",
                })
        # Check for missing required fields
        for field in REQUIRED_FIELDS:
            if field not in claim:
                # Only flag if there's no drifted equivalent
                has_equivalent = any(
                    wrong in claim for wrong, correct in FIELD_CORRECTIONS.items()
                    if correct == field
                )
                if not has_equivalent:
                    issues.append({
                        "type": "schema_missing_field",
                        "entry": i + 1,
                        "detail": f"Entry {i+1}: missing required field '{field}'",
                    })
        # Check for missing recommended fields (warn, not error)
        for field in RECOMMENDED_FIELDS:
            if field not in claim:
                issues.append({
                    "type": "schema_recommended",
                    "entry": i + 1,
                    "detail": f"Entry {i+1}: missing recommended field '{field}'",
                })
        # Check that source_findings is a list (not a string)
        sf = claim.get("source_findings", claim.get("source_finding"))
        if sf is not None and isinstance(sf, str):
            issues.append({
                "type": "schema_type_error",
                "entry": i + 1,
                "detail": f"Entry {i+1}: 'source_findings' must be a list, got string '{sf}'",
            })
    return issues


def validate_manifest_schema(manifest: dict) -> list[dict]:
    """Check report_manifest.json sections for schema drift.

    Canonical section fields: title, finding_ids (list)
    Optional: id, figures, narrative/narratives
    """
    FIELD_CORRECTIONS = {
        "findings": "finding_ids",
        "finding_list": "finding_ids",
    }
    issues = []
    for i, section in enumerate(manifest.get("sections", [])):
        for wrong, correct in FIELD_CORRECTIONS.items():
            if wrong in section and correct not in section:
                issues.append({
                    "type": "schema_drift",
                    "section": section.get("title", f"section {i+1}"),
                    "detail": f"Manifest section '{section.get('title', i+1)}': uses '{wrong}' instead of '{correct}'",
                })
        if "finding_ids" not in section and "findings" not in section:
            issues.append({
                "type": "schema_missing_field",
                "section": section.get("title", f"section {i+1}"),
                "detail": f"Manifest section '{section.get('title', i+1)}': missing 'finding_ids'",
            })
    return issues


def validate(dataset_dir: Path) -> str:
    """Run all provenance checks and return audit report."""
    # Load all data
    survey_findings = load_jsonl(dataset_dir / "survey" / "findings.jsonl")
    explore_findings = load_jsonl(dataset_dir / "exploration" / "findings.jsonl")
    all_findings = survey_findings + explore_findings

    claims = load_jsonl(dataset_dir / "MANUSCRIPT_CLAIMS.jsonl")
    citations = load_jsonl(dataset_dir / "literature_citations.jsonl")
    manifest = load_json(dataset_dir / "reports" / "report_manifest.json")
    manuscript_path = dataset_dir / "MANUSCRIPT.md"

    # === SCHEMA VALIDATION ===
    schema_issues = []
    if claims:
        schema_issues.extend(validate_claims_schema(claims))
    if manifest:
        schema_issues.extend(validate_manifest_schema(manifest))

    # Normalize claims: resolve drifted field names for downstream checks
    for claim in claims:
        if "id" in claim and "claim_id" not in claim:
            claim["claim_id"] = claim["id"]
        if "claim" in claim and "claim_text" not in claim:
            claim["claim_text"] = claim["claim"]
        if "type" in claim and "claim_type" not in claim:
            claim["claim_type"] = claim["type"]
        if "source_finding" in claim and "source_findings" not in claim:
            sf = claim["source_finding"]
            claim["source_findings"] = [sf] if isinstance(sf, str) else sf

    # Normalize manifest: resolve drifted field names for downstream checks
    if manifest:
        for section in manifest.get("sections", []):
            if "findings" in section and "finding_ids" not in section:
                section["finding_ids"] = section["findings"]

    # Build indexes
    finding_ids = {f.get("id") for f in all_findings if f.get("id")}
    citation_ids = {c.get("citation_id") for c in citations if c.get("citation_id")}

    # IDs referenced by claims
    claimed_finding_ids = set()
    claimed_citation_ids = set()
    claimed_figure_paths = set()
    for c in claims:
        for fid in c.get("source_findings", []):
            claimed_finding_ids.add(fid)
        for cid in c.get("source_citations", []):
            claimed_citation_ids.add(cid)
        for fig in c.get("figures", []):
            claimed_figure_paths.add(fig)

    # IDs referenced by manifest
    manifest_finding_ids = set()
    manifest_figure_paths = set()
    if manifest:
        for section in manifest.get("sections", []):
            for fid in section.get("finding_ids", []):
                manifest_finding_ids.add(fid)
            for fig in section.get("figures", []):
                manifest_figure_paths.add(fig)
        for ex in manifest.get("excluded_findings", []):
            manifest_finding_ids.add(ex.get("id"))

    all_figure_paths = claimed_figure_paths | manifest_figure_paths

    # === CHECKS ===

    issues = []

    # 1. Missing findings (referenced by claims but not in findings.jsonl)
    missing_findings = claimed_finding_ids - finding_ids
    for fid in sorted(missing_findings):
        # Find which claim references it
        for c in claims:
            if fid in c.get("source_findings", []):
                issues.append({
                    "type": "missing_finding",
                    "claim_id": c.get("claim_id"),
                    "section": c.get("section"),
                    "finding_id": fid,
                    "detail": f"Claim {c.get('claim_id')} references finding '{fid}' not found in any findings.jsonl",
                })
                break

    # 2. Missing citations (referenced by claims but not in literature_citations.jsonl)
    missing_citations = claimed_citation_ids - citation_ids
    for cid in sorted(missing_citations):
        for c in claims:
            if cid in c.get("source_citations", []):
                issues.append({
                    "type": "missing_citation",
                    "claim_id": c.get("claim_id"),
                    "section": c.get("section"),
                    "citation_id": cid,
                    "detail": f"Claim {c.get('claim_id')} references citation '{cid}' not found in literature_citations.jsonl",
                })
                break

    # 3. Missing figures (referenced but not on disk)
    for fig in sorted(all_figure_paths):
        fig_path = dataset_dir / fig
        if not fig_path.exists():
            issues.append({
                "type": "missing_figure",
                "path": fig,
                "detail": f"Figure '{fig}' referenced but not found on disk",
            })

    # 4. Orphan findings (in findings.jsonl but not in any claim or manifest)
    all_referenced = claimed_finding_ids | manifest_finding_ids
    orphan_ids = finding_ids - all_referenced
    orphans = []
    for f in all_findings:
        fid = f.get("id")
        if fid and fid in orphan_ids:
            orphans.append({"id": fid, "title": f.get("title", ""), "phase": f.get("phase", "")})

    # 5. Quantitative claims in manuscript without claim entries
    quant_claims = find_quantitative_claims(manuscript_path)
    claim_texts = {c.get("claim_text", "") for c in claims}
    uncovered_quant = []
    for qc in quant_claims:
        # Simple heuristic: check if any claim_text overlaps with this line
        covered = any(
            ct and (ct[:40] in qc["text"] or qc["text"][:40] in ct)
            for ct in claim_texts
        )
        if not covered:
            uncovered_quant.append(qc)

    # === BUILD REPORT ===

    lines = ["# Provenance Audit", ""]

    # Summary — count claims with at least one valid finding reference
    n_traceable = sum(
        1 for c in claims
        if any(fid in finding_ids for fid in c.get("source_findings", []))
    )
    n_total_claims = len(claims)
    pct = round(100 * n_traceable / n_total_claims, 1) if n_total_claims else 0

    lines.append("## Summary")
    lines.append("")
    lines.append(f"- **Findings loaded:** {len(all_findings)} ({len(survey_findings)} survey, {len(explore_findings)} exploration)")
    lines.append(f"- **Manuscript claims:** {n_total_claims}")
    lines.append(f"- **Traceable to findings:** {n_traceable} ({pct}%)")
    lines.append(f"- **Missing finding references:** {len(missing_findings)}")
    lines.append(f"- **Missing citation references:** {len(missing_citations)}")
    lines.append(f"- **Missing figures:** {sum(1 for i in issues if i['type'] == 'missing_figure')}")
    lines.append(f"- **Orphan findings:** {len(orphans)}")
    lines.append(f"- **Literature citations loaded:** {len(citations)}")
    if manifest:
        lines.append(f"- **Report manifest:** {len(manifest.get('sections', []))} sections")
    else:
        lines.append("- **Report manifest:** not found")
    lines.append("")

    # Schema drift warnings — separate hard errors from soft recommendations
    schema_errors = [i for i in schema_issues if i["type"] != "schema_recommended"]
    schema_recs = [i for i in schema_issues if i["type"] == "schema_recommended"]

    if schema_errors:
        lines.append("## Schema Errors")
        lines.append("")
        lines.append(f"*{len(schema_errors)} schema error(s). Fix these — they break the provenance chain.*")
        lines.append("")
        lines.append("| Type | Detail |")
        lines.append("|------|--------|")
        for issue in schema_errors:
            lines.append(f"| {issue['type']} | {issue['detail']} |")
        lines.append("")

    if schema_recs:
        # Deduplicate — e.g. "missing 'section'" appears N times, just summarize
        rec_fields = {}
        for r in schema_recs:
            field = r["detail"].split("'")[1] if "'" in r["detail"] else r["detail"]
            rec_fields[field] = rec_fields.get(field, 0) + 1
        lines.append("## Schema Recommendations")
        lines.append("")
        for field, count in rec_fields.items():
            lines.append(f"- {count} claim(s) missing recommended field `{field}`")
        lines.append("")

    # Missing provenance
    if issues:
        lines.append("## Issues")
        lines.append("")
        lines.append("| Type | Detail |")
        lines.append("|------|--------|")
        for issue in issues:
            lines.append(f"| {issue['type']} | {issue['detail']} |")
        lines.append("")

    # Orphan findings
    if orphans:
        lines.append("## Orphan Findings")
        lines.append("")
        lines.append(f"*{len(orphans)} findings not referenced by any manuscript claim or report manifest section.*")
        lines.append("")
        lines.append("| ID | Title | Phase |")
        lines.append("|----|-------|-------|")
        for o in sorted(orphans, key=lambda x: x["id"]):
            lines.append(f"| {o['id']} | {o['title'][:80]} | {o['phase']} |")
        lines.append("")

    # Uncovered quantitative claims
    if uncovered_quant and claims:
        lines.append("## Potentially Uncovered Quantitative Claims")
        lines.append("")
        lines.append(f"*{len(uncovered_quant)} lines in MANUSCRIPT.md contain quantitative data "
                      "but may not have corresponding MANUSCRIPT_CLAIMS.jsonl entries (heuristic check).*")
        lines.append("")
        lines.append("| Line | Text |")
        lines.append("|------|------|")
        for qc in uncovered_quant[:30]:  # Cap at 30
            text = qc["text"].replace("|", "\\|")
            lines.append(f"| {qc['line']} | {text} |")
        if len(uncovered_quant) > 30:
            lines.append(f"| ... | *{len(uncovered_quant) - 30} more* |")
        lines.append("")

    # File inventory
    lines.append("## File Inventory")
    lines.append("")
    files = {
        "survey/findings.jsonl": (dataset_dir / "survey" / "findings.jsonl").exists(),
        "exploration/findings.jsonl": (dataset_dir / "exploration" / "findings.jsonl").exists(),
        "MANUSCRIPT.md": manuscript_path.exists(),
        "MANUSCRIPT_CLAIMS.jsonl": (dataset_dir / "MANUSCRIPT_CLAIMS.jsonl").exists(),
        "literature_citations.jsonl": (dataset_dir / "literature_citations.jsonl").exists(),
        "reports/report_manifest.json": (dataset_dir / "reports" / "report_manifest.json").exists(),
    }
    for fname, exists in files.items():
        status = "found" if exists else "**MISSING**"
        lines.append(f"- `{fname}`: {status}")
    lines.append("")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Validate provenance chain from findings to manuscript."
    )
    parser.add_argument(
        "--dataset",
        required=True,
        help="Path to the dataset directory",
    )
    parser.add_argument(
        "--output",
        default="PROVENANCE_AUDIT.md",
        help="Output filename (default: PROVENANCE_AUDIT.md)",
    )
    args = parser.parse_args()

    dataset_dir = Path(args.dataset)
    if not dataset_dir.exists():
        print(f"Error: dataset directory not found: {dataset_dir}")
        return 1

    report = validate(dataset_dir)
    output_path = dataset_dir / args.output
    output_path.write_text(report)

    print(f"Provenance audit written to {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
