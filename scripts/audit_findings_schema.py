#!/usr/bin/env python3
"""Audit and safely migrate dataset findings archives."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from sharur.core.analysis_record_audit import (
    audit_dataset_findings,
    render_findings_audit_markdown,
)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Audit findings schema drift and optionally apply safe rewrites.",
    )
    parser.add_argument(
        "--dataset",
        required=True,
        help="Path to the dataset directory.",
    )
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Apply safe rewrites in-place.",
    )
    parser.add_argument(
        "--convert-legacy-evidence",
        action="store_true",
        help=(
            "Convert legacy free-text evidence strings into object wrappers "
            "that preserve the original text and clause-level statements."
        ),
    )
    parser.add_argument(
        "--report",
        help="Optional path to write a Markdown audit report.",
    )
    parser.add_argument(
        "--json",
        dest="json_output",
        help="Optional path to write the raw audit summary as JSON.",
    )
    args = parser.parse_args()

    dataset_dir = Path(args.dataset)
    if not dataset_dir.exists():
        raise SystemExit(f"Dataset not found: {dataset_dir}")

    summary = audit_dataset_findings(
        dataset_dir,
        apply_changes=args.apply,
        convert_legacy_evidence=args.convert_legacy_evidence,
    )
    report_text = render_findings_audit_markdown(summary)
    print(report_text)

    if args.report:
        report_path = Path(args.report)
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(report_text)
        print(f"\nWrote Markdown report to {report_path}")

    if args.json_output:
        json_path = Path(args.json_output)
        json_path.parent.mkdir(parents=True, exist_ok=True)
        json_path.write_text(json.dumps(summary, indent=2))
        print(f"Wrote JSON summary to {json_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
