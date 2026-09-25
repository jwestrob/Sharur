#!/usr/bin/env python3
"""Deprecated wrapper for the normalized exploration PDF generator."""

from __future__ import annotations

import argparse

from sharur.reports.findings_pdf import generate_phase_report


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Deprecated: use scripts/render_exploration_pdf.py instead.",
    )
    parser.add_argument(
        "--dataset",
        default="data/DATASET",
        help="Path to the dataset directory.",
    )
    parser.add_argument("--output", help="Optional output PDF path.")
    args = parser.parse_args()

    print("Deprecated: use scripts/render_exploration_pdf.py instead.")
    output = generate_phase_report(
        args.dataset,
        "exploration",
        output_path=args.output,
    )
    print(f"Rendered exploration PDF: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
