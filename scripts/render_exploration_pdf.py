#!/usr/bin/env python3
"""Render a normalized exploration findings PDF report."""

from __future__ import annotations

import argparse

from sharur.reports.findings_pdf import generate_phase_report


def main() -> int:
    parser = argparse.ArgumentParser(description="Render an exploration findings PDF.")
    parser.add_argument("--dataset", required=True, help="Path to the dataset directory.")
    parser.add_argument("--output", help="Optional output PDF path.")
    args = parser.parse_args()

    output = generate_phase_report(
        args.dataset,
        "exploration",
        output_path=args.output,
    )
    print(f"Rendered exploration PDF: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
