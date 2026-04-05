"""
Sharur report generation module.

Provides templated report generators that consume analysis manifests
to produce consistent, publication-ready PDF reports.
"""

from sharur.reports.template import (
    SharurReport,
    generate_report_from_manifest,
    clean_text,
)
from sharur.reports.findings_pdf import generate_phase_report

__all__ = [
    "SharurReport",
    "generate_report_from_manifest",
    "generate_phase_report",
    "clean_text",
]
