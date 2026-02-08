"""
Bennu report generation module.

Provides templated report generators that consume analysis manifests
to produce consistent, publication-ready PDF reports.
"""

from bennu.reports.template import (
    BennuReport,
    generate_report_from_manifest,
    clean_text,
)

__all__ = [
    "BennuReport",
    "generate_report_from_manifest",
    "clean_text",
]
