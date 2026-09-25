#!/usr/bin/env python3
"""Regenerate docs/reference/cli.md from the installed console scripts.

Typer apps are rendered with ``typer ... utils docs``; argparse entry points
contribute their ``--help`` text.

Usage:
    python scripts/gen_cli_docs.py
"""

from __future__ import annotations

import subprocess
from pathlib import Path


OUT = Path(__file__).resolve().parents[1] / "docs/reference/cli.md"
TYPER_APPS = [
    ("sharur", "sharur.cli"),
    ("sharur-ingest", "sharur.ingest_cli"),
    ("sharur-atlas", "sharur.atlas_cli"),
    ("sharur-review", "sharur.review_cli"),
    ("sharur-worker", "sharur.workers_cli"),
]
ARGPARSE_APPS = ["sharur-ops", "sharur-query"]


def typer_docs(name: str, module: str) -> str:
    import importlib  # noqa: PLC0415

    import typer  # noqa: PLC0415
    from typer.cli import get_docs_for_click  # noqa: PLC0415

    command = typer.main.get_command(importlib.import_module(module).app)
    ctx = typer.Context(command, info_name=name)
    text = get_docs_for_click(obj=command, ctx=ctx, name=name, title=None)
    # Demote headings one level so each command sits under the page title.
    return "\n".join("#" + line if line.startswith("#") else line for line in text.splitlines())


def argparse_docs(name: str) -> str:
    text = subprocess.run([name, "--help"], check=True, capture_output=True, text=True).stdout
    return f"## `{name}`\n\n```text\n{text.rstrip()}\n```\n"


def main() -> int:
    parts = [
        "# Command-line reference\n",
        "Generated from the installed commands by `scripts/gen_cli_docs.py`. "
        "Every command also prints this text with `--help`.\n",
    ]
    parts += [typer_docs(name, module) for name, module in TYPER_APPS]
    parts += [argparse_docs(name) for name in ARGPARSE_APPS]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n\n".join(parts).rstrip() + "\n")
    print(f"Wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
