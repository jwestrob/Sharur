#!/usr/bin/env python3
"""Repo-local shim for the packaged `sharur-ingest` CLI."""

from sharur.ingest_cli import app, main, run


__all__ = ["app", "main", "run"]


if __name__ == "__main__":
    main()
