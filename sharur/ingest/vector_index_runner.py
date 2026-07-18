"""Build persistent protein-vector sidecars in a Torch-free process."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from sharur.storage.vector_store import build_vector_index


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--embeddings", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=None)
    parser.add_argument("--force", action="store_true")
    return parser


def main() -> int:
    args = _parser().parse_args()
    inspection = build_vector_index(
        args.embeddings,
        force=args.force,
        threads=args.threads,
    )
    print(json.dumps(inspection.to_dict(), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
