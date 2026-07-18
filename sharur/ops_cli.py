"""Dependency-aware launcher for the optional Sharur Ops HTTP service."""

from __future__ import annotations


def main() -> None:
    try:
        from sharur.ops.server import main as server_main  # noqa: PLC0415 - optional extra

        server_main()
    except ModuleNotFoundError as exc:
        if exc.name not in {"fastapi", "uvicorn"}:
            raise
        raise SystemExit(
            "Sharur Ops HTTP dependencies are not installed. "
            'Install them with: pip install "sharur[ops]"'
        ) from exc


__all__ = ["main"]


if __name__ == "__main__":
    main()
