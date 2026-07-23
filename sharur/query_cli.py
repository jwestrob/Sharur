"""Dependency-aware launcher for the optional Sharur Query HTTP service."""

from __future__ import annotations


def main() -> None:
    try:
        from sharur.query.server import main as server_main  # noqa: PLC0415

        server_main()
    except ModuleNotFoundError as exc:
        if exc.name not in {"fastapi", "uvicorn"}:
            raise
        raise SystemExit(
            "Sharur Query HTTP dependencies are not installed. "
            'Install them with: pip install "sharur[ops]"'
        ) from exc


__all__ = ["main"]


if __name__ == "__main__":
    main()
