"""Atomic, compact evidence-bundle export for biological cases."""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING, Any

from sharur import __version__


if TYPE_CHECKING:
    from sharur.core.case_models import ContextComparison
    from sharur.core.claim_compiler import FindingDraft
    from sharur.operators.cases import BiologicalCase


EVIDENCE_BUNDLE_SCHEMA_VERSION = "1.0"


def _sha256_file(path: Path, *, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def _write_json(path: Path, payload: Any) -> None:
    path.write_text(
        json.dumps(payload, indent=2, default=str, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _source_fingerprints() -> dict[str, str]:
    """Hash the small code surfaces needed to interpret/replay a case."""
    package_root = Path(__file__).resolve().parent
    relative_paths = (
        "assembly_evidence.py",
        "core/case_models.py",
        "core/claim_compiler.py",
        "operators/cases.py",
    )
    return {
        relative_path: _sha256_file(package_root / relative_path)
        for relative_path in relative_paths
        if (package_root / relative_path).is_file()
    }


def _sequence_fasta(case: BiologicalCase) -> str:
    protein_ids = list(dict.fromkeys(case.record.components))
    if not protein_ids:
        return ""
    placeholders = ", ".join(["?"] * len(protein_ids))
    rows = case.store.execute(
        f"""
        SELECT protein_id, sequence
        FROM proteins
        WHERE protein_id IN ({placeholders})
        ORDER BY protein_id
        """,
        protein_ids,
    )
    chunks = []
    for protein_id, sequence in rows:
        if not sequence:
            continue
        chunks.append(f">{protein_id}\n")
        text = str(sequence)
        chunks.extend(text[index : index + 60] + "\n" for index in range(0, len(text), 60))
    return "".join(chunks)


def _verification_sql(finding_payload: dict[str, Any] | None) -> str:
    if finding_payload is None:
        return "-- No finding draft was supplied.  See case.json and comparison_recipe.json.\n"
    blocks = []
    for index, verification in enumerate(
        finding_payload.get("verification", []),
        start=1,
    ):
        query = str(verification.get("query", "")).strip()
        blocks.extend(
            [
                f"-- Verification {index}: {verification.get('claim', '')}",
                "-- Expected: " + json.dumps(verification.get("expected"), default=str),
            ]
        )
        if query.upper().startswith(("SELECT", "WITH")):
            blocks.append(query.rstrip(";") + ";")
        else:
            blocks.append(f"-- Replay instruction: {query}")
        blocks.append("")
    return "\n".join(blocks).rstrip() + "\n"


def _comparison_verifier(
    case: BiologicalCase,
    comparison: ContextComparison,
) -> str:
    recipe = dict(comparison.recipe)
    for key in ("method", "case_entity_id", "case_entity_type", "source_table"):
        recipe.pop(key, None)
    expected = {
        "foreground_positive": comparison.composite.foreground_positive,
        "foreground_total": comparison.composite.foreground_total,
        "background_positive": comparison.composite.background_positive,
        "background_total": comparison.composite.background_total,
        "fisher_p": comparison.composite.fisher_p,
    }
    entity = case.record.entity
    return f'''#!/usr/bin/env python3
"""Replay the context comparison captured in this evidence bundle."""

import argparse
import json
import math

from sharur.operators import Sharur


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", required=True, help="Path to sharur.duckdb")
    args = parser.parse_args()

    sharur = Sharur(args.db, read_only=True)
    case = sharur.inspect(
        {entity.entity_id!r},
        entity_type={entity.entity_type.value!r},
        source_table={entity.source_table!r},
    )
    observed = case.compare_context(**{recipe!r}).composite
    expected = {expected!r}
    actual = {{
        "foreground_positive": observed.foreground_positive,
        "foreground_total": observed.foreground_total,
        "background_positive": observed.background_positive,
        "background_total": observed.background_total,
        "fisher_p": observed.fisher_p,
    }}
    for key in expected:
        if key == "fisher_p":
            assert math.isclose(actual[key], expected[key], rel_tol=1e-12, abs_tol=1e-15)
        else:
            assert actual[key] == expected[key]
    print(json.dumps(actual, indent=2))


if __name__ == "__main__":
    main()
'''


def _bundle_readme(
    case: BiologicalCase,
    *,
    has_comparison: bool,
    has_finding: bool,
    has_sequences: bool,
    has_plot: bool,
) -> str:
    lines = [
        f"# Evidence bundle: {case.record.entity.entity_id}",
        "",
        "This directory is a compact, replayable Sharur evidence packet. It "
        "does not contain a copy of the dataset DuckDB or any embedding index.",
        "",
        "## Contents",
        "",
        "- `case.json`: typed entity, neighborhood, provenance, and limitations.",
        "- `case.md`: human-readable case summary.",
    ]
    if has_comparison:
        lines.extend(
            [
                "- `comparison.json`: complete foreground/background result.",
                "- `comparison_recipe.json`: arguments needed to replay it.",
                "- `verify_comparison.py`: executable result check.",
            ]
        )
    if has_finding:
        lines.extend(
            [
                "- `finding.json`: claim-compiled canonical finding draft.",
                "- `verification.sql`: directly executable SQL checks and replay instructions.",
            ]
        )
    if has_sequences:
        lines.append("- `components.faa`: anchor-component sequences only.")
    if has_plot:
        lines.append("- `case.png`: plot rendered from the resolved case.")
    lines.extend(
        [
            "- `manifest.json`: file hashes, dataset identity, and software version.",
            "",
            "## Replay",
            "",
        ]
    )
    if has_comparison:
        lines.append("Run `python verify_comparison.py --db /path/to/sharur.duckdb`.")
    else:
        lines.append("Inspect `case.json` against the dataset recorded in `manifest.json`.")
    return "\n".join(lines) + "\n"


def export_evidence_bundle(
    case: BiologicalCase,
    output_dir: str | Path,
    *,
    comparison: ContextComparison | None = None,
    finding: FindingDraft | None = None,
    include_sequences: bool = True,
    include_plot: bool = False,
    overwrite: bool = False,
) -> Path:
    """Write a complete evidence packet through a temporary sibling directory."""
    output_dir = Path(output_dir).expanduser().resolve()
    if output_dir.exists() and not overwrite:
        raise FileExistsError(output_dir)
    if output_dir.exists() and not output_dir.is_dir():
        raise ValueError(f"Evidence-bundle destination is not a directory: {output_dir}")
    finding_payload = finding.compile(strict=True) if finding is not None else None
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_dir.with_name(f".{output_dir.name}.{uuid.uuid4().hex}.tmp")
    temporary.mkdir()
    backup: Path | None = None
    try:
        _write_json(temporary / "case.json", case.to_dict())
        (temporary / "case.md").write_text(
            case.to_markdown() + "\n",
            encoding="utf-8",
        )

        if comparison is not None:
            _write_json(
                temporary / "comparison.json",
                comparison.model_dump(mode="json"),
            )
            _write_json(temporary / "comparison_recipe.json", comparison.recipe)
            verifier = temporary / "verify_comparison.py"
            verifier.write_text(
                _comparison_verifier(case, comparison),
                encoding="utf-8",
            )
            verifier.chmod(0o755)

        if finding_payload is not None:
            _write_json(temporary / "finding.json", finding_payload)
            (temporary / "verification.sql").write_text(
                _verification_sql(finding_payload),
                encoding="utf-8",
            )

        fasta = _sequence_fasta(case) if include_sequences else ""
        if fasta:
            (temporary / "components.faa").write_text(fasta, encoding="utf-8")

        if include_plot:
            case.plot(temporary / "case.png")

        (temporary / "README.md").write_text(
            _bundle_readme(
                case,
                has_comparison=comparison is not None,
                has_finding=finding_payload is not None,
                has_sequences=bool(fasta),
                has_plot=include_plot,
            ),
            encoding="utf-8",
        )

        db_path = getattr(case.store, "db_path", None)
        dataset: dict[str, Any] = {
            "path": str(db_path) if db_path is not None else "memory",
        }
        if db_path is not None and Path(db_path).is_file():
            stat = Path(db_path).stat()
            dataset.update(
                {
                    "size": stat.st_size,
                    "mtime_ns": stat.st_mtime_ns,
                }
            )
            seal_path = Path(db_path).parent / "dataset.seal.json"
            if seal_path.is_file():
                dataset["seal"] = {
                    "path": str(seal_path),
                    "sha256": _sha256_file(seal_path),
                }

        files = []
        for path in sorted(temporary.iterdir()):
            if not path.is_file() or path.name == "manifest.json":
                continue
            files.append(
                {
                    "path": path.name,
                    "size": path.stat().st_size,
                    "sha256": _sha256_file(path),
                }
            )
        manifest = {
            "schema_version": EVIDENCE_BUNDLE_SCHEMA_VERSION,
            "created_at": datetime.now(timezone.utc).isoformat(),
            "sharur_version": __version__,
            "source_sha256": _source_fingerprints(),
            "entity": case.record.entity.model_dump(mode="json"),
            "dataset": dataset,
            "assembly_evidence_sidecar": (
                str(case.assembly_evidence_path)
                if case.assembly_evidence_path is not None
                else None
            ),
            "files": files,
        }
        _write_json(temporary / "manifest.json", manifest)

        if output_dir.exists():
            backup = output_dir.with_name(f".{output_dir.name}.{uuid.uuid4().hex}.backup")
            os.replace(output_dir, backup)
        try:
            os.replace(temporary, output_dir)
        except Exception:
            if backup is not None and backup.exists() and not output_dir.exists():
                os.replace(backup, output_dir)
            raise
        if backup is not None and backup.exists():
            shutil.rmtree(backup)
    except Exception:
        if temporary.exists():
            shutil.rmtree(temporary)
        if backup is not None and backup.exists() and not output_dir.exists():
            os.replace(backup, output_dir)
        raise
    return output_dir


__all__ = [
    "EVIDENCE_BUNDLE_SCHEMA_VERSION",
    "export_evidence_bundle",
]
