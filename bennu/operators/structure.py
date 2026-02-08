"""
Structure prediction operators using ESM3 API.

Provides structure prediction for proteins using EvolutionaryScale's ESM3 model.
Requires ESM_API_KEY environment variable to be set.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import TYPE_CHECKING, Optional

from bennu.operators.base import BennuResult, OperatorContext

if TYPE_CHECKING:
    from bennu.storage.duckdb_store import DuckDBStore


def _get_esm_client():
    """Get ESM3 API client."""
    api_key = os.environ.get("ESM_API_KEY")
    if not api_key:
        return None, "ESM_API_KEY not set. Run: export ESM_API_KEY=your_key"

    try:
        import esm.sdk
        client = esm.sdk.client(model="esm3-open-2024-03", token=api_key)
        return client, None
    except Exception as e:
        return None, f"Failed to create ESM client: {e}"


def predict_structure(
    store: "DuckDBStore",
    protein_id: str,
    output_path: Optional[str] = None,
) -> BennuResult:
    """
    Predict protein structure using ESM3.

    Args:
        store: DuckDB store
        protein_id: Protein ID to predict structure for
        output_path: Optional path to save PDB file

    Returns:
        BennuResult with structure prediction info and PDB path
    """
    params = {"protein_id": protein_id, "output_path": output_path}

    with OperatorContext("predict_structure", params) as ctx:
        # Get client
        client, error = _get_esm_client()
        if error:
            return ctx.make_result(data=error, rows=0)

        # Get sequence
        rows = store.execute(
            "SELECT sequence, sequence_length FROM proteins WHERE protein_id = ?",
            [protein_id],
        )
        if not rows or not rows[0][0]:
            return ctx.make_result(
                data=f"No sequence found for {protein_id}",
                rows=0,
            )

        sequence, length = rows[0]

        # Clean sequence - remove stop codons and non-standard characters
        standard_aa = set("ACDEFGHIKLMNPQRSTVWY")
        clean_sequence = "".join(c for c in sequence.upper() if c in standard_aa)

        if len(clean_sequence) < len(sequence):
            # Report cleaning
            removed = len(sequence) - len(clean_sequence)
            sequence = clean_sequence
            length = len(sequence)

        # Check length (ESM3 has limits)
        if length > 1024:
            return ctx.make_result(
                data=f"Sequence too long ({length} aa). ESM3 open model supports up to 1024 aa.",
                rows=0,
            )

        try:
            from esm.sdk.api import ESMProtein, GenerationConfig
            import warnings
            warnings.filterwarnings("ignore", message="Entity ID not found")

            # Create protein from cleaned sequence
            protein = ESMProtein(sequence=clean_sequence)

            # Generate structure - use fewer steps for speed
            config = GenerationConfig(track="structure", num_steps=8)
            result = client.generate(protein, config)

            # Get confidence scores
            plddt = None
            if result.plddt is not None:
                plddt = float(result.plddt.mean())
            ptm = float(result.ptm) if result.ptm is not None else None

            # Save PDB
            if output_path is None:
                # Use a safe filename
                safe_name = "".join(c if c.isalnum() or c in "_-" else "_" for c in protein_id[:40])
                output_path = f"/tmp/bennu_structure_{safe_name}.pdb"

            result.to_pdb(output_path)

            # Build result summary
            lines = [
                f"# Structure Prediction: {protein_id[:50]}...",
                f"**Length:** {length} aa",
                f"**Model:** ESM3-open-2024-03",
                "",
                "## Confidence Scores",
                f"- **pLDDT (mean):** {plddt:.2f}" if plddt else "- pLDDT: N/A",
            ]
            if ptm:
                lines.append(f"- **pTM:** {ptm:.2f}")

            lines.extend([
                "",
                f"## Output",
                f"**PDB file:** {output_path}",
                f"",
                f"View with: `open {output_path}` or upload to https://molstar.org/viewer/",
            ])

            return ctx.make_result(
                data="\n".join(lines),
                rows=1,
                raw={
                    "protein_id": protein_id,
                    "length": length,
                    "plddt_mean": float(plddt) if plddt else None,
                    "ptm": float(ptm) if ptm else None,
                    "pdb_path": output_path,
                },
            )

        except Exception as e:
            return ctx.make_result(
                data=f"Structure prediction failed: {e}",
                rows=0,
            )


def predict_structure_from_sequence(
    sequence: str,
    output_path: Optional[str] = None,
    name: str = "unknown",
) -> BennuResult:
    """
    Predict structure from a raw sequence (no database lookup).

    Args:
        sequence: Amino acid sequence
        output_path: Optional path to save PDB file
        name: Name for the protein

    Returns:
        BennuResult with structure prediction info
    """
    params = {"sequence_length": len(sequence), "name": name}

    with OperatorContext("predict_structure_from_sequence", params) as ctx:
        client, error = _get_esm_client()
        if error:
            return ctx.make_result(data=error, rows=0)

        length = len(sequence)
        if length > 1024:
            return ctx.make_result(
                data=f"Sequence too long ({length} aa). ESM3 open model supports up to 1024 aa.",
                rows=0,
            )

        try:
            from esm.sdk.api import ESMProtein, GenerationConfig

            protein = ESMProtein(sequence=sequence)
            config = GenerationConfig(track="structure", num_steps=8)
            result = client.generate(protein, config)

            plddt = result.plddt.mean() if result.plddt is not None else None

            if output_path is None:
                output_path = f"/tmp/bennu_structure_{name[:30]}.pdb"

            result.to_pdb(output_path)

            lines = [
                f"# Structure Prediction: {name}",
                f"**Length:** {length} aa",
                f"**pLDDT (mean):** {plddt:.2f}" if plddt else "**pLDDT:** N/A",
                f"**PDB file:** {output_path}",
            ]

            return ctx.make_result(
                data="\n".join(lines),
                rows=1,
                raw={"pdb_path": output_path, "plddt_mean": float(plddt) if plddt else None},
            )

        except Exception as e:
            return ctx.make_result(
                data=f"Structure prediction failed: {e}",
                rows=0,
            )


def batch_predict_structures(
    store: "DuckDBStore",
    protein_ids: list[str],
    output_dir: Optional[str] = None,
    max_length: int = 1024,
) -> BennuResult:
    """
    Predict structures for multiple proteins.

    Args:
        store: DuckDB store
        protein_ids: List of protein IDs
        output_dir: Directory to save PDB files
        max_length: Skip proteins longer than this

    Returns:
        BennuResult with summary of predictions
    """
    params = {"protein_ids": protein_ids, "max_length": max_length}

    with OperatorContext("batch_predict_structures", params) as ctx:
        client, error = _get_esm_client()
        if error:
            return ctx.make_result(data=error, rows=0)

        if output_dir is None:
            output_dir = "/tmp/bennu_structures"
        Path(output_dir).mkdir(parents=True, exist_ok=True)

        # Get sequences
        placeholders = ",".join(["?"] * len(protein_ids))
        rows = store.execute(
            f"""
            SELECT protein_id, sequence, sequence_length
            FROM proteins
            WHERE protein_id IN ({placeholders})
            """,
            protein_ids,
        )

        results = []
        skipped = []

        for protein_id, sequence, length in rows:
            if not sequence:
                skipped.append((protein_id, "no sequence"))
                continue
            if length > max_length:
                skipped.append((protein_id, f"too long ({length} aa)"))
                continue

            try:
                from esm.sdk.api import ESMProtein, GenerationConfig

                protein = ESMProtein(sequence=sequence)
                config = GenerationConfig(track="structure", num_steps=8)
                result = client.generate(protein, config)

                plddt = result.plddt.mean() if result.plddt is not None else None
                pdb_path = f"{output_dir}/{protein_id[:50]}.pdb"
                result.to_pdb(pdb_path)

                results.append({
                    "protein_id": protein_id,
                    "length": length,
                    "plddt_mean": float(plddt) if plddt else None,
                    "pdb_path": pdb_path,
                })

            except Exception as e:
                skipped.append((protein_id, str(e)))

        # Build summary
        lines = [
            f"# Batch Structure Prediction",
            f"**Predicted:** {len(results)} structures",
            f"**Skipped:** {len(skipped)} proteins",
            f"**Output directory:** {output_dir}",
            "",
            "## Results",
            "| Protein | Length | pLDDT |",
            "|---------|--------|-------|",
        ]

        for r in results:
            pid_short = r["protein_id"][:30]
            plddt_str = f"{r['plddt_mean']:.2f}" if r['plddt_mean'] else "N/A"
            lines.append(f"| {pid_short}... | {r['length']} aa | {plddt_str} |")

        if skipped:
            lines.extend(["", "## Skipped"])
            for pid, reason in skipped[:5]:
                lines.append(f"- {pid[:40]}...: {reason}")
            if len(skipped) > 5:
                lines.append(f"- ... and {len(skipped) - 5} more")

        return ctx.make_result(
            data="\n".join(lines),
            rows=len(results),
            raw={"predictions": results, "skipped": skipped},
        )


__all__ = [
    "predict_structure",
    "predict_structure_from_sequence",
    "batch_predict_structures",
]
