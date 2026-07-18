"""Generate tiny deterministic stage artifacts for ``sharur-ingest --mode fast``."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import TYPE_CHECKING


if TYPE_CHECKING:
    from collections.abc import Iterator


def _parse_contigs(fasta_path: Path) -> Iterator[tuple[str, str]]:
    header = None
    sequence_parts: list[str] = []
    with fasta_path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence_parts)
                header = line[1:].split()[0]
                sequence_parts = []
            else:
                sequence_parts.append(line)
    if header is not None:
        yield header, "".join(sequence_parts)


def synthesize_pipeline(input_dir: Path, data_dir: Path) -> None:
    """Write the minimal Stage 02-05 artifacts consumed by Stage 07."""
    stage02 = data_dir / "stage02_dfast_qc"
    stage03 = data_dir / "stage03_prodigal"
    stage04 = data_dir / "stage04_astra"
    stage05a = data_dir / "stage05a_gecco"
    stage05b = data_dir / "stage05b_dbcan"
    stage05c = data_dir / "stage05c_crispr"
    for directory in (stage02, stage03, stage04, stage05a, stage05b, stage05c):
        directory.mkdir(parents=True, exist_ok=True)

    fasta_paths = [
        path
        for pattern in ("*.fna", "*.fa", "*.fasta")
        for path in sorted(input_dir.glob(pattern))
    ]
    if not fasta_paths:
        raise FileNotFoundError(f"No FASTA files found in {input_dir}")

    genomes_manifest: dict[str, list[dict]] = {"genomes": []}
    pfam_rows: list[dict] = []
    gecco_clusters: list[dict] = []
    for fasta_path in fasta_paths:
        raw_bin_id = fasta_path.stem.replace(".contigs", "")
        bin_id = raw_bin_id.replace("_", "")
        contigs = list(_parse_contigs(fasta_path))
        if not contigs:
            raise ValueError(f"No FASTA records found in {fasta_path}")
        genomes_manifest["genomes"].append(
            {
                "genome_id": bin_id,
                "taxonomy": {
                    "name": "unknown",
                    "completeness": 90.0,
                    "contamination": 1.0,
                },
                "n_contigs": len(contigs),
                "total_length": sum(len(sequence) for _, sequence in contigs),
            }
        )

        bin_dir = stage03 / "genomes" / bin_id
        bin_dir.mkdir(parents=True, exist_ok=True)
        faa_path = bin_dir / f"{bin_id}.faa"
        with faa_path.open("w") as faa:
            for index, (_contig_id, sequence) in enumerate(contigs):
                protein_id = f"{bin_id}_ctg{index:05d}_00001"
                faa.write(
                    f">{protein_id} # 1 # {min(len(sequence), 300)} # 1 "
                    f"# ID={protein_id};partial=00\n"
                )
                faa.write("M" * min(len(sequence), 100) + "\n")
                if index == 0:
                    pfam_rows.append(
                        {
                            "sequence_id": protein_id,
                            "hmm_name": "PF00001",
                            "e_value": 1e-5,
                            "bit_score": 42.0,
                            "ali_from": 1,
                            "ali_to": min(30, len(sequence)),
                            "name": "MockDomain",
                            "description": "Synthetic annotation",
                        }
                    )
                    gecco_clusters.append(
                        {
                            "cluster_id": f"{bin_id}_cluster1",
                            "contig": protein_id.rsplit("_", 1)[0],
                            "start": 1,
                            "end": 100,
                            "bgc_type": "nrps",
                            "protein_list": [protein_id],
                        }
                    )

    (stage02 / "processing_manifest.json").write_text(json.dumps(genomes_manifest))
    header = [
        "sequence_id",
        "hmm_name",
        "e_value",
        "bit_score",
        "ali_from",
        "ali_to",
        "name",
        "description",
    ]
    lines = ["\t".join(header)]
    lines.extend("\t".join(str(row[column]) for column in header) for row in pfam_rows)
    (stage04 / "synthetic_hits_df.tsv").write_text("\n".join(lines) + "\n")
    (stage05a / "combined_bgc_data.json").write_text(
        json.dumps({"clusters": gecco_clusters})
    )
    (stage05b / "synthetic_cazyme_results.json").write_text(
        json.dumps({"annotations": []})
    )
    (stage05c / "synthetic_crispr_arrays.json").write_text(
        json.dumps({"arrays": []})
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--data-dir", required=True, type=Path)
    args = parser.parse_args()
    synthesize_pipeline(
        args.input_dir.expanduser().resolve(),
        args.data_dir.expanduser().resolve(),
    )


__all__ = ["main", "synthesize_pipeline"]


if __name__ == "__main__":
    main()
