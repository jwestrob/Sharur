"""Integration test for ingest pipeline using dummy_dataset genomes."""

import json
from importlib.machinery import SourceFileLoader
from pathlib import Path
from typing import Iterable, Tuple

import duckdb
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
KB_PATH = REPO_ROOT / "src" / "ingest" / "07_build_knowledge_base.py"
kb_module = SourceFileLoader("kb_build_full", str(KB_PATH)).load_module()
KnowledgeBaseBuilder = kb_module.KnowledgeBaseBuilder
PipelineOutputs = kb_module.PipelineOutputs


def _parse_contigs(fasta_path: Path) -> Iterable[Tuple[str, str]]:
    """Yield (contig_id, sequence) pairs from a FASTA file."""
    header = None
    seq_parts = []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
        if header is not None:
            yield header, "".join(seq_parts)


def _write(path: Path, text: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def test_ingest_pipeline_with_dummy_dataset(tmp_path):
    dataset_dir = REPO_ROOT / "dummy_dataset"
    if not dataset_dir.exists():
        pytest.skip("dummy_dataset not available")

    data_dir = tmp_path / "data"
    stage02 = data_dir / "stage02_dfast_qc"
    stage03 = data_dir / "stage03_prodigal"
    stage04 = data_dir / "stage04_astra"
    stage05a = data_dir / "stage05a_gecco"
    stage05b = data_dir / "stage05b_dbcan"
    stage05c = data_dir / "stage05c_crispr"

    genomes_manifest = {"genomes": []}
    pfam_rows = []
    cazy_rows = []
    gecco_clusters = []
    total_proteins = 0
    total_contigs = 0

    for fasta_path in sorted(dataset_dir.glob("*.fna")):
        raw_bin_id = fasta_path.name.replace(".contigs.fna", "")
        # Sanitize bin_id so split('_')[0] matches (mirrors _row_from_header behavior)
        bin_id = raw_bin_id.replace("_", "")
        contigs = list(_parse_contigs(fasta_path))
        total_length = sum(len(seq) for _, seq in contigs)
        total_contigs += len(contigs)

        genomes_manifest["genomes"].append(
            {
                "genome_id": bin_id,
                "taxonomy": {"name": "unknown", "completeness": 90.0, "contamination": 1.0},
                "n_contigs": len(contigs),
                "total_length": total_length,
            }
        )

        bin_dir = stage03 / "genomes" / bin_id
        bin_dir.mkdir(parents=True, exist_ok=True)
        faa_path = bin_dir / f"{bin_id}.faa"

        with open(faa_path, "w") as faa:
            for idx, (contig_id, seq) in enumerate(contigs):
                protein_id = f"{bin_id}_ctg{idx:05d}_00001"
                faa.write(f">{protein_id} # 1 # {min(len(seq), 300)} # 1 # ID={protein_id};partial=00\n")
                faa.write("M" * 30 + "\n")
                total_proteins += 1

                # Annotate the first protein per bin
                if idx == 0:
                    pfam_rows.append(
                        {
                            "sequence_id": protein_id,
                            "hmm_name": "PF00001",
                            "e_value": 1e-5,
                            "bit_score": 42.0,
                            "ali_from": 1,
                            "ali_to": 30,
                            "name": "MockDomain",
                            "description": "Synthetic annotation",
                        }
                    )
                    cazy_rows.append({"protein_id": protein_id, "cazyme_family": "GH1"})
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

    # Stage02 manifest
    stage02.mkdir(parents=True, exist_ok=True)
    _write(stage02 / "processing_manifest.json", json.dumps(genomes_manifest))

    # Stage04 PFAM/KOFAM TSV (minimal)
    if pfam_rows:
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
        for row in pfam_rows:
            lines.append("\t".join(str(row[h]) for h in header))
        stage04.mkdir(parents=True, exist_ok=True)
        _write(stage04 / "synthetic_pfam_hits_df.tsv", "\n".join(lines) + "\n")

    # Stage05b CAZy annotations
    if cazy_rows:
        stage05b.mkdir(parents=True, exist_ok=True)
        _write(stage05b / "synthetic_cazyme_results.json", json.dumps({"annotations": cazy_rows}))

    # Stage05a GECCO clusters
    stage05a.mkdir(parents=True, exist_ok=True)
    _write(stage05a / "combined_bgc_data.json", json.dumps({"clusters": gecco_clusters}))

    # Stage05c CRISPR arrays (empty)
    stage05c.mkdir(parents=True, exist_ok=True)
    _write(stage05c / "synthetic_crispr_arrays.json", json.dumps({"arrays": []}))

    outputs = PipelineOutputs(
        stage00_dir=data_dir / "stage00_prepared",
        stage01_dir=data_dir / "stage01_quast",
        stage02_dir=stage02,
        stage03_dir=stage03,
        stage04_dir=stage04,
        stage05a_dir=stage05a,
        stage05b_dir=stage05b,
        stage05c_dir=stage05c,
        stage06_dir=data_dir / "stage06_embeddings",
    )

    db_path = data_dir / "sharur.duckdb"
    builder = KnowledgeBaseBuilder(outputs, db_path, force=True)
    stats = builder.build()

    assert stats["bins"] == len(genomes_manifest["genomes"])
    assert stats["contigs"] == total_contigs
    assert stats["proteins"] == total_proteins
    assert stats["annotations"] >= len(pfam_rows) + len(cazy_rows)
    assert stats["loci"] >= len(gecco_clusters)

    # Feature store should have one row per protein
    conn = duckdb.connect(str(db_path))
    feature_count = conn.execute("SELECT COUNT(*) FROM feature_store").fetchone()[0]
    assert feature_count == total_proteins
