import json
from importlib.machinery import SourceFileLoader
from pathlib import Path

import duckdb

REPO_ROOT = Path(__file__).resolve().parents[1]
KB_PATH = REPO_ROOT / "src" / "ingest" / "07_build_knowledge_base.py"
kb_module = SourceFileLoader("kb_build", str(KB_PATH)).load_module()
KnowledgeBaseBuilder = kb_module.KnowledgeBaseBuilder
PipelineOutputs = kb_module.PipelineOutputs


def _write(path: Path, text: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def test_ingest_smoke(tmp_path):
    data_dir = tmp_path / "data"
    # Stage02 manifest with one genome
    stage02 = data_dir / "stage02_dfast_qc"
    stage02.mkdir(parents=True)
    manifest = {
        "genomes": [
            {
                "genome_id": "bin1",
                "taxonomy": {"name": "Bac", "completeness": 90.0, "contamination": 1.0},
                "n_contigs": 1,
                "total_length": 1000,
            }
        ]
    }
    _write(stage02 / "processing_manifest.json", json.dumps(manifest))

    # Stage03 Prodigal .faa
    stage03 = data_dir / "stage03_prodigal" / "genomes" / "bin1"
    stage03.mkdir(parents=True)
    _write(stage03 / "bin1.faa", ">bin1_contig1_00001 # 10 # 40 # 1\nMAAA\n")

    # Stage04 Astra annotation TSV
    stage04 = data_dir / "stage04_astra"
    stage04.mkdir(parents=True)
    _write(
        stage04 / "bin1_hits_df.tsv",
        "sequence_id\thmm_name\te_value\tbit_score\tali_from\tali_to\tname\tdescription\n"
        "bin1_contig1_00001\tPF00001\t1e-5\t50\t1\t10\tName\tDesc\n",
    )

    # Stage05a GECCO
    stage05a = data_dir / "stage05a_gecco"
    stage05a.mkdir(parents=True)
    _write(
        stage05a / "combined_bgc_data.json",
        json.dumps(
            {
                "clusters": [
                    {
                        "cluster_id": "cluster1",
                        "contig": "bin1_contig1",
                        "start": 5,
                        "end": 60,
                        "bgc_type": "nrps",
                        "protein_list": ["bin1_contig1_00001"],
                    }
                ]
            }
        ),
    )

    # Stage05c CRISPR (optional)
    stage05c = data_dir / "stage05c_crispr"
    stage05c.mkdir(parents=True)
    _write(
        stage05c / "bin1_crispr_arrays.json",
        json.dumps({"arrays": []}),
    )

    outputs = PipelineOutputs(
        stage00_dir=data_dir / "stage00_prepared",
        stage01_dir=data_dir / "stage01_quast",
        stage02_dir=stage02,
        stage03_dir=data_dir / "stage03_prodigal",
        stage04_dir=stage04,
        stage05a_dir=stage05a,
        stage05b_dir=data_dir / "stage05b_dbcan",
        stage05c_dir=stage05c,
        stage06_dir=data_dir / "stage06_embeddings",
    )

    db_path = data_dir / "bennu.duckdb"
    builder = KnowledgeBaseBuilder(outputs, db_path, force=True)
    stats = builder.build()

    assert stats["proteins"] == 1
    assert stats["annotations"] >= 1
    assert stats["loci"] >= 1

    # Quick SQL check
    conn = duckdb.connect(str(db_path))
    count = conn.execute("SELECT COUNT(*) FROM proteins").fetchone()[0]
    assert count == 1
