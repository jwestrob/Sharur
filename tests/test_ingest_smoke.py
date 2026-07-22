import json
from importlib.machinery import SourceFileLoader
from pathlib import Path

import duckdb
import pytest
from pandas.errors import EmptyDataError

import sharur.colocation


REPO_ROOT = Path(__file__).resolve().parents[1]
KB_PATH = REPO_ROOT / "src" / "ingest" / "07_build_knowledge_base.py"


def _load_kb_module():
    """Load the knowledge base module lazily to avoid slowing pytest collection."""
    kb_module = SourceFileLoader("kb_build", str(KB_PATH)).load_module()
    return kb_module.KnowledgeBaseBuilder, kb_module.PipelineOutputs


def _pipeline_outputs(outputs_cls, data_dir: Path):
    return outputs_cls(
        stage00_dir=data_dir / "stage00_prepared",
        stage01_dir=data_dir / "stage01_quast",
        stage02_dir=data_dir / "stage02_dfast_qc",
        stage03_dir=data_dir / "stage03_prodigal",
        stage04_dir=data_dir / "stage04_astra",
        stage05a_dir=data_dir / "stage05a_gecco",
        stage05b_dir=data_dir / "stage05b_dbcan",
        stage05c_dir=data_dir / "stage05c_crispr",
        stage06_dir=data_dir / "stage06_embeddings",
    )


def _write(path: Path, text: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def test_ingest_smoke(tmp_path):
    builder_cls, outputs_cls = _load_kb_module()
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

    outputs = outputs_cls(
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

    db_path = data_dir / "sharur.duckdb"
    builder = builder_cls(outputs, db_path, force=True)
    stats = builder.build()

    assert stats["proteins"] == 1
    assert stats["annotations"] >= 1
    assert stats["loci"] >= 1
    assert stats["predicates"] == 1
    assert stats["semantic_state"] == 1
    assert stats["semantic_atoms"] >= 1

    # Quick SQL check
    conn = duckdb.connect(str(db_path))
    count = conn.execute("SELECT COUNT(*) FROM proteins").fetchone()[0]
    assert count == 1
    semantic_state = conn.execute("SELECT COUNT(*) FROM semantic_state").fetchone()[0]
    semantic_atoms = conn.execute("SELECT COUNT(*) FROM semantic_atoms").fetchone()[0]
    legacy = conn.execute("SELECT COUNT(*) FROM protein_predicates").fetchone()[0]
    assert semantic_state == 1
    assert semantic_atoms >= 1
    assert legacy == 1
    assert (data_dir / "reports" / "predicates_v2_review_queue.tsv").exists()


def test_stage07_uses_slurm_threads_for_every_db_connection(tmp_path, monkeypatch):
    builder_cls, outputs_cls = _load_kb_module()
    monkeypatch.setenv("SLURM_CPUS_ON_NODE", "7")
    builder = builder_cls(
        _pipeline_outputs(outputs_cls, tmp_path),
        tmp_path / "sharur.duckdb",
        force=True,
    )

    assert builder.threads == 7
    builder._init_db()
    assert int(builder.conn.execute("SELECT current_setting('threads')").fetchone()[0]) == 7
    builder._release_db()
    builder._reacquire_db()
    assert int(builder.conn.execute("SELECT current_setting('threads')").fetchone()[0]) == 7
    builder._release_db()


def test_stage07_explicit_threads_override_slurm(tmp_path, monkeypatch):
    builder_cls, outputs_cls = _load_kb_module()
    monkeypatch.setenv("SLURM_CPUS_ON_NODE", "48")

    builder = builder_cls(
        _pipeline_outputs(outputs_cls, tmp_path),
        tmp_path / "sharur.duckdb",
        threads=5,
    )

    assert builder.threads == 5


def test_stage07_rejects_invalid_slurm_thread_count(tmp_path, monkeypatch):
    builder_cls, outputs_cls = _load_kb_module()
    monkeypatch.setenv("SLURM_CPUS_ON_NODE", "many")

    with pytest.raises(ValueError, match="SLURM_CPUS_ON_NODE must be a positive integer"):
        builder_cls(
            _pipeline_outputs(outputs_cls, tmp_path),
            tmp_path / "sharur.duckdb",
        )


def test_stage07_colocation_failure_aborts_build(tmp_path, monkeypatch):
    builder_cls, outputs_cls = _load_kb_module()
    builder = builder_cls(
        _pipeline_outputs(outputs_cls, tmp_path),
        tmp_path / "sharur.duckdb",
        threads=1,
    )

    class FakeResult:
        @staticmethod
        def fetchone():
            return (1,)

    class FakeConnection:
        @staticmethod
        def execute(*_args, **_kwargs):
            return FakeResult()

    def fail_validation(**_kwargs):
        raise RuntimeError("synthetic caller failure")

    builder.conn = FakeConnection()
    monkeypatch.setattr(builder, "_release_db", lambda: None)
    monkeypatch.setattr(builder, "_reacquire_db", lambda: None)
    monkeypatch.setattr(sharur.colocation, "validate_systems", fail_validation)

    with pytest.raises(RuntimeError, match="synthetic caller failure"):
        builder._run_colocation("defensefinder", "defense")


def test_stage07_annotation_chunk_failure_aborts_build(tmp_path):
    builder_cls, outputs_cls = _load_kb_module()
    outputs = _pipeline_outputs(outputs_cls, tmp_path)
    outputs.stage04_dir.mkdir(parents=True)
    _write(outputs.stage04_dir / "PFAM_broken_hits_df.tsv", "")
    builder = builder_cls(outputs, tmp_path / "sharur.duckdb", threads=1)
    builder._init_db()

    try:
        with pytest.raises(EmptyDataError):
            builder._load_annotations()
    finally:
        builder._release_db()
