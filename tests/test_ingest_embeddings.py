from importlib.machinery import SourceFileLoader
from pathlib import Path
import json

REPO_ROOT = Path(__file__).resolve().parents[1]
KB_PATH = REPO_ROOT / "src" / "ingest" / "07_build_knowledge_base.py"
kb_module = SourceFileLoader("kb_build_embed", str(KB_PATH)).load_module()
KnowledgeBaseBuilder = kb_module.KnowledgeBaseBuilder
PipelineOutputs = kb_module.PipelineOutputs


def test_embeddings_manifest_updates_stats(tmp_path):
    # Prepare stage dirs
    stage06 = tmp_path / "stage06_embeddings"
    stage06.mkdir(parents=True)
    manifest = {
        "total_proteins": 5,
        "output_files": {"lancedb": str(stage06 / "lancedb")},
    }
    (stage06 / "embedding_manifest.json").write_text(json.dumps(manifest))

    outputs = PipelineOutputs(
        stage00_dir=tmp_path / "stage00_prepared",
        stage01_dir=tmp_path / "stage01_quast",
        stage02_dir=tmp_path / "stage02_dfast_qc",
        stage03_dir=tmp_path / "stage03_prodigal",
        stage04_dir=tmp_path / "stage04_astra",
        stage05a_dir=tmp_path / "stage05a_gecco",
        stage05b_dir=tmp_path / "stage05b_dbcan",
        stage05c_dir=tmp_path / "stage05c_crispr",
        stage06_dir=stage06,
    )
    db_path = tmp_path / "bennu.duckdb"
    builder = KnowledgeBaseBuilder(outputs, db_path, force=True)
    stats = builder.build()
    assert stats["embeddings"] == 5
    assert builder.embeddings_path.endswith("lancedb")
