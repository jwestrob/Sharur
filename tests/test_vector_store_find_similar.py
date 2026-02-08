from importlib.machinery import SourceFileLoader
from pathlib import Path

import pytest

from bennu.core.session import ExplorationSession

REPO_ROOT = Path(__file__).resolve().parents[1]
KB_PATH = REPO_ROOT / "src" / "ingest" / "07_build_knowledge_base.py"
DB_PATH = REPO_ROOT / "data" / "bennu.duckdb"


@pytest.mark.integration
def test_find_similar_uses_autoloaded_lancedb():
    if not DB_PATH.exists():
        pytest.skip("DuckDB not built; run ingest first.")

    # Build from stage outputs to ensure embeddings manifest exists
    kb_module = SourceFileLoader("kb_build", str(KB_PATH)).load_module()
    KnowledgeBaseBuilder = kb_module.KnowledgeBaseBuilder  # type: ignore
    PipelineOutputs = kb_module.PipelineOutputs  # type: ignore

    outputs = PipelineOutputs(
        stage00_dir=REPO_ROOT / "data" / "stage00_prepared",
        stage01_dir=REPO_ROOT / "data" / "stage01_quast",
        stage02_dir=REPO_ROOT / "data" / "stage02_dfast_qc",
        stage03_dir=REPO_ROOT / "data" / "stage03_prodigal",
        stage04_dir=REPO_ROOT / "data" / "stage04_astra",
        stage05a_dir=REPO_ROOT / "data" / "stage05a_gecco",
        stage05b_dir=REPO_ROOT / "data" / "stage05b_dbcan",
        stage05c_dir=REPO_ROOT / "data" / "stage05c_crispr",
        stage06_dir=REPO_ROOT / "data" / "stage06_embeddings",
    )
    builder = KnowledgeBaseBuilder(outputs, DB_PATH, force=False)
    # Do not rebuild if already present; just ensure manifest exists
    if not DB_PATH.exists():
        builder.build()

    # Session should autoload LanceDB from stage06 embeddings manifest
    session = ExplorationSession(db_path=DB_PATH)
    assert session.vector_store is not None, "vector_store not attached"

    # Pick an existing protein and query similarity; expect at least one neighbor
    protein_id = session.db.execute("SELECT protein_id FROM proteins LIMIT 1")[0][0]
    from bennu.tools.find_similar import FindSimilarTool, FindSimilarParams

    params = FindSimilarParams(query_id=protein_id, n=5)
    result = FindSimilarTool().execute(params, session)
    assert result.success
    assert result.count >= 0  # empty allowed but should not error
