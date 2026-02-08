from importlib.machinery import SourceFileLoader
from pathlib import Path
import shutil

import duckdb
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "data"
KB_PATH = REPO_ROOT / "src" / "ingest" / "07_build_knowledge_base.py"


@pytest.mark.integration
def test_crispr_minced_ingests_into_loci(tmp_path):
    minced_path = shutil.which("minced")
    if not minced_path:
        pytest.skip("minced not installed")

    # Use existing stage outputs for 00,03,04,05a,05b,06; generate fresh 05c via minced
    stage05c = tmp_path / "stage05c_crispr"
    stage05c.mkdir(parents=True, exist_ok=True)

    # Run minced_crispr
    cmd = [
        str(REPO_ROOT / "src" / "ingest" / "minced_crispr.py"),
        "--input-dir",
        str(DATA_DIR / "stage00_prepared"),
        "--output-dir",
        str(stage05c),
        "--force",
    ]
    import subprocess

    subprocess.run(["python"] + cmd, check=True)

    # Build KB pointing to the minced outputs
    kb_module = SourceFileLoader("kb_build", str(KB_PATH)).load_module()
    KnowledgeBaseBuilder = kb_module.KnowledgeBaseBuilder  # type: ignore
    PipelineOutputs = kb_module.PipelineOutputs  # type: ignore

    outputs = PipelineOutputs(
        stage00_dir=DATA_DIR / "stage00_prepared",
        stage01_dir=DATA_DIR / "stage01_quast",
        stage02_dir=DATA_DIR / "stage02_dfast_qc",
        stage03_dir=DATA_DIR / "stage03_prodigal",
        stage04_dir=DATA_DIR / "stage04_astra",
        stage05a_dir=DATA_DIR / "stage05a_gecco",
        stage05b_dir=DATA_DIR / "stage05b_dbcan",
        stage05c_dir=stage05c,
        stage06_dir=DATA_DIR / "stage06_embeddings",
    )

    db_path = tmp_path / "bennu.duckdb"
    builder = KnowledgeBaseBuilder(outputs, db_path, force=True)
    stats = builder.build()

    conn = duckdb.connect(str(db_path))
    crispr_count = conn.execute("SELECT COUNT(*) FROM loci WHERE locus_type='crispr'").fetchone()[0]
    # If minced produced arrays, ensure we ingested them; otherwise allow zero.
    expected = 0
    for jf in stage05c.glob("*_crispr_arrays.json"):
        data = json.loads(jf.read_text())
        expected += len(data.get("arrays", []))
    assert crispr_count == expected
