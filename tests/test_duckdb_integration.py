from importlib.machinery import SourceFileLoader
from pathlib import Path

import duckdb
import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "data"
KB_PATH = REPO_ROOT / "src" / "ingest" / "07_build_knowledge_base.py"


def _stage_dirs_exist() -> bool:
    expected = [
        "stage00_prepared",
        "stage01_quast",
        "stage03_prodigal",
        "stage04_astra",
        "stage05a_gecco",
        "stage05b_dbcan",
        "stage06_embeddings",
    ]
    return all((DATA_DIR / d).exists() for d in expected)


@pytest.mark.integration
def test_build_and_query_duckdb(tmp_path):
    if not _stage_dirs_exist():
        pytest.skip("Stage outputs not present; run ingest first.")

    # Load builder
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
        stage05c_dir=DATA_DIR / "stage05c_crispr",
        stage06_dir=DATA_DIR / "stage06_embeddings",
    )

    db_path = tmp_path / "sharur.duckdb"
    builder = KnowledgeBaseBuilder(outputs, db_path, force=True)
    stats = builder.build()

    # Expected counts for dummy_dataset run
    assert stats["bins"] == 8
    assert stats["contigs"] == 593
    assert stats["proteins"] == 10102
    assert stats["annotations"] >= 1700  # pfam+kofam+cazy combined
    assert stats["loci"] == 2  # GECCO clusters
    assert stats["embeddings"] == 10102

    conn = duckdb.connect(str(db_path))

    # No duplicate contig IDs
    dup_contigs = conn.execute(
        "SELECT COUNT(*) FROM contigs GROUP BY contig_id HAVING COUNT(*)>1"
    ).fetchall()
    assert dup_contigs == []

    # All protein bin_ids exist in bins
    missing_bin = conn.execute(
        """
        SELECT COUNT(*) FROM proteins p
        LEFT JOIN bins b ON p.bin_id = b.bin_id
        WHERE b.bin_id IS NULL
        """
    ).fetchone()[0]
    assert missing_bin == 0

    # View matches filter
    view_count = conn.execute("SELECT COUNT(*) FROM bgc_loci").fetchone()[0]
    filter_count = conn.execute("SELECT COUNT(*) FROM loci WHERE locus_type='bgc'").fetchone()[0]
    assert view_count == filter_count == 2

    # PFAM/KOFAM annotations present
    pfam_kofam = conn.execute(
        "SELECT COUNT(*) FROM annotations WHERE source IN ('pfam','kofam')"
    ).fetchone()[0]
    assert pfam_kofam > 0

    # CAZy annotations present
    cazy = conn.execute(
        "SELECT COUNT(*) FROM annotations WHERE source='cazy'"
    ).fetchone()[0]
    assert cazy > 0
