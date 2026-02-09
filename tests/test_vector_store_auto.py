from pathlib import Path
import json

from sharur.core.session import ExplorationSession


def test_session_auto_loads_lancedb_if_manifest_present(tmp_path):
    db_path = tmp_path / "sharur.duckdb"
    # create empty duckdb file to satisfy path existence expectation in store
    db_path.touch()
    stage06 = tmp_path / "stage06_embeddings"
    stage06.mkdir(parents=True)
    manifest = {
        "total_proteins": 3,
        "output_files": {"lancedb": str(stage06 / "lancedb")},
    }
    (stage06 / "embedding_manifest.json").write_text(json.dumps(manifest))

    session = ExplorationSession(db_path=db_path)
    assert session.vector_store is not None
