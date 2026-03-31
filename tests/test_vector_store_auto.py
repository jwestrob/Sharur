from pathlib import Path
import json

import h5py
import numpy as np

from sharur.core.session import ExplorationSession


def test_session_auto_loads_faiss_if_h5_present(tmp_path):
    db_path = tmp_path / "sharur.duckdb"
    db_path.touch()
    embeddings_dir = tmp_path / "embeddings"
    embeddings_dir.mkdir(parents=True)

    # Write a small H5 embeddings file
    h5_path = embeddings_dir / "protein_embeddings.h5"
    with h5py.File(h5_path, "w") as f:
        f.create_dataset("protein_ids", data=["prot_1", "prot_2", "prot_3"])
        f.create_dataset("embeddings", data=np.random.randn(3, 320).astype(np.float32))

    session = ExplorationSession(db_path=db_path)
    assert session.vector_store is not None
