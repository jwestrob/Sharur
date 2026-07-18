import json
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
import pytest

from sharur.core.session import ExplorationSession
from sharur.storage import vector_store
from sharur.storage.vector_store import FAISSStore, inspect_vector_index


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
    assert session.vector_store_available is True
    assert session.vector_store_path == h5_path
    assert session.vector_store is not None


def test_session_does_not_construct_faiss_until_vector_store_is_accessed(tmp_path, monkeypatch):
    db_path = tmp_path / "sharur.duckdb"
    embeddings_dir = tmp_path / "embeddings"
    embeddings_dir.mkdir(parents=True)
    h5_path = embeddings_dir / "protein_embeddings.h5"
    h5_path.touch()

    constructed = []

    class StubFAISSStore:
        def __init__(self, path):
            constructed.append(path)

    monkeypatch.setattr("sharur.core.session.FAISSStore", StubFAISSStore)

    session = ExplorationSession(db_path=db_path)
    assert session.vector_store_available is True
    assert constructed == []

    _ = session.db
    assert constructed == []

    assert session.vector_store is not None
    assert constructed == [str(h5_path)]


def test_persistent_index_is_reused_and_id_mapping_is_stable(tmp_path):
    h5_path = tmp_path / "protein_embeddings.h5"
    vectors = np.asarray(
        [
            [1.0, 0.0, 0.0],
            [0.9, 0.1, 0.0],
            [0.0, 1.0, 0.0],
        ],
        dtype=np.float32,
    )
    with h5py.File(h5_path, "w") as handle:
        handle.create_dataset(
            "protein_ids",
            data=np.asarray(["query", "nearest", "other"], dtype="S"),
        )
        handle.create_dataset("embeddings", data=vectors)

    built = FAISSStore.build_persistent_index(h5_path)
    assert built.state == "available"
    assert built.index_path is not None
    assert built.id_map_path is not None

    index_mtime = built.index_path and Path(built.index_path).stat().st_mtime_ns
    store = FAISSStore(h5_path)
    assert store.query("query", k=1, include_distances=True)[0][0] == "nearest"
    store.close()

    reused = FAISSStore.build_persistent_index(h5_path)
    assert reused.index_path == built.index_path
    assert Path(reused.index_path).stat().st_mtime_ns == index_mtime


def test_index_becomes_stale_when_canonical_h5_changes(tmp_path):
    h5_path = tmp_path / "protein_embeddings.h5"
    with h5py.File(h5_path, "w") as handle:
        handle.create_dataset("protein_ids", data=np.asarray(["a", "b"], dtype="S"))
        handle.create_dataset(
            "embeddings",
            data=np.asarray([[1.0, 0.0], [0.0, 1.0]], dtype=np.float32),
        )
    FAISSStore.build_persistent_index(h5_path)

    with h5py.File(h5_path, "a") as handle:
        handle.attrs["created_at"] = "changed"

    inspection = inspect_vector_index(h5_path)
    assert inspection.state == "stale"


def test_ivf_generation_opens_with_mmap_and_preserves_id_mapping(tmp_path, monkeypatch):
    monkeypatch.setattr(vector_store, "FLAT_INDEX_LIMIT", 2)
    rng = np.random.default_rng(7)
    vectors = rng.normal(size=(2_000, 16)).astype(np.float32)
    vectors[1] = vectors[0] + 1e-4
    protein_ids = np.asarray([f"protein-{index:04d}" for index in range(len(vectors))], dtype="S")
    h5_path = tmp_path / "protein_embeddings.h5"
    with h5py.File(h5_path, "w") as handle:
        handle.create_dataset("protein_ids", data=protein_ids)
        handle.create_dataset("embeddings", data=vectors)

    inspection = FAISSStore.build_persistent_index(
        h5_path,
        chunk_size=211,
        nprobe=32,
    )
    assert inspection.index_kind == "ivf_flat_ip"

    store = FAISSStore(h5_path)
    try:
        result = store.query("protein-0000", k=1, include_distances=True)
    finally:
        store.close()

    assert result[0][0] == "protein-0001"


def test_query_vector_does_not_mutate_caller_array(tmp_path):
    h5_path = tmp_path / "protein_embeddings.h5"
    with h5py.File(h5_path, "w") as handle:
        handle.create_dataset("protein_ids", data=np.asarray(["a", "b"], dtype="S"))
        handle.create_dataset(
            "embeddings",
            data=np.asarray([[1.0, 0.0], [0.0, 1.0]], dtype=np.float32),
        )
    store = FAISSStore(h5_path)
    query = np.asarray([2.0, 0.0], dtype=np.float32)
    original = query.copy()
    try:
        store.query_vector(query, k=1)
    finally:
        store.close()

    np.testing.assert_array_equal(query, original)


@pytest.mark.parametrize(
    ("bad_vector", "message"),
    [
        (np.asarray([[0.0, 0.0]], dtype=np.float32), "zero vectors"),
        (np.asarray([[np.nan, 1.0]], dtype=np.float32), "NaN or infinite"),
    ],
)
def test_index_build_rejects_invalid_vectors(tmp_path, bad_vector, message):
    h5_path = tmp_path / "protein_embeddings.h5"
    with h5py.File(h5_path, "w") as handle:
        handle.create_dataset("protein_ids", data=np.asarray(["bad"], dtype="S"))
        handle.create_dataset("embeddings", data=bad_vector)

    with pytest.raises(ValueError, match=message):
        FAISSStore.build_persistent_index(h5_path)


def test_torch_free_index_runner_emits_machine_readable_result(tmp_path):
    h5_path = tmp_path / "protein_embeddings.h5"
    with h5py.File(h5_path, "w") as handle:
        handle.create_dataset("protein_ids", data=np.asarray(["a", "b"], dtype="S"))
        handle.create_dataset(
            "embeddings",
            data=np.asarray([[1.0, 0.0], [0.0, 1.0]], dtype=np.float32),
        )

    process = subprocess.run(
        [
            sys.executable,
            "-m",
            "sharur.ingest.vector_index_runner",
            "--embeddings",
            str(h5_path),
            "--threads",
            "1",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    payload = json.loads(process.stdout)

    assert payload["state"] == "available"
    assert Path(payload["index_path"]).is_file()
