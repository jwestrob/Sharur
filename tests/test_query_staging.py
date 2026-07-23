"""Sealed, atomic campaign-local DuckDB replica staging."""

from __future__ import annotations

import json
import stat

import pytest

from sharur.dataset_seal import (
    DatasetSealError,
    build_dataset_seal,
    write_dataset_seal,
)
from sharur.query.staging import stage_database
from sharur.storage.duckdb_store import DuckDBStore


def _sealed_database(root):
    database = root / "sharur.duckdb"
    store = DuckDBStore(database)
    store.execute(
        "INSERT INTO bins(bin_id, n_contigs, total_length) VALUES ('bin1', 1, 300)"
    )
    store.execute(
        "INSERT INTO contigs(contig_id, bin_id, length) VALUES ('contig1', 'bin1', 300)"
    )
    store.execute(
        """
        INSERT INTO proteins(
            protein_id, contig_id, bin_id, start, end_coord, strand,
            sequence_length
        ) VALUES ('p1', 'contig1', 'bin1', 1, 300, '+', 100)
        """
    )
    store.close()
    seal_path = root / "dataset.seal.json"
    seal = build_dataset_seal(database, hash_large_files=True)
    write_dataset_seal(seal, seal_path)
    return database, seal_path, seal


def test_stage_database_publishes_verified_read_only_replica_and_reuses_it(tmp_path):
    source_root = tmp_path / "source"
    source_root.mkdir()
    database, seal_path, seal = _sealed_database(source_root)
    stage_dir = tmp_path / "campaign-local"

    first = stage_database(
        database,
        stage_dir,
        seal_path=seal_path,
        reserve_bytes=0,
    )

    assert first.path.is_file()
    assert first.path.parent == stage_dir.resolve()
    assert first.dataset_id == seal["dataset_id"]
    assert first.reused is False
    assert stat.S_IMODE(first.path.stat().st_mode) == 0o444
    metadata = json.loads(
        first.path.with_suffix(".duckdb.stage.json").read_text()
    )
    assert metadata["dataset_id"] == seal["dataset_id"]
    reader = DuckDBStore(first.path, read_only=True)
    assert reader.execute("SELECT bin_id FROM bins") == [("bin1",)]
    reader.close()

    second = stage_database(
        database,
        stage_dir,
        seal_path=seal_path,
        reserve_bytes=0,
    )
    assert second.path == first.path
    assert second.reused is True
    assert second.staged_at == first.staged_at


def test_stage_database_repairs_a_corrupt_cached_replica(tmp_path):
    source_root = tmp_path / "source"
    source_root.mkdir()
    database, seal_path, _seal = _sealed_database(source_root)
    stage_dir = tmp_path / "campaign-local"
    staged = stage_database(
        database,
        stage_dir,
        seal_path=seal_path,
        reserve_bytes=0,
    )
    staged.path.chmod(0o644)
    staged.path.write_bytes(b"corrupt")

    repaired = stage_database(
        database,
        stage_dir,
        seal_path=seal_path,
        reserve_bytes=0,
    )

    assert repaired.reused is False
    reader = DuckDBStore(repaired.path, read_only=True)
    assert reader.execute("SELECT COUNT(*) FROM bins") == [(1,)]
    reader.close()


def test_stage_database_fails_closed_for_unsealed_or_drifted_sources(tmp_path):
    unsealed_root = tmp_path / "unsealed"
    unsealed_root.mkdir()
    unsealed = unsealed_root / "sharur.duckdb"
    store = DuckDBStore(unsealed)
    _ = store.conn
    store.close()
    with pytest.raises(DatasetSealError, match="require"):
        stage_database(unsealed, tmp_path / "stage-a", reserve_bytes=0)

    sealed_root = tmp_path / "sealed"
    sealed_root.mkdir()
    database, seal_path, _seal = _sealed_database(sealed_root)
    writer = DuckDBStore(database)
    writer.execute("INSERT INTO bins(bin_id) VALUES ('drift')")
    writer.close()
    with pytest.raises(DatasetSealError, match="verification failed"):
        stage_database(
            database,
            tmp_path / "stage-b",
            seal_path=seal_path,
            reserve_bytes=0,
        )
