"""Browser startup summaries replayed from a per-database cache match the live computation."""

import os
import threading

import numpy as np
import pytest

from sharur.browser import startup_cache
from sharur.browser.catalog import load_background, load_catalog
from sharur.browser.routes_compare import PfamPairs
from sharur.predicates_v2.persistence import generate_and_persist_v2
from sharur.storage.duckdb_store import DuckDBStore

SEQ = "MKTAYIAKQRQISFVKSHFSRQ" * 10


@pytest.fixture
def db(tmp_path):
    path = tmp_path / "sharur.duckdb"
    store = DuckDBStore(str(path))
    c = store.conn
    for b in ("g1", "g2"):
        c.execute("INSERT INTO bins (bin_id, completeness, contamination, taxonomy) VALUES (?, 90, 1, 'd__Archaea')", [b])
        c.execute("INSERT INTO contigs (contig_id, bin_id, length, length_source) VALUES (?, ?, 5000, 'assembly')",
                  [f"{b}_c", b])
        c.executemany("INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand, gene_index, "
                      "sequence, sequence_length) VALUES (?, ?, ?, ?, ?, '+', ?, ?, 220)",
                      [(f"{b}_{i}", f"{b}_c", b, 1000 * i + 1, 1000 * i + 600, i, SEQ) for i in range(3)])
    c.executemany("INSERT INTO annotations (annotation_id, protein_id, source, accession, name, evalue, score) "
                  "VALUES (?, ?, ?, ?, ?, 1e-30, 100)",
                  [(1, "g1_0", "pfam", "PF00005.1", "ABC_tran"), (2, "g2_0", "pfam", "PF00005.1", "ABC_tran"),
                   (3, "g2_1", "pfam", "PF00004.2", "AAA"), (4, "g1_1", "kofam", "K02003", "K02003"),
                   (5, "g2_1", "kofam", "K02004", "K02004")])
    generate_and_persist_v2(store, chunk_size=10, return_states=False, update_legacy_predicates=True)
    store.close()
    return path


def _summaries(source):
    lock = threading.Lock()
    catalog = load_catalog(source)
    load_background(source, catalog, lock)
    pfam = PfamPairs()
    pfam.load(source, catalog, lock)
    return catalog, pfam


def test_replayed_summaries_match_live_ones_and_go_stale_with_the_database(db, tmp_path):
    cache_dir = tmp_path / "cache"
    startup_cache.build(db, startup_cache.cache_file(db, cache_dir))
    rows, artifacts = startup_cache.load(db, cache_dir)
    store = DuckDBStore(str(db), read_only=True)
    try:
        live, live_pfam = _summaries(store)
        replayer = startup_cache.Replayer(store, rows, artifacts)
        cached, cached_pfam = _summaries(replayer)
    finally:
        store.close()

    assert replayer.misses == 0
    assert cached.totals == live.totals and cached.predicates == live.predicates
    assert cached.ko_sets == live.ko_sets == {"g1": {"K02003"}, "g2": {"K02004"}}
    assert cached.category_share == live.category_share
    order = np.lexsort((live.pair_pred, live.pair_bin)), np.lexsort((cached.pair_pred, cached.pair_bin))
    for a in ("pair_bin", "pair_pred", "pair_count"):
        assert np.array_equal(getattr(live, a)[order[0]], getattr(cached, a)[order[1]])
    assert cached_pfam.accessions == live_pfam.accessions == ["PF00004", "PF00005"]
    assert sorted(zip(cached_pfam.bin_idx, cached_pfam.acc_idx)) == sorted(zip(live_pfam.bin_idx, live_pfam.acc_idx))

    st = db.stat()
    os.utime(db, ns=(st.st_atime_ns, st.st_mtime_ns + 1_000_000_000))   # the database changed
    assert startup_cache.load(db, cache_dir) is None
