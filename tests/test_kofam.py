"""KOfam score_type: which bitscore a KO row stores, and the record of it."""

import duckdb
import pandas as pd

from sharur.kofam import (
    SCORE_BASIS_FULL,
    SCORE_BASIS_KOFAM,
    KoThreshold,
    compares_as_kofam,
    ko_score_bases,
    kofam_scores,
    read_ko_list,
    record_ko_score_basis,
)


KO_LIST = (
    "knum\tthreshold\tscore_type\tprofile_type\tdefinition\n"
    "K00001\t100.0\tfull\tall\tfull-type KO\n"
    "K00002\t200.0\tdomain\tall\tdomain-type KO\n"
    "K00003\t-\t-\t-\tno threshold\n"
)


def _domains(*rows):
    """Astra/Aksha hits rows: (protein, KO, full-sequence bitscore, domain bitscore, ali_from)."""
    return pd.DataFrame(
        [{"protein_id": p, "accession": ko, "bitscore": full, "dom_bitscore": dom, "ali_from": start,
          "evalue": 1e-50} for p, ko, full, dom, start in rows]
    )


def test_read_ko_list(tmp_path):
    path = tmp_path / "ko_list"
    path.write_text(KO_LIST)
    assert read_ko_list(path) == {
        "K00001": KoThreshold(100.0, "full"),
        "K00002": KoThreshold(200.0, "domain"),
        "K00003": KoThreshold(None, None),
    }
    assert set(read_ko_list(path, ["K00002"])) == {"K00002"}


def test_kofam_scores_follow_score_type():
    df = _domains(
        ("p1", "K00001", 150.0, 90.0, 10),    # full KO: full-sequence score, though no domain reaches 100
        ("p1", "K00001", 150.0, 60.0, 300),
        ("p1", "K00002", 260.0, 180.0, 400),  # domain KO: best domain, though the full score reaches 200
        ("p1", "K00002", 260.0, 120.0, 20),
        ("p2", "K00003", 80.0, 79.0, 5),      # no threshold: full-sequence score
        ("p2", "K09999", 70.0, 75.0, 5),      # not in ko_list: full-sequence score
    )
    score_types = {"K00001": "full", "K00002": "domain", "K00003": None}
    out, basis = kofam_scores(df, score_types)
    assert basis == SCORE_BASIS_KOFAM
    got = {(r.protein_id, r.accession): (r.score, r.ali_from) for r in out.itertuples()}
    # One row per (protein, KO), at the best domain's coordinates.
    assert got == {("p1", "K00001"): (150.0, 10), ("p1", "K00002"): (180.0, 400),
                   ("p2", "K00003"): (80.0, 5), ("p2", "K09999"): (70.0, 5)}

    without_ko_list, basis = kofam_scores(df, None)
    assert basis == SCORE_BASIS_FULL
    assert {(r.protein_id, r.accession): r.score for r in without_ko_list.itertuples()} == {
        ("p1", "K00001"): 150.0, ("p1", "K00002"): 260.0, ("p2", "K00003"): 80.0, ("p2", "K09999"): 70.0}


def test_compares_as_kofam():
    table = {(basis, score_type): compares_as_kofam(basis, score_type)
             for basis in (SCORE_BASIS_KOFAM, SCORE_BASIS_FULL, None) for score_type in ("full", "domain", None)}
    assert table == {
        (SCORE_BASIS_KOFAM, "full"): True, (SCORE_BASIS_KOFAM, "domain"): True, (SCORE_BASIS_KOFAM, None): True,
        (SCORE_BASIS_FULL, "full"): True, (SCORE_BASIS_FULL, "domain"): False, (SCORE_BASIS_FULL, None): True,
        (None, "full"): False, (None, "domain"): False, (None, None): False,
    }


def test_score_basis_record_round_trip(tmp_path):
    path = tmp_path / "ko_list"
    path.write_text(KO_LIST)
    conn = duckdb.connect()
    assert ko_score_bases(conn) == {}
    record_ko_score_basis(conn, "kegg", SCORE_BASIS_KOFAM, ko_list=path, hits_paths=["KOFAM_hits_df.tsv"],
                          raw_rows=6, loaded_rows=4)
    assert ko_score_bases(conn) == {"kegg": SCORE_BASIS_KOFAM}
    row = conn.execute("SELECT model_set, hits_path, raw_rows, loaded_rows FROM annotation_source_versions").fetchone()
    assert row == ("KOfam", "KOFAM_hits_df.tsv", 6, 4)
    record_ko_score_basis(conn, "kegg", None)  # basis unknown: no record
    assert ko_score_bases(conn) == {}
