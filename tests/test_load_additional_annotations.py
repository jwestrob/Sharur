"""scripts/load_additional_annotations.py: KO rows carry the score KOfam's score_type names."""

from importlib.machinery import SourceFileLoader
from pathlib import Path

import duckdb

from sharur.kofam import SCORE_BASIS_KOFAM, ko_score_bases
from sharur.storage.schema import SCHEMA


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "load_additional_annotations.py"
KO_LIST = (
    "knum\tthreshold\tscore_type\tprofile_type\tdefinition\n"
    "K02588\t362.67\tfull\tall\tnitrogenase iron protein NifH\n"
    "K02591\t403.20\tdomain\tall\tnitrogenase molybdenum-iron protein beta chain\n"
)
HEADER = "sequence_id\thmm_name\tbitscore\tevalue\tc_evalue\ti_evalue\tenv_from\tenv_to\tdom_bitscore\n"


def _db(tmp_path, proteins):
    path = tmp_path / "sharur.duckdb"
    conn = duckdb.connect(str(path))
    conn.execute(SCHEMA)
    conn.execute("INSERT INTO bins (bin_id) VALUES ('b')")
    conn.execute("INSERT INTO contigs (contig_id, bin_id, length) VALUES ('c', 'b', 9000)")
    for i, pid in enumerate(proteins):
        conn.execute("INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand, gene_index) "
                     "VALUES (?, 'c', 'b', ?, ?, '+', ?)", [pid, i * 1000 + 1, i * 1000 + 900, i])
    conn.close()
    return path


def test_ko_rows_store_kofam_scores_and_record_the_basis(tmp_path, monkeypatch):
    module = SourceFileLoader("load_additional_annotations", str(SCRIPT)).load_module()
    ko_list = tmp_path / "ko_list"
    ko_list.write_text(KO_LIST)
    monkeypatch.setenv("SHARUR_KOFAM_KO_LIST", str(ko_list))  # found as sharur.health.kofam_ko_list() would
    tsv = tmp_path / "KOFAM_hits_df.tsv"
    tsv.write_text(HEADER
                   + "pH\tK02588\t400.0\t1e-100\t1e-100\t1e-100\t1\t150\t300.0\n"     # full KO, two domains
                   + "pH\tK02588\t400.0\t1e-100\t1e-100\t1e-100\t160\t280\t150.0\n"
                   + "pK\tK02591\t500.0\t1e-100\t1e-100\t1e-100\t1\t200\t390.0\n"     # domain KO, two domains
                   + "pK\tK02591\t500.0\t1e-100\t1e-100\t1e-100\t220\t400\t300.0\n")
    db = _db(tmp_path, ["pH", "pK"])
    module.load_annotations(str(db), str(tsv), "kegg")

    conn = duckdb.connect(str(db), read_only=True)
    rows = conn.execute("SELECT protein_id, score, start_aa FROM annotations ORDER BY 1").fetchall()
    bases = ko_score_bases(conn)
    conn.close()
    # One row per (protein, KO) at the best domain; nifH keeps its full-sequence score,
    # nifK (score_type domain) its best domain's.
    assert rows == [("pH", 400.0, 1), ("pK", 390.0, 1)]
    assert bases == {"kegg": SCORE_BASIS_KOFAM}


def test_ko_rows_appended_to_unrecorded_rows_leave_the_basis_unknown(tmp_path):
    module = SourceFileLoader("load_additional_annotations", str(SCRIPT)).load_module()
    ko_list = tmp_path / "ko_list"
    ko_list.write_text(KO_LIST)
    db = _db(tmp_path, ["pH", "pK"])
    conn = duckdb.connect(str(db))
    conn.execute("INSERT INTO annotations (annotation_id, protein_id, source, accession, score, evalue) "
                 "VALUES (1, 'pK', 'kegg', 'K02591', 410.0, 1e-100)")  # loaded earlier, basis unrecorded
    conn.close()
    tsv = tmp_path / "KOFAM_hits_df.tsv"
    tsv.write_text(HEADER + "pH\tK02588\t400.0\t1e-100\t1e-100\t1e-100\t1\t150\t300.0\n")
    module.load_annotations(str(db), str(tsv), "kegg", ko_list=str(ko_list))

    conn = duckdb.connect(str(db), read_only=True)
    assert ko_score_bases(conn) == {}
    conn.close()
