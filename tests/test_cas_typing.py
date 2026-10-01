"""CRISPRCasTyper port: thresholds, operons, subtype scoring, linking, storage."""

import json

import duckdb
import pytest

from sharur import cas_typing as ct

TYPES = ["I-A", "I-B", "III-B", "V-A"]
SCORING = {
    "Cas1_0_I": {"I-A": 2, "I-B": 2},
    "Cas2_0_I": {"I-A": 1, "I-B": 1},
    "Cas3_0_IA": {"I-A": 4},
    "Cas5_0_IB": {"I-B": 4},
    "Cas7_0_IB": {"I-B": 3},
    "Cas8b_0_IB": {"I-B": 4},
    "Cas10_0_IIIB": {"III-B": 4},
    "Cmr4_0_IIIB": {"III-B": 3},
    "Cmr5_0_IIIB": {"III-B": 3},
    "Cas12a_0_VA": {"V-A": 4},
}


@pytest.fixture
def data(tmp_path):
    lines = ["Hmm," + ",".join(TYPES)]
    for hmm, row in SCORING.items():
        lines.append(hmm + "," + ",".join(str(row.get(t, "")) for t in TYPES))
    (tmp_path / "CasScoring.csv").write_text("\n".join(lines) + "\n")
    (tmp_path / "cutoffs.tab").write_text("cas12a:1e-5,0.9,0.6\n")
    (tmp_path / "interference.json").write_text(json.dumps(
        {"I-B": [["Cas5"], ["Cas7"], ["Cas8"]], "III-B": [["Cas10"], ["Cmr4"], ["Cmr5"]], "V-A": [["Cas12"]]}))
    (tmp_path / "adaptation.json").write_text(json.dumps({"I-B": [["Cas1_"], ["Cas2"]]}))
    (tmp_path / ct.PROFILES).write_text("")
    return ct.TypingData.load(tmp_path)


def hit(hmm, pos, score=100.0, strand=1, pid=None, evalue=1e-20, cov=0.95):
    return {"hmm": hmm, "protein_id": pid or f"p{pos}", "pos": pos, "score": score, "evalue": evalue,
            "cov_seq": cov, "cov_hmm": cov, "start": pos * 1000, "end": pos * 1000 + 900, "strand": strand}


def test_signatures_come_from_specific_cutoffs(data):
    assert data.signature == {"Cas12a"}
    assert data.single_effector == {"V-A"}


def test_best_profile_is_chosen_before_thresholds(data):
    hits = [hit("Cas5_0_IB", 1, score=50, cov=0.1), hit("Cas7_0_IB", 1, score=40),
            hit("Cas12a_0_VA", 2, cov=0.8), hit("Cas1_0_I", 3, evalue=0.5)]
    kept = ct.filter_hits(hits, data)
    # p1's best profile fails coverage (its weaker passing profile is ignored); the Cas12a
    # hit fails its own sequence-coverage cutoff (0.8 < 0.9); p3 fails the overall E-value
    assert kept == []


def test_cluster_operons_allows_three_genes_between():
    assert ct.cluster_operons([1, 2, 6, 11, 20], dist=3) == [[1, 2, 6], [11], [20]]


def test_type_operon_scores_and_completeness(data):
    rows = [hit("Cas1_0_I", 1), hit("Cas2_0_I", 2), hit("Cas5_0_IB", 3), hit("Cas7_0_IB", 4, strand=-1)]
    op = ct.type_operon(rows, data)
    assert op["prediction"] == "I-B" and op["best_score"] == 10
    assert op["complete_interference"] == "67%" and op["complete_adaptation"] == "100%"
    assert op["strand_interference"] == 0 and op["strand_adaptation"] == 1


def test_low_score_operon_without_signature_is_false(data):
    rows = [hit("Cas1_0_I", 1), hit("Cas2_0_I", 2), hit("Cas1_0_I", 3, pid="dup")]
    assert ct.type_operon(rows, data)["prediction"] == "False"


def test_single_signature_gene(data):
    assert ct.type_operon([hit("Cas12a_0_VA", 5)], data)["prediction"] == "V-A"
    assert ct.type_operon([hit("Cas10_0_IIIB", 5)], data)["prediction"] == "False"


def test_adjacent_systems_make_a_hybrid(data):
    rows = [hit(h, i) for i, h in enumerate(["Cas3_0_IA", "Cas5_0_IB", "Cas7_0_IB", "Cas8b_0_IB",
                                              "Cas10_0_IIIB", "Cmr4_0_IIIB", "Cmr5_0_IIIB"], 1)]
    op = ct.type_operon(rows, data)
    assert op["prediction"] == "Hybrid(I-B,III-B)"
    assert op["best_type"] == ["I-B", "III-B"]


def test_reconcile_operon_and_array_calls():
    assert ct.reconcile("I-B", "I-B", "I-B") == "I-B"
    assert ct.reconcile("Ambiguous", ["I-B", "I-C"], "I-C") == "I-C"
    assert ct.reconcile("Ambiguous", ["I-B", "I-C"], "I-E") == "I"
    assert ct.reconcile("False", "I-B", "I-B") == "I-B(Putative)"
    assert ct.reconcile("False", "I-B", "III-A") == "Unknown"
    assert ct.reconcile("III-B", "III-B", "Unknown") == "III-B"


def test_link_arrays_within_distance():
    ops = [{"contig_id": "c1", "start": 1000, "end": 5000, "prediction": "I-B", "best_type": "I-B"},
           {"contig_id": "c1", "start": 90000, "end": 95000, "prediction": "False", "best_type": "I-A"}]
    arrays = [{"locus_id": "a1", "contig_id": "c1", "start": 12000, "end": 13000, "prediction": "I-B",
               "trusted": False},
              {"locus_id": "a2", "contig_id": "c1", "start": 96000, "end": 97000, "prediction": "I-A",
               "trusted": False}]
    ct.link_arrays(ops, arrays)
    assert ops[0]["status"] == "crispr_cas" and ops[0]["distances"] == [7000]
    assert arrays[0]["trusted"] and arrays[0]["near_cas"]
    assert ops[1]["status"] == "cas_putative" and not arrays[1]["trusted"]


def test_repeat_features_and_identity():
    feats = ct.repeat_features("ACGTACGT")
    assert len(ct.canonical_kmers()) == 136 and feats["Length"] == 8 and feats["GC"] == 0.5
    assert feats["ACGT"] == 2  # palindrome counted on its canonical strand
    assert ct._identity("ACGT", "AGT") == 75.0
    assert ct.consensus_repeat(["AAA", "CCC", "CCC"]) == "CCC"


def test_write_results_replaces_earlier_calls(tmp_path):
    from sharur.storage.schema import SCHEMA, SCHEMA_VERSION

    path = tmp_path / "d.duckdb"
    conn = duckdb.connect(str(path))
    conn.execute(SCHEMA)
    conn.execute("CREATE TABLE IF NOT EXISTS schema_version (version INTEGER, description VARCHAR, "
                 "applied_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP)")
    conn.execute("INSERT INTO schema_version (version, description) VALUES (?, 'test')", [SCHEMA_VERSION])
    conn.execute("INSERT INTO bins (bin_id) VALUES ('g')")
    conn.execute("INSERT INTO contigs (contig_id, bin_id, length) VALUES ('c', 'g', 5000)")
    conn.execute("INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand) "
                 "VALUES ('p1', 'c', 'g', 1000, 1900, '+')")
    conn.close()
    op = {"genome_id": "g", "contig_id": "c", "start": 1, "end": 900, "status": "cas", "joint_prediction": "V-A",
          "prediction": "V-A", "best_type": "V-A", "best_score": 4.0, "complete_interference": "100%",
          "complete_adaptation": "NA", "strand_interference": 1, "strand_adaptation": "NA",
          "genes": [hit("Cas12a_0_VA", 1)], "arrays": [], "distances": []}
    result = {"systems": [op], "arrays": []}
    assert ct.write_results(path, result) == {"systems": 1, "arrays": 0, "members": 1}
    assert ct.write_results(path, result)["systems"] == 1
    conn = duckdb.connect(str(path), read_only=True)
    assert conn.execute("SELECT system_id, prediction FROM crispr_cas_systems").fetchall() == [
        ("cctyper:g:c:1-900", "V-A")]
    assert conn.execute("SELECT COUNT(*) FROM system_proteins WHERE system_source = 'cctyper'").fetchone()[0] == 1
