"""Protein-page annotation display: KO names, links, threshold wording, lanes."""

import duckdb

from sharur.browser import annotations as ann


def _store(tmp_path):
    conn = duckdb.connect(str(tmp_path / "a.duckdb"))
    conn.execute("CREATE TABLE annotations (protein_id VARCHAR, source VARCHAR, accession VARCHAR, name VARCHAR, "
                 "description VARCHAR, evalue DOUBLE, start_aa INTEGER, end_aa INTEGER)")
    conn.execute("""INSERT INTO annotations VALUES
        ('p', 'kegg', 'K03737', 'pfor', 'pyruvate-ferredoxin oxidoreductase [EC:1.2.7.1 1.2.7.-]', 1e-30, 5, 300),
        ('p', 'kofam', 'K03737', 'K03737', 'evalue_1e-15', 1e-20, NULL, NULL),
        ('p', 'pfam', 'PF01558.21', 'POR', 'Pyruvate ferredoxin oxidoreductase', 1e-10, 10, 200),
        ('p', 'defensefinder', 'X__Y', 'X__Y', NULL, 1e-9, 20, 80)""")
    return conn


def test_ko_names_links_and_threshold_text(tmp_path):
    names = ann.KoNames(_store(tmp_path), kegg_dir=tmp_path)
    view = ann.hit_view("kofam", {"accession": "K03737", "name": "K03737", "description": "evalue_1e-15"}, names)
    assert view["external"] == "https://www.kegg.jp/entry/K03737"
    assert "pyruvate-ferredoxin" in view["name"] and 'href="https://www.kegg.jp/entry/ec:1.2.7.1"' in view["name"]
    assert "1.2.7.-" in view["name"] and "ec:1.2.7.-" not in view["name"]
    assert view["detail"] == "E ≤ 1e-15 (KO has no KOfam threshold)"
    assert ann.threshold_text("GA") == "passes KOfam score threshold"
    assert ann.external_url("pfam", "PF01558.21") == "https://www.ebi.ac.uk/interpro/entry/pfam/PF01558/"
    assert ann.describe_ko("kofam:K03737", names).startswith("kofam:K03737 (pyruvate-ferredoxin oxidoreductase")


def test_snapshot_symbols_override(tmp_path):
    (tmp_path / "kegg_evidence_snapshot.tsv").write_text("# built\n# ko\tsymbols\tname\tbrite\tmodules\n"
                                                         "K03737\tpor, nifJ\tPFOR [EC:1.2.7.1]\t\t\n")
    names = ann.KoNames(_store(tmp_path), kegg_dir=tmp_path)
    assert names.get("K03737") == ("por, nifJ", "PFOR [EC:1.2.7.1]")
    assert ann.describe_ko("kofam:K03737", names) == "kofam:K03737 por (PFOR)"


def test_domain_lanes_one_per_placed_source(tmp_path):
    lanes = ann.domain_lanes(_store(tmp_path), "p")
    assert [s for s, _ in lanes] == ["pfam", "defensefinder", "kegg"]
    assert lanes[0][1][0]["start_aa"] == 10
