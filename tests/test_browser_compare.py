"""Browser comparison views: genomes or clades side by side, and the function × clade heatmap."""

import re

import pytest
from fastapi.testclient import TestClient

from sharur.browser import create_app
from sharur.browser import routes_compare as rc
from sharur.predicates_v2.persistence import generate_and_persist_v2
from sharur.storage.duckdb_store import DuckDBStore


SEQ = "MKTAYIAKQRQISFVKSHFSRQ" * 10
# three genomes in class Alpha carry ABC transporters; three in Beta carry short-chain dehydrogenases
GENOMES = {f"a{i}": ("Alpha", "PF00005", "ABC_tran", "ABC transporter") for i in range(3)} | {
    f"b{i}": ("Beta", "PF00106", "adh_short", "short chain dehydrogenase") for i in range(3)}


@pytest.fixture
def db(tmp_path):
    path = tmp_path / "sharur.duckdb"
    store = DuckDBStore(str(path))
    bins, contigs, proteins, annotations = [], [], [], []
    for n, (bin_id, (cls, acc, name, desc)) in enumerate(GENOMES.items()):
        bins.append(f"('{bin_id}', 'd__Bacteria;p__Testota;c__{cls};o__{cls}ales;f__F;g__G{n}', 2000)")
        contigs.append(f"('{bin_id}_c1', '{bin_id}', 9000, 'assembly')")
        for k in range(2):
            pid = f"{bin_id}_c1_{k + 1}"
            proteins.append(f"('{pid}', '{bin_id}_c1', '{bin_id}', {1 + 1000 * k}, {700 + 1000 * k}, '+', {k}, '{SEQ}', 220)")
        annotations.append(f"({n + 1}, '{bin_id}_c1_1', 'pfam', '{acc}', '{name}', '{desc}', 1e-40, 150.0, 10, 160)")
    store.conn.execute(f"""
        INSERT INTO bins (bin_id, taxonomy, total_length) VALUES {", ".join(bins)};
        INSERT INTO contigs (contig_id, bin_id, length, length_source) VALUES {", ".join(contigs)};
        INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand, gene_index, sequence,
                              sequence_length) VALUES {", ".join(proteins)};
        INSERT INTO annotations (annotation_id, protein_id, source, accession, name, description, evalue, score,
                                 start_aa, end_aa) VALUES {", ".join(annotations)};
        INSERT INTO defense_systems (system_id, genome_id, system_type, system_subtype, genes_count, protein_ids)
        VALUES ('sys1', 'a0', 'Viperin', 'Viperin', 1, 'a0_c1_2'), ('sys2', 'a1', 'Viperin', 'Viperin', 1, 'a1_c1_2');
    """)
    generate_and_persist_v2(store, chunk_size=2, return_states=False, update_legacy_predicates=True)
    store.close()
    return path


@pytest.fixture
def client(db):
    return TestClient(create_app(db, background=False))


def test_sides_resolve_from_genomes_and_clades(client):
    catalog = client.app.state.catalog
    url = lambda *parts: "/" + "/".join(parts)  # noqa: E731
    genome = rc.resolve_side(catalog, "a0", url)
    clade = rc.resolve_side(catalog, "class:Beta", url)
    assert genome.is_genome and [g.bin_id for g in genome.genomes] == ["a0"]
    assert clade.kind == "class" and {g.bin_id for g in clade.genomes} == {"b0", "b1", "b2"}
    assert rc.resolve_side(catalog, "class:Nope", url) is None and rc.resolve_side(catalog, "zzz", url) is None
    assert genome.filter_param == "genome=a0" and clade.filter_param == "clade=class%3ABeta"


def test_clade_comparison_reports_differences_with_links(client):
    page = client.get("/compare", params={"a": "class:Alpha", "b": "class:Beta"})
    assert page.status_code == 200
    text = page.text
    assert "Alpha" in text and "Beta" in text and "swap" in text
    # ABC transporter only in Alpha, linked to Alpha's proteins; dehydrogenase only in Beta
    assert re.search(r'href="/function/abc_transporter\?clade=class%3AAlpha"', text)
    assert re.search(r'href="/function/[a-z_]*dehydrogenase[a-z_]*\?clade=class%3ABeta"', text)
    assert 'href="/domain/PF00005"' in text and 'href="/domain/PF00106"' in text
    assert 'href="/system/defense/Viperin"' in text  # 2 of 3 Alpha genomes vs none in Beta
    assert "MKTAYIAKQ" not in text


def test_genome_comparisons_and_mixed_sides(client):
    pair = client.get("/compare", params={"a": "a0", "b": "b0"})
    assert pair.status_code == 200 and "not detected" in pair.text and "proteins" in pair.text
    assert 'href="/function/abc_transporter?genome=a0"' in pair.text
    mixed = client.get("/compare", params={"a": "a0", "b": "class:Beta"})
    assert mixed.status_code == 200 and "medians across genomes" in mixed.text
    assert client.get("/compare", params={"a": "a0", "b": "nope"}).status_code == 404
    picker = client.get("/compare", params={"a": "a0"})
    assert picker.status_code == 200 and 'name="a" value="a0"' in picker.text


def test_compare_entry_points(client):
    assert 'name="a" value="a0"' in client.get("/genome/a0").text
    assert 'name="a" value="class:Alpha"' in client.get("/taxa/class/Alpha").text
    root = client.get("/taxa").text
    assert 'href="/compare"' in root and 'href="/heatmap"' in root
    assert 'href="/heatmap"' in client.get("/functions").text


def test_function_pages_filter_to_a_clade_or_genome(client):
    whole = client.get("/function/abc_transporter")
    alpha = client.get("/function/abc_transporter", params={"clade": "class:Alpha"})
    one = client.get("/function/abc_transporter", params={"genome": "a1"})
    assert whole.status_code == alpha.status_code == one.status_code == 200
    assert "show the whole dataset" in alpha.text and "Alpha (class)" in alpha.text
    assert "a1_c1_1" in one.text and "a0_c1_1" not in one.text
    assert "show the whole dataset" not in whole.text


def test_heatmap_cells_presets_and_small_clades(client):
    page = client.get("/heatmap", params={"functions": "abc_transporter", "rank": "class"})
    assert page.status_code == 200
    cells = re.findall(r'href="/function/abc_transporter\?clade=class%3A(\w+)" style="--v:([0-9.]+)"', page.text)
    assert dict(cells) == {"Alpha": "1.000", "Beta": "0.000"}
    # genera have one genome each: hidden by default, shown on request
    genus = client.get("/heatmap", params={"functions": "abc_transporter", "rank": "genus"})
    assert "6 clades with fewer than 3 genomes hidden" in genus.text and "clade=genus" not in genus.text
    shown = client.get("/heatmap", params={"functions": "abc_transporter", "rank": "genus", "small": 1})
    assert shown.text.count("clade=genus%3AG") == 6
    assert client.get("/heatmap").status_code == 200  # default preset or empty state
    data = rc.heatmap(client.app.state.catalog, ["abc_transporter", "not_a_label"], "class", False)
    assert [r["carriers"] for r in data["rows"]] == [3, 0]


def test_suggestions_filter_by_kind(client):
    taxa = client.get("/api/suggest", params={"q": "alph", "kind": "taxon,genome"}).json()
    assert taxa and all(r["kind"] in ("taxon", "genome") for r in taxa)
    functions = client.get("/api/suggest", params={"q": "abc", "kind": "function"}).json()
    assert functions and all(r["kind"] == "function" for r in functions)
    assert client.get("/api/suggest", params={"q": "abc"}).status_code == 200


def test_presets_and_descendants_follow_the_vocabulary():
    hydrogenases = rc.descendants("hydrogenase")
    assert hydrogenases[0] == "hydrogenase" and "nife_group3" in hydrogenases and "fefe_groupA" in hydrogenases
    assert rc.descendants("not_a_predicate") == []
