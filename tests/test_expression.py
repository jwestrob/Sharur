"""Metatranscriptome gene-count import, genome transcript share, top genes and feature expression."""

import math

import pytest

from sharur.expression import (
    feature_expression,
    feature_expression_markdown,
    genome_expression,
    genome_genes,
    genome_genes_markdown,
    import_counts,
    import_gene_flags,
    samples,
)
from sharur.storage.duckdb_store import DuckDBStore


def _core() -> DuckDBStore:
    store = DuckDBStore()
    store.conn.execute("""
        INSERT INTO bins (bin_id) VALUES ('A'), ('B');
        INSERT INTO contigs (contig_id, bin_id, length) VALUES ('a1', 'A', 10000), ('b1', 'B', 10000);
        INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand, gene_index, sequence_length)
        VALUES ('a1_1', 'a1', 'A', 1, 900, '+', 0, 300), ('a1_2', 'a1', 'A', 1001, 1300, '-', 1, 100),
               ('a1_3', 'a1', 'A', 2001, 2600, '+', 2, 200), ('b1_1', 'b1', 'B', 1, 900, '+', 0, 300);
        INSERT INTO annotations (annotation_id, protein_id, source, accession, name, evalue)
        VALUES (1, 'b1_1', 'kegg', 'K00370', 'narG', 1e-50), (2, 'a1_2', 'kegg', 'K00399', 'mcrA', 1e-80);
        INSERT INTO protein_predicates (protein_id, predicates) VALUES ('a1_2', ['mcr']), ('b1_1', ['nitrate_reductase']);
    """)
    return store


LONG = """sample_id\tprotein_id\tcount_s1\tcount_s2
r1\ta1_1\t2\t90
r1\ta1_2\t1\t300
r1\ta1_3\t0\t0
r1\tb1_1\t1\t30
r2\ta1_1\t40\t0
r2\tb1_1\t160\t4
"""


@pytest.fixture
def loaded(tmp_path):
    store = _core()
    path = tmp_path / "counts.tsv"
    path.write_text(LONG)
    sidecar = tmp_path / "expression.duckdb"
    result = import_counts(store, sidecar, [path])
    return store, sidecar, result


def test_orientation_is_the_strand_holding_more_fragments(loaded):
    _store, sidecar, result = loaded
    assert result["sense_orientation"] == {"r1": "s2", "r2": "s1"}
    assert result["rows"] == 5  # the all-zero row is dropped
    meta = {s["sample_id"]: s for s in samples(sidecar)}
    assert meta["r1"]["sense_fragments"] == 420
    assert math.isclose(meta["r1"]["strand_specificity"], 420 / 424)
    assert math.isclose(meta["r2"]["strand_specificity"], 200 / 204)


def test_gene_vote_resists_a_one_strand_hotspot(tmp_path):
    store = DuckDBStore()
    rows = ", ".join(f"('g{i}', 'c1', 'A', {i * 1000 + 1}, {i * 1000 + 300}, '+', {i}, 100)" for i in range(26))
    store.conn.execute(f"""
        INSERT INTO bins (bin_id) VALUES ('A');
        INSERT INTO contigs (contig_id, bin_id, length) VALUES ('c1', 'A', 30000);
        INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand, gene_index, sequence_length)
        VALUES {rows};""")
    # 25 mRNAs read on s2; one ORF inside an abundant non-coding RNA collects 100,000 fragments on s1
    lines = ["sample_id\tprotein_id\tcount_s1\tcount_s2", "r\tg0\t100000\t0"]
    lines += [f"r\tg{i}\t2\t50" for i in range(1, 26)]
    path = tmp_path / "counts.tsv"
    path.write_text("\n".join(lines) + "\n")
    result = import_counts(store, tmp_path / "e.duckdb", [path])
    assert result["sense_orientation"] == {"r": "s2"}
    meta = samples(tmp_path / "e.duckdb")[0]
    assert math.isclose(meta["strand_specificity"], 25 / 26)
    assert meta["sense_fragments"] == 25 * 50


def test_metadata_orientation_overrides_auto(tmp_path):
    store = _core()
    path = tmp_path / "counts.tsv"
    path.write_text(LONG)
    meta = tmp_path / "samples.tsv"
    meta.write_text("sample_id\tsense_orientation\tdepth_cm\nr1\ts1\t40\n")
    result = import_counts(store, tmp_path / "e.duckdb", [path], sample_metadata=meta)
    assert result["sense_orientation"]["r1"] == "s1"
    assert {s["sample_id"]: s for s in samples(tmp_path / "e.duckdb")}["r1"]["metadata"] == {
        "sense_orientation": "s1", "depth_cm": "40"}


def test_unknown_proteins_are_rejected_unless_allowed(tmp_path):
    store = _core()
    path = tmp_path / "counts.tsv"
    path.write_text(LONG + "r1\tzz_9\t5\t5\n")
    with pytest.raises(ValueError, match="absent from the dataset"):
        import_counts(store, tmp_path / "e.duckdb", [path])
    result = import_counts(store, tmp_path / "e2.duckdb", [path], allow_unknown_proteins=True)
    assert result["unknown_proteins"] == 1
    assert result["rows"] == 5


def test_featurecounts_tables_import_with_given_orientation(tmp_path):
    store = _core()
    path = tmp_path / "fc.txt"
    path.write_text("# Program:featureCounts v2.0.6\n"
                    "Geneid\tChr\tStart\tEnd\tStrand\tLength\t/x/libA.bam\n"
                    "a1_1\ta1\t1\t900\t+\t900\t7\n"
                    "a1_2\ta1\t1001\t1300\t-\t300\t0\n"
                    "b1_1\tb1\t1\t900\t+\t900\t3\n")
    result = import_counts(store, tmp_path / "e.duckdb", [path], fmt="featurecounts")
    assert result["samples"] == ["libA"]
    assert result["rows"] == 2
    rows = genome_expression(tmp_path / "e.duckdb")
    assert [(r["bin_id"], r["fragments"]) for r in rows] == [("A", 7), ("B", 3)]


def test_reimport_replaces_the_whole_sample(tmp_path):
    store = _core()
    sidecar = tmp_path / "e.duckdb"
    first, recount = tmp_path / "v1.tsv", tmp_path / "v2.tsv"
    first.write_text("sample_id\tprotein_id\tsense\tantisense\nr\ta1_1\t50\t1\nr\ta1_2\t300\t2\nq\ta1_2\t9\t0\n")
    recount.write_text("sample_id\tprotein_id\tsense\tantisense\nr\ta1_1\t48\t1\n")  # a1_2 lost its reads
    import_counts(store, sidecar, [first])
    result = import_counts(store, sidecar, [recount])
    assert result["cleared_rows"] == 1
    meta = {s["sample_id"]: s for s in samples(sidecar)}
    assert meta["r"]["sense_fragments"] == 48
    assert meta["q"]["sense_fragments"] == 9  # samples outside the import keep their counts
    assert [(r["sample_id"], r["fragments"]) for r in genome_expression(sidecar)] == [("q", 9), ("r", 48)]


def test_featurecounts_orientations_merge_across_imports(tmp_path):
    store = _core()
    sidecar = tmp_path / "e.duckdb"
    header = "# Program:featureCounts v2.0.6\nGeneid\tChr\tStart\tEnd\tStrand\tLength\tlibA.bam\n"
    sense, anti = tmp_path / "s.txt", tmp_path / "a.txt"
    sense.write_text(header + "a1_1\ta1\t1\t900\t+\t900\t50\na1_2\ta1\t1001\t1300\t-\t300\t20\n")
    anti.write_text(header + "a1_1\ta1\t1\t900\t+\t900\t3\nb1_1\tb1\t1\t900\t+\t900\t4\n")
    import_counts(store, sidecar, [sense], fmt="featurecounts", orientation="sense")
    import_counts(store, sidecar, [anti], fmt="featurecounts", orientation="antisense")
    (meta,) = samples(sidecar)
    assert meta["sense_orientation"] == "sense"
    assert meta["sense_fragments"] == 70
    assert math.isclose(meta["strand_specificity"], 70 / 77)  # fewer than 20 voting genes: fragment totals
    assert {(r["bin_id"], r["fragments"]) for r in genome_expression(sidecar)} == {("A", 70), ("B", 0)}


def test_genome_transcript_share_sums_to_one(loaded):
    _, sidecar, _ = loaded
    rows = genome_expression(sidecar)
    for sample in ("r1", "r2"):
        assert math.isclose(sum(r["transcript_share"] for r in rows if r["sample_id"] == sample), 1.0)
    r1 = {r["bin_id"]: r for r in rows if r["sample_id"] == "r1"}
    assert r1["A"]["fragments"] == 390
    assert r1["A"]["genes_with_reads"] == 2


def test_genome_genes_ranks_by_length_normalized_rate(loaded):
    store, sidecar, _ = loaded
    result = genome_genes(store, sidecar, "A", sample_ids=["r1"], top=3)
    (s,) = result["samples"]
    assert s["genes"] == 3
    assert s["genes_with_reads"] == 2
    top = s["top"]
    assert top[0]["protein_id"] == "a1_2"
    assert top[0]["percentile"] == 1.0
    assert "mcrA" in top[0]["label"]
    assert math.isclose(sum(g["genome_tpm"] for g in top), 1e6)
    assert "a1_2" in genome_genes_markdown(result)


def test_feature_expression_share_and_within_genome_rank(loaded):
    store, sidecar, _ = loaded
    result = feature_expression(store, sidecar, annotation="K00399", min_genome_fragments=10)
    shares = {s["sample_id"]: s for s in result["samples"]}
    assert shares["r1"]["fragments"] == 300
    assert math.isclose(shares["r1"]["share_of_sample"], 300 / 420)
    assert shares["r2"]["fragments"] == 0
    genome = result["genomes"][0]
    assert genome["bin_id"] == "A"
    assert genome["samples"] == 2
    assert "K00399" in feature_expression_markdown(result)
    with pytest.raises(ValueError, match="exactly one"):
        feature_expression(store, sidecar)


def test_genes_without_reads_rank_at_zero(loaded):
    store, sidecar, _ = loaded
    (s,) = genome_genes(store, sidecar, "A", sample_ids=["r2"], top=3)["samples"]
    pct = {g["protein_id"]: g["percentile"] for g in s["top"]}
    assert pct["a1_1"] == 1.0
    assert pct["a1_2"] == 0.0 and pct["a1_3"] == 0.0  # two of three genes have no reads
    mcr = feature_expression(store, sidecar, annotation="K00399", min_genome_fragments=10)
    by_sample = {g["bin_id"]: g for g in mcr["genomes"]}
    assert by_sample["A"]["median_best_percentile"] == 0.5  # 1.0 in r1, 0.0 in r2


def test_flagged_genes_leave_every_query(loaded, tmp_path):
    store, sidecar, _ = loaded
    flags = tmp_path / "overlap.tsv"
    flags.write_text("protein_id\trfam_acc\trna_name\toverlap_bp\toverlap_fraction\tsame_strand\n"
                     "a1_2\tRF00177\tSSU_rRNA_bacteria\t300\t1.0\ttrue\n"
                     "a1_1\tRF00005\ttRNA\t20\t0.02\tfalse\n")
    result = import_gene_flags(sidecar, flags)
    assert result["flagged_proteins"] == 1  # a1_1 overlaps by 20 bp, below the 60 bp default
    r1 = {r["bin_id"]: r for r in genome_expression(sidecar, sample_ids=["r1"])}
    assert r1["A"]["fragments"] == 90
    assert math.isclose(r1["A"]["transcript_share"], 90 / 120)
    assert genome_expression(sidecar, sample_ids=["r1"], include_flagged=True)[0]["fragments"] == 390
    (s,) = genome_genes(store, sidecar, "A", sample_ids=["r1"])["samples"]
    assert s["genes"] == 2
    assert all(g["protein_id"] != "a1_2" for g in s["top"])
    mcr = feature_expression(store, sidecar, annotation="K00399", min_genome_fragments=10)
    assert {x["sample_id"]: x["fragments"] for x in mcr["samples"]} == {"r1": 0, "r2": 0}
    meta = {m["sample_id"]: m for m in samples(sidecar)}
    assert meta["r1"]["flagged_sense_fragments"] == 300
