"""Source-ID ambiguity and peptide boundaries must survive historical remapping."""
import importlib.util
from pathlib import Path
import pytest


def mapper():
    path = Path(__file__).resolve().parents[2] / "analyses/scripts/remap_faa_to_prodigal.py"
    if not path.is_file():
        pytest.skip("Historical analysis script is available in the local research checkout")
    spec = importlib.util.spec_from_file_location("historical_coordinate_remap", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def gene(contig, start, sequence="MPEPTIDE", strand="+"):
    return {"contig_id": contig, "start": start, "end": start + 26,
            "strand": strand, "sequence": sequence, "gene_index": 0}


def test_unique_exact_peptide_preserves_locus_and_strand():
    updates, matched, unmatched = mapper().remap_genome(
        "genome", [gene("contig", 100, strand="-")], [("protein", "MPEPTIDE*")])
    assert (matched, unmatched) == (1, 0)
    assert updates == [{"protein_id": "protein", "contig_id": "contig", "gene_index": 0,
                        "start": 100, "end": 126, "strand": "-"}]


def test_identical_peptides_at_distinct_loci_remain_unresolved():
    module = mapper()
    for candidates in ([gene("a", 100), gene("b", 500)],
                       [gene("b", 500), gene("a", 100)]):
        assert module.remap_genome("genome", candidates, [("source", "MPEPTIDE")]) == ([], 0, 1)


def test_distinct_source_ids_sharing_one_peptide_remain_unresolved():
    assert mapper().remap_genome("genome", [gene("a", 100)],
                               [("source1", "MPEPTIDE"), ("source2", "MPEPTIDE*")]) == ([], 0, 2)


def test_initial_methionine_removal_cannot_reuse_shorter_cds_bounds():
    assert mapper().remap_genome("genome", [gene("a", 100, sequence="PEPTIDE")],
                               [("source", "MPEPTIDE")]) == ([], 0, 1)


def test_empty_source_payload_remains_unresolved():
    assert mapper().remap_genome("genome", [gene("a", 100, sequence="")],
                               [("source", None)]) == ([], 0, 1)


def test_parsed_prodigal_strand_uses_schema_symbols():
    parsed = mapper()._parse_header_and_seq("contig_1 # 100 # 126 # -1 # ID=1_1;partial=00", "MPEPTIDE*")
    assert parsed["strand"] == "-"


def test_legacy_writer_refuses_modern_protein_columns():
    import duckdb
    c = duckdb.connect()
    c.execute("CREATE TABLE proteins(protein_id VARCHAR PRIMARY KEY, partial VARCHAR)")
    with pytest.raises(RuntimeError, match="before mutation.*protein columns"):
        mapper().validate_legacy_schema(c)
    assert c.execute("SELECT count(*) FROM proteins").fetchone()[0] == 0


def test_legacy_writer_refuses_unpreserved_dependencies():
    import duckdb
    c = duckdb.connect()
    c.execute("CREATE TABLE proteins(protein_id VARCHAR PRIMARY KEY, contig_id VARCHAR, "
              "bin_id VARCHAR, gene_index INTEGER, start INTEGER, end_coord INTEGER, "
              "strand VARCHAR, sequence VARCHAR, sequence_length INTEGER, gc_content FLOAT)")
    c.execute("CREATE TABLE system_proteins(protein_id VARCHAR REFERENCES proteins(protein_id))")
    with pytest.raises(RuntimeError, match="before mutation.*system_proteins"):
        mapper().validate_legacy_schema(c)
