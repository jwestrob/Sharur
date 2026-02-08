from bennu.core.session import ExplorationSession
from bennu.tools.export import ExportTool, ExportParams


def test_export_tsv_no_sequences():
    session = ExplorationSession()
    # Seed minimal rows
    session.db.execute(
        "INSERT INTO bins (bin_id, completeness, contamination, taxonomy, n_contigs, total_length) VALUES "
        "('b1', 90, 1, 'd__Bacteria', 1, 1000)"
    )
    session.db.execute(
        "INSERT INTO contigs (contig_id, bin_id, length, gc_content, is_circular, taxonomy) VALUES "
        "('c1','b1',1000,0.5,FALSE,'d__Bacteria')"
    )
    session.db.execute(
        "INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand, sequence) VALUES "
        "('p1','c1','b1',10,50,'+',''),"
        "('p2','c1','b1',60,90,'-','')"
    )
    session.create_set("s1", ["p1", "p2"])

    tool = ExportTool()
    res = tool.execute(ExportParams(format="tsv", source="set", source_id="s1", include_sequences=False), session)
    assert res.success
    assert "p1" in res.data
    assert res.summary.startswith("Exported 2 proteins")


def test_export_json_includes_annotations():
    session = ExplorationSession()
    session.db.execute(
        "INSERT INTO bins (bin_id, completeness, contamination, taxonomy, n_contigs, total_length) VALUES "
        "('b1', 90, 1, 'd__Bacteria', 1, 1000)"
    )
    session.db.execute(
        "INSERT INTO contigs (contig_id, bin_id, length, gc_content, is_circular, taxonomy) VALUES "
        "('c1','b1',1000,0.5,FALSE,'d__Bacteria')"
    )
    session.db.execute(
        "INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand, sequence) VALUES "
        "('p1','c1','b1',10,50,'+','MAAA'),"
        "('p2','c1','b1',60,90,'-','MBBB')"
    )
    session.db.execute(
        "INSERT INTO annotations (annotation_id, protein_id, source, accession) VALUES "
        "(1,'p1','pfam','PF00001')"
    )
    session.create_set("s1", ["p1", "p2"])

    tool = ExportTool()
    res = tool.execute(ExportParams(format="json", source="set", source_id="s1"), session)
    assert res.success
    assert '"accession": "PF00001"' in res.data


def test_export_fasta_includes_sequences():
    session = ExplorationSession()
    session.db.execute(
        "INSERT INTO bins (bin_id, completeness, contamination, taxonomy, n_contigs, total_length) VALUES "
        "('b1', 90, 1, 'd__Bacteria', 1, 1000)"
    )
    session.db.execute(
        "INSERT INTO contigs (contig_id, bin_id, length, gc_content, is_circular, taxonomy) VALUES "
        "('c1','b1',1000,0.5,FALSE,'d__Bacteria')"
    )
    session.db.execute(
        "INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand, sequence) VALUES "
        "('p1','c1','b1',10,50,'+','MAAA'),"
        "('p2','c1','b1',60,90,'-','MBBB')"
    )
    session.create_set("s1", ["p1", "p2"])

    tool = ExportTool()
    res = tool.execute(ExportParams(format="fasta", source="set", source_id="s1"), session)
    assert res.success
    assert ">p1" in res.data and "MAAA" in res.data
