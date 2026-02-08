import uuid

from bennu.core.session import ExplorationSession


def test_working_set_create_and_resolve():
    session = ExplorationSession()
    ws = session.create_set("foo", ["p1", "p2"])
    assert ws.name == "foo"
    assert len(ws.protein_ids) == 2
    resolved = session.resolve_reference("foo")
    assert resolved == ws.set_id.hex


def test_focus_pronoun_resolution():
    session = ExplorationSession()
    session.push_focus("protein", "gene_001")
    resolved = session.resolve_reference("it")
    assert resolved == "gene_001"


def test_provenance_logging():
    session = ExplorationSession()
    entry = session.log_query("q", [], "summary", 10)
    assert entry.query == "q"
    assert len(session.get_provenance()) == 1


def test_tool_log_updates_provenance():
    session = ExplorationSession()
    from bennu.tools.find_proteins import FindProteinsTool, FindProteinsParams

    tool = FindProteinsTool()
    result = tool.execute(FindProteinsParams(limit=1), session)
    assert result.success
    assert len(session.get_provenance()) == 1


def test_multiple_tools_append_provenance():
    session = ExplorationSession()
    # Seed minimal protein for context tool
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
        "('p1','c1','b1',10,50,'+','MAAA')"
    )

    from bennu.tools.get_context import GetContextTool, GetContextParams
    from bennu.tools.find_proteins import FindProteinsTool, FindProteinsParams

    fp = FindProteinsTool()
    fp.execute(FindProteinsParams(limit=1), session)
    gc = GetContextTool()
    gc.execute(GetContextParams(protein_id="p1", window_genes=1), session)

    assert len(session.get_provenance()) == 2
