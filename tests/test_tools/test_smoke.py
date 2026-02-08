from bennu.core.session import ExplorationSession
from bennu.tools.find_proteins import FindProteinsTool, FindProteinsParams
from bennu.tools.manage_sets import ManageSetsTool, ManageSetsParams
from bennu.tools.find_similar import FindSimilarTool, FindSimilarParams


def test_manage_sets_create_and_list():
    session = ExplorationSession()
    tool = ManageSetsTool()
    res = tool.execute(
        ManageSetsParams(action="create", set_name="s1", protein_ids=["p1", "p2"]), session
    )
    assert res.success
    assert session.get_set("s1") is not None

    res_list = tool.execute(ManageSetsParams(action="list"), session)
    assert res_list.count == 1


def test_find_proteins_runs():
    session = ExplorationSession()
    tool = FindProteinsTool()
    res = tool.execute(FindProteinsParams(limit=10), session)
    assert res.success
    assert res.count == 0


def test_find_proteins_exclude_set():
    session = ExplorationSession()
    # Seed proteins
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
        "('p_include','c1','b1',10,50,'+','MAAA'),"
        "('p_exclude','c1','b1',60,90,'-','MBBB')"
    )
    session.create_set("all", ["p_include", "p_exclude"])
    session.create_set("exclude", ["p_exclude"])

    tool = FindProteinsTool()
    params = FindProteinsParams(limit=10, in_set="all", exclude_set="exclude")
    res = tool.execute(params, session)
    assert res.success
    ids = {p.protein_id for p in res.data}
    assert ids == {"p_include"}


def test_find_similar_defaults_without_vector_store():
    session = ExplorationSession()
    tool = FindSimilarTool()
    res = tool.execute(FindSimilarParams(query_id="foo"), session)
    assert res.success
    assert res.count == 0  # No vector store configured
