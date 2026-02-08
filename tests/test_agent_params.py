import types

import dspy
import pytest

from bennu.agent.orchestrator import BennuAgent
from bennu.core.session import ExplorationSession


class DummySigObj:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


class DummyPredictor:
    def __init__(self, sig_cls):
        self.sig_cls = sig_cls

    def __call__(self, **kwargs):
        # Return window_genes and taxonomy to verify mapping
        if self.sig_cls.__name__ == "FindProteinsParams":
            return DummySigObj(domains=["PF00001"], taxonomy="d__Bacteria", limit=5)
        if self.sig_cls.__name__ == "GetContextParams":
            return DummySigObj(protein_id="p1", window_genes=3, window_bp=0, include_annotations=True, expand_operon=False)
        if self.sig_cls.__name__ == "FindAnomaliesParams":
            return DummySigObj(scope="set", scope_id="weird", signals=["fusion"], min_score=0.5, limit=7)
        if self.sig_cls.__name__ == "DetectLociParams":
            return DummySigObj(locus_type="cas", scope="contig", scope_id="ctg123")
        if self.sig_cls.__name__ == "ManageSetsParams":
            return DummySigObj(action="create", set_name="demo", protein_ids=["pA"], description="desc", other_set="aux")
        if self.sig_cls.__name__ == "ExportParams":
            return DummySigObj(format="json", source="set", source_id="export_me", include_sequences=False, include_annotations=False)
        return DummySigObj()


def test_param_extraction_populates_fields(monkeypatch):
    session = ExplorationSession()
    agent = BennuAgent(session)
    agent.use_llm = True  # force TypedPredictor path
    agent._predictor_factory = lambda sig_cls: DummyPredictor(sig_cls)

    params = agent._extract_params("find_proteins", "find stuff")
    assert params.domains == ["PF00001"]
    assert params.taxonomy == "d__Bacteria"
    assert params.limit == 5

    ctx_params = agent._extract_params("get_context", "what's around it")
    assert ctx_params.protein_id == "p1"
    assert ctx_params.window_genes == 3


def test_param_extraction_maps_additional_tools(monkeypatch):
    session = ExplorationSession()
    agent = BennuAgent(session)
    agent.use_llm = True  # force TypedPredictor path
    agent._predictor_factory = lambda sig_cls: DummyPredictor(sig_cls)

    an_params = agent._extract_params("find_anomalies", "look for anomalies in that set")
    assert an_params.scope == "set"
    assert an_params.scope_id == "weird"
    assert an_params.signals == ["fusion"]
    assert an_params.limit == 7

    loci_params = agent._extract_params("detect_loci", "find cas loci on contig")
    assert loci_params.locus_type == "cas"
    assert loci_params.scope == "contig"
    assert loci_params.scope_id == "ctg123"

    ms_params = agent._extract_params("manage_sets", "create a demo set")
    assert ms_params.action == "create"
    assert ms_params.set_name == "demo"
    assert ms_params.protein_ids == ["pA"]

    ex_params = agent._extract_params("export", "export that set to json")
    assert ex_params.format == "json"
    assert ex_params.source == "set"
    assert ex_params.source_id == "export_me"
    assert ex_params.include_sequences is False


def test_param_extraction_uses_focus_when_missing():
    session = ExplorationSession()
    session.push_focus("protein", "gene_focus")
    agent = BennuAgent(session)
    agent.use_llm = False  # force fallback path

    ctx_params = agent._extract_params("get_context", "context please")
    assert ctx_params.protein_id == "gene_focus"

    sim_params = agent._extract_params("find_similar", "similar proteins")
    assert sim_params.query_id == "gene_focus"


def test_param_extraction_populates_in_set_from_working_sets():
    session = ExplorationSession()
    session.create_set("placeholder_set", ["x1", "x2"])
    agent = BennuAgent(session)
    agent.use_llm = False

    params = agent._extract_params("find_proteins", "find stuff in my set")
    assert params.in_set == "placeholder_set"

    from bennu.tools.find_proteins import FindProteinsTool

    res = FindProteinsTool().execute(params, session)
    assert res.success


def test_manage_sets_and_export_fill_from_focus_and_sets():
    session = ExplorationSession()
    session.create_set("placeholder_set", ["x1", "x2"])
    session.push_focus("protein", "gene_focus")
    agent = BennuAgent(session)
    agent.use_llm = False

    ms_params = agent._extract_params("manage_sets", "add this protein")
    assert ms_params.protein_ids == ["gene_focus"]
    assert ms_params.set_name == "placeholder_set"

    session.push_focus("set", "placeholder_set")
    ex_params = agent._extract_params("export", "export it")
    assert ex_params.source == "set"
    assert ex_params.source_id == "placeholder_set"


def test_scope_fills_from_focus_bin():
    session = ExplorationSession()
    session.push_focus("bin", "bin123")
    agent = BennuAgent(session)
    agent.use_llm = False

    det_params = agent._extract_params("detect_loci", "find loci")
    assert det_params.scope_id == "bin123"
    assert det_params.scope == "genome" or det_params.scope in {"bin", "genome"}

    an_params = agent._extract_params("find_anomalies", "any anomalies")
    assert an_params.scope in {"genome", "bin"}
