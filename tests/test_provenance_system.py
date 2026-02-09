"""Tests for provenance DAG, evidence linking, hypothesis registry, and renderer."""

import json
import uuid
from pathlib import Path

import pytest

from bennu.core.types import Evidence, Hypothesis, HypothesisStatus, ProvenanceEntry
from bennu.core.session import ExplorationSession
from bennu.core.hypothesis_registry import HypothesisRegistry
from bennu.core.provenance_renderer import (
    _sanitize_label,
    render_provenance_mermaid,
    render_provenance_summary,
)


# ------------------------------------------------------------------ #
# Change 1: ProvenanceEntry parent_ids
# ------------------------------------------------------------------ #


class TestProvenanceParentIds:
    def test_default_empty(self):
        entry = ProvenanceEntry(
            query="q", tool_calls=[], results_summary="r", duration_ms=10
        )
        assert entry.parent_ids == []

    def test_with_parent_ids(self):
        pid = uuid.uuid4()
        entry = ProvenanceEntry(
            query="q",
            tool_calls=[],
            results_summary="r",
            duration_ms=10,
            parent_ids=[pid],
        )
        assert entry.parent_ids == [pid]

    def test_json_round_trip(self):
        pid = uuid.uuid4()
        entry = ProvenanceEntry(
            query="q",
            tool_calls=[],
            results_summary="r",
            duration_ms=10,
            parent_ids=[pid],
        )
        data = json.loads(entry.model_dump_json())
        restored = ProvenanceEntry(**data)
        assert restored.parent_ids == [pid]

    def test_log_query_without_parent_ids(self):
        session = ExplorationSession()
        entry = session.log_query("q", [], "summary", 10)
        assert entry.parent_ids == []

    def test_log_query_with_parent_ids(self):
        session = ExplorationSession()
        e1 = session.log_query("q1", [], "s1", 10)
        e2 = session.log_query("q2", [], "s2", 20, parent_ids=[e1.entry_id])
        assert e2.parent_ids == [e1.entry_id]

    def test_session_save_load_preserves_parent_ids(self, tmp_path):
        session = ExplorationSession()
        e1 = session.log_query("q1", [], "s1", 10)
        e2 = session.log_query("q2", [], "s2", 20, parent_ids=[e1.entry_id])

        path = tmp_path / "session.json"
        session.save(path)
        loaded = ExplorationSession.load(path)

        prov = loaded.get_provenance()
        assert len(prov) == 2
        assert prov[0].parent_ids == []
        assert prov[1].parent_ids == [e1.entry_id]


# ------------------------------------------------------------------ #
# Change 2: Evidence provenance_id linking
# ------------------------------------------------------------------ #


class TestEvidenceProvenanceLink:
    def test_evidence_default_no_provenance_id(self):
        ev = Evidence(query="q", result_summary="r", supports=True, confidence=0.8)
        assert ev.provenance_id is None

    def test_evidence_with_provenance_id(self):
        pid = uuid.uuid4()
        ev = Evidence(
            query="q", result_summary="r", supports=True, confidence=0.8,
            provenance_id=pid,
        )
        assert ev.provenance_id == pid

    def test_evidence_json_round_trip(self):
        pid = uuid.uuid4()
        ev = Evidence(
            query="q", result_summary="r", supports=True, confidence=0.8,
            provenance_id=pid,
        )
        data = json.loads(ev.model_dump_json())
        restored = Evidence(**data)
        assert restored.provenance_id == pid

    def test_add_evidence_from_provenance(self):
        session = ExplorationSession()
        hypo = session.propose_hypothesis("Test hypothesis")
        entry = session.log_query("SELECT count(*)", [], "42 rows", 100)

        session.add_evidence_from_provenance(
            hypo.hypothesis_id, entry.entry_id, True, 0.9
        )

        evidence = session.list_hypotheses()[0].evidence
        assert len(evidence) == 1
        assert evidence[0].provenance_id == entry.entry_id
        assert evidence[0].query == "SELECT count(*)"
        assert evidence[0].result_summary == "42 rows"
        assert evidence[0].supports is True
        assert evidence[0].confidence == 0.9

    def test_add_evidence_from_provenance_unknown_entry(self):
        session = ExplorationSession()
        hypo = session.propose_hypothesis("Test")
        with pytest.raises(KeyError, match="Provenance entry"):
            session.add_evidence_from_provenance(
                hypo.hypothesis_id, uuid.uuid4(), True, 0.8
            )

    def test_add_evidence_from_provenance_unknown_hypothesis(self):
        session = ExplorationSession()
        entry = session.log_query("q", [], "r", 10)
        with pytest.raises(KeyError, match="Hypothesis"):
            session.add_evidence_from_provenance(
                uuid.uuid4(), entry.entry_id, True, 0.8
            )

    def test_session_save_load_preserves_provenance_id(self, tmp_path):
        session = ExplorationSession()
        hypo = session.propose_hypothesis("Test")
        entry = session.log_query("q", [], "r", 10)
        session.add_evidence_from_provenance(
            hypo.hypothesis_id, entry.entry_id, True, 0.7
        )

        path = tmp_path / "session.json"
        session.save(path)
        loaded = ExplorationSession.load(path)

        evidence = loaded.list_hypotheses()[0].evidence
        assert evidence[0].provenance_id == entry.entry_id

    def test_old_add_evidence_still_works(self):
        """Backward compat: original add_evidence() keeps working."""
        session = ExplorationSession()
        hypo = session.propose_hypothesis("Test")
        session.add_evidence(hypo.hypothesis_id, "q", "r", True, 0.8)
        ev = session.list_hypotheses()[0].evidence[0]
        assert ev.provenance_id is None


# ------------------------------------------------------------------ #
# Change 3: HypothesisRegistry
# ------------------------------------------------------------------ #


class TestHypothesisRegistry:
    def test_creates_file_on_save(self, tmp_path):
        path = tmp_path / "exploration" / "hypotheses.json"
        reg = HypothesisRegistry(path)
        h = Hypothesis(statement="Test")
        reg.register(h)
        assert path.exists()

    def test_load_existing(self, tmp_path):
        path = tmp_path / "hypotheses.json"
        reg = HypothesisRegistry(path)
        h = Hypothesis(statement="Test")
        reg.register(h)

        reg2 = HypothesisRegistry(path)
        assert len(reg2.list_all()) == 1
        assert reg2.get(h.hypothesis_id).statement == "Test"

    def test_register_persists(self, tmp_path):
        path = tmp_path / "hypotheses.json"
        reg = HypothesisRegistry(path)
        h = Hypothesis(statement="Test")
        reg.register(h)

        data = json.loads(path.read_text())
        assert len(data) == 1

    def test_add_evidence_persists(self, tmp_path):
        path = tmp_path / "hypotheses.json"
        reg = HypothesisRegistry(path)
        h = Hypothesis(statement="Test")
        reg.register(h)

        ev = Evidence(query="q", result_summary="r", supports=True, confidence=0.9)
        reg.add_evidence(h.hypothesis_id, ev)

        reg2 = HypothesisRegistry(path)
        assert len(reg2.get(h.hypothesis_id).evidence) == 1

    def test_update_status_persists(self, tmp_path):
        path = tmp_path / "hypotheses.json"
        reg = HypothesisRegistry(path)
        h = Hypothesis(statement="Test")
        reg.register(h)
        reg.update_status(h.hypothesis_id, HypothesisStatus.SUPPORTED)

        reg2 = HypothesisRegistry(path)
        assert reg2.get(h.hypothesis_id).status == HypothesisStatus.SUPPORTED

    def test_list_active(self, tmp_path):
        path = tmp_path / "hypotheses.json"
        reg = HypothesisRegistry(path)
        h1 = Hypothesis(statement="Active")
        h2 = Hypothesis(statement="Refuted", status=HypothesisStatus.REFUTED)
        h3 = Hypothesis(statement="Supported", status=HypothesisStatus.SUPPORTED)
        reg.register(h1)
        reg.register(h2)
        reg.register(h3)

        active = reg.list_active()
        assert len(active) == 2
        statements = {h.statement for h in active}
        assert statements == {"Active", "Supported"}

    def test_summary(self, tmp_path):
        path = tmp_path / "hypotheses.json"
        reg = HypothesisRegistry(path)
        assert reg.summary() == "No hypotheses registered."

        h = Hypothesis(statement="Test hypothesis")
        reg.register(h)
        ev = Evidence(query="q", result_summary="r", supports=True, confidence=0.8)
        reg.add_evidence(h.hypothesis_id, ev)

        summary = reg.summary()
        assert "Test hypothesis" in summary
        assert "+1/-0" in summary

    def test_get_unknown_raises(self, tmp_path):
        path = tmp_path / "hypotheses.json"
        reg = HypothesisRegistry(path)
        with pytest.raises(KeyError, match="not found"):
            reg.get(uuid.uuid4())

    def test_round_trip_with_evidence_provenance_id(self, tmp_path):
        path = tmp_path / "hypotheses.json"
        reg = HypothesisRegistry(path)
        h = Hypothesis(statement="Test")
        reg.register(h)

        prov_id = uuid.uuid4()
        ev = Evidence(
            query="q", result_summary="r", supports=True, confidence=0.8,
            provenance_id=prov_id,
        )
        reg.add_evidence(h.hypothesis_id, ev)

        reg2 = HypothesisRegistry(path)
        loaded_ev = reg2.get(h.hypothesis_id).evidence[0]
        assert loaded_ev.provenance_id == prov_id


# ------------------------------------------------------------------ #
# Change 4: Provenance renderer
# ------------------------------------------------------------------ #


class TestProvenanceRenderer:
    def test_sanitize_label_basic(self):
        assert _sanitize_label("hello") == "hello"

    def test_sanitize_label_truncates(self):
        result = _sanitize_label("a" * 50, max_len=10)
        assert len(result) == 10
        assert result.endswith("...")

    def test_sanitize_label_escapes_quotes(self):
        assert '"' not in _sanitize_label('say "hello"')

    def test_render_empty_session(self):
        session = ExplorationSession()
        result = render_provenance_mermaid(session)
        assert result.startswith("graph TD")

    def test_render_with_entries(self):
        session = ExplorationSession()
        e1 = session.log_query("Count proteins", [], "100", 10)
        e2 = session.log_query("Filter large", [], "20", 20, parent_ids=[e1.entry_id])

        result = render_provenance_mermaid(session)
        assert "Count proteins" in result
        assert "Filter large" in result
        assert "-->" in result  # parent edge

    def test_render_with_hypothesis(self):
        session = ExplorationSession()
        hypo = session.propose_hypothesis("Test hypothesis")
        entry = session.log_query("q", [], "r", 10)
        session.add_evidence_from_provenance(
            hypo.hypothesis_id, entry.entry_id, True, 0.9
        )

        result = render_provenance_mermaid(session)
        assert "hypothesis" in result  # classDef
        assert "supports" in result

    def test_render_refuting_evidence(self):
        session = ExplorationSession()
        hypo = session.propose_hypothesis("Test")
        entry = session.log_query("q", [], "r", 10)
        session.add_evidence_from_provenance(
            hypo.hypothesis_id, entry.entry_id, False, 0.6
        )

        result = render_provenance_mermaid(session)
        assert "against" in result
        assert "-.->" in result

    def test_render_with_title(self):
        session = ExplorationSession()
        result = render_provenance_mermaid(session, title="My Analysis")
        assert "My Analysis" in result
        assert "title_node" in result

    def test_render_with_error_node(self):
        session = ExplorationSession()
        session.log_query("bad query", [], "", 10, error="failed")
        result = render_provenance_mermaid(session)
        assert "errorNode" in result

    def test_render_with_registry(self):
        session = ExplorationSession()
        session.log_query("q", [], "r", 10)

        from bennu.core.hypothesis_registry import HypothesisRegistry
        import tempfile

        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
            path = Path(f.name)

        try:
            reg = HypothesisRegistry(path)
            h = Hypothesis(statement="Registry hypothesis")
            reg.register(h)

            result = render_provenance_mermaid(session, registry=reg)
            assert "Registry hypothesis" in result
        finally:
            path.unlink(missing_ok=True)

    def test_summary_empty_session(self):
        session = ExplorationSession()
        result = render_provenance_summary(session)
        assert result == "No provenance entries."

    def test_summary_with_dag(self):
        session = ExplorationSession()
        e1 = session.log_query("q1", [], "r1", 10)
        e2 = session.log_query("q2", [], "r2", 20, parent_ids=[e1.entry_id])
        e3 = session.log_query("q3", [], "r3", 30, parent_ids=[e1.entry_id, e2.entry_id])

        result = render_provenance_summary(session)
        assert "3 entries" in result
        assert "Root entries: 1" in result
        assert "Leaf entries: 1" in result
        assert "Max chain depth: 2" in result

    def test_summary_reports_linked_hypotheses(self):
        session = ExplorationSession()
        hypo = session.propose_hypothesis("Test")
        entry = session.log_query("q", [], "r", 10)
        session.add_evidence_from_provenance(
            hypo.hypothesis_id, entry.entry_id, True, 0.9
        )

        result = render_provenance_summary(session)
        assert "1 with provenance links" in result
