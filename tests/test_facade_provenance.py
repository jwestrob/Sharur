"""Tests for provenance/hypothesis facade methods on Bennu class."""

import uuid
from pathlib import Path

import pytest

from bennu.operators import Bennu


@pytest.fixture
def bennu_instance(tmp_path):
    """Create a Bennu instance with a temporary DuckDB database."""
    db_path = tmp_path / "dataset" / "bennu.duckdb"
    db_path.parent.mkdir(parents=True, exist_ok=True)
    b = Bennu(db_path)
    # Force session creation so the DB is initialized
    _ = b.store
    return b


# ------------------------------------------------------------------ #
# hypothesis_registry property
# ------------------------------------------------------------------ #


class TestHypothesisRegistryProperty:
    def test_lazy_creation(self, bennu_instance):
        b = bennu_instance
        assert b._hypothesis_registry is None
        reg = b.hypothesis_registry
        assert reg is not None
        assert b._hypothesis_registry is reg

    def test_path_is_correct(self, bennu_instance):
        b = bennu_instance
        expected = b._db_path.parent / "exploration" / "hypotheses.json"
        assert b.hypothesis_registry.path == expected


# ------------------------------------------------------------------ #
# propose_hypothesis
# ------------------------------------------------------------------ #


class TestProposeHypothesis:
    def test_returns_hypothesis(self, bennu_instance):
        b = bennu_instance
        h = b.propose_hypothesis("Test statement")
        assert h.statement == "Test statement"
        assert h.hypothesis_id is not None

    def test_persists_to_disk(self, bennu_instance):
        b = bennu_instance
        h = b.propose_hypothesis("Persistent hypothesis")

        # New registry instance reads from same path
        from bennu.core.hypothesis_registry import HypothesisRegistry

        reg2 = HypothesisRegistry(b.hypothesis_registry.path)
        loaded = reg2.get(h.hypothesis_id)
        assert loaded.statement == "Persistent hypothesis"


# ------------------------------------------------------------------ #
# add_evidence
# ------------------------------------------------------------------ #


class TestAddEvidence:
    def test_with_uuid(self, bennu_instance):
        b = bennu_instance
        h = b.propose_hypothesis("Test")
        b.add_evidence(h.hypothesis_id, "query", "result", True, 0.8)
        loaded = b.hypothesis_registry.get(h.hypothesis_id)
        assert len(loaded.evidence) == 1
        assert loaded.evidence[0].supports is True

    def test_with_string_uuid(self, bennu_instance):
        b = bennu_instance
        h = b.propose_hypothesis("Test")
        b.add_evidence(str(h.hypothesis_id), "query", "result", False, 0.6)
        loaded = b.hypothesis_registry.get(h.hypothesis_id)
        assert len(loaded.evidence) == 1
        assert loaded.evidence[0].supports is False

    def test_unknown_hypothesis_raises(self, bennu_instance):
        b = bennu_instance
        with pytest.raises(KeyError):
            b.add_evidence(uuid.uuid4(), "q", "r", True, 0.5)


# ------------------------------------------------------------------ #
# list_hypotheses / hypothesis_summary
# ------------------------------------------------------------------ #


class TestListAndSummary:
    def test_list_empty(self, bennu_instance):
        b = bennu_instance
        assert b.list_hypotheses() == []

    def test_list_returns_all(self, bennu_instance):
        b = bennu_instance
        b.propose_hypothesis("H1")
        b.propose_hypothesis("H2")
        assert len(b.list_hypotheses()) == 2

    def test_summary_empty(self, bennu_instance):
        b = bennu_instance
        assert b.hypothesis_summary() == "No hypotheses registered."

    def test_summary_with_evidence(self, bennu_instance):
        b = bennu_instance
        h = b.propose_hypothesis("NiFe Group 4 is energy-conserving")
        b.add_evidence(h.hypothesis_id, "q", "r", True, 0.9)
        summary = b.hypothesis_summary()
        assert "NiFe Group 4" in summary
        assert "+1/-0" in summary


# ------------------------------------------------------------------ #
# render_provenance
# ------------------------------------------------------------------ #


class TestRenderProvenance:
    def test_returns_mermaid(self, bennu_instance):
        b = bennu_instance
        mermaid = b.render_provenance()
        assert "graph TD" in mermaid

    def test_with_title(self, bennu_instance):
        b = bennu_instance
        mermaid = b.render_provenance(title="My Analysis")
        assert "My Analysis" in mermaid

    def test_includes_hypotheses(self, bennu_instance):
        b = bennu_instance
        b.propose_hypothesis("Test hypothesis for rendering")
        mermaid = b.render_provenance()
        assert "Test hypothesis for rendering" in mermaid

    def test_writes_to_file(self, bennu_instance, tmp_path):
        b = bennu_instance
        out = tmp_path / "provenance.mermaid"
        mermaid = b.render_provenance(output_path=str(out))
        assert out.exists()
        assert out.read_text() == mermaid


# ------------------------------------------------------------------ #
# log_provenance
# ------------------------------------------------------------------ #


class TestLogProvenance:
    def test_basic(self, bennu_instance):
        b = bennu_instance
        entry = b.log_provenance("Count proteins", "42 found")
        assert entry.query == "Count proteins"
        assert entry.results_summary == "42 found"
        assert entry.entry_id is not None

    def test_appears_in_session(self, bennu_instance):
        b = bennu_instance
        b.log_provenance("q1", "r1")
        b.log_provenance("q2", "r2")
        assert len(b.session.get_provenance()) == 2

    def test_parent_chaining_with_uuids(self, bennu_instance):
        b = bennu_instance
        e1 = b.log_provenance("Step 1", "result 1")
        e2 = b.log_provenance("Step 2", "result 2", parent_ids=[e1.entry_id])
        assert e2.parent_ids == [e1.entry_id]

    def test_parent_chaining_with_strings(self, bennu_instance):
        b = bennu_instance
        e1 = b.log_provenance("Step 1", "result 1")
        e2 = b.log_provenance("Step 2", "result 2", parent_ids=[str(e1.entry_id)])
        assert e2.parent_ids == [e1.entry_id]

    def test_no_parent_ids(self, bennu_instance):
        b = bennu_instance
        entry = b.log_provenance("q", "r")
        assert entry.parent_ids == []


# ------------------------------------------------------------------ #
# Manifest integration
# ------------------------------------------------------------------ #


class TestManifestHypothesisIntegration:
    def test_empty_manifest_has_hypotheses_section(self, bennu_instance):
        b = bennu_instance
        assert "hypotheses" in b.manifest.data

    def test_resume_includes_hypotheses(self, bennu_instance):
        b = bennu_instance
        h = b.propose_hypothesis("Dockerins indicate syntrophy")
        b.add_evidence(h.hypothesis_id, "Count dockerins", "15 genomes", True, 0.8)

        status = b.resume()
        assert "Hypotheses" in status
        assert "Dockerins indicate syntrophy" in status
        assert "+1/-0" in status

    def test_resume_omits_hypotheses_when_none(self, bennu_instance):
        b = bennu_instance
        status = b.resume()
        assert "### Hypotheses" not in status

    def test_manifest_counts_updated(self, bennu_instance):
        b = bennu_instance
        b.propose_hypothesis("H1")
        b.propose_hypothesis("H2")
        _ = b.resume()  # triggers count update
        assert b.manifest.data["hypotheses"]["count"] == 2
        assert b.manifest.data["hypotheses"]["active"] == 2
