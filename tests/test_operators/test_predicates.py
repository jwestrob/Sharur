"""Tests for predicate system."""

import pytest

from bennu.predicates.registry import PredicateRegistry, PredicateDefinition, get_registry, reset_registry
from bennu.predicates.evaluator import evaluate_predicate, compute_predicates_for_protein, compute_all_predicates


class TestPredicateRegistry:
    """Tests for predicate registry."""

    def test_register_and_get(self):
        """Registry should store and retrieve predicates."""
        registry = PredicateRegistry()
        pred = PredicateDefinition(
            predicate_id="test_pred",
            name="Test Predicate",
            description="A test predicate",
            category="custom",
            eval_query="SELECT protein_id FROM proteins WHERE 1=0",
        )
        registry.register(pred)

        retrieved = registry.get("test_pred")
        assert retrieved is not None
        assert retrieved.predicate_id == "test_pred"

    def test_get_nonexistent(self):
        """Registry should return None for unknown predicate."""
        registry = PredicateRegistry()
        assert registry.get("nonexistent") is None

    def test_list_predicates(self):
        """Registry should list all predicates."""
        registry = PredicateRegistry()
        registry.register(PredicateDefinition(
            predicate_id="pred1",
            name="Pred 1",
            description="First",
            category="size",
            eval_query="SELECT 1",
        ))
        registry.register(PredicateDefinition(
            predicate_id="pred2",
            name="Pred 2",
            description="Second",
            category="annotation",
            eval_query="SELECT 1",
        ))

        preds = registry.list_predicates()
        assert len(preds) == 2

    def test_list_by_category(self):
        """Registry should filter by category."""
        registry = PredicateRegistry()
        registry.register(PredicateDefinition(
            predicate_id="size_pred",
            name="Size",
            description="Size predicate",
            category="size",
            eval_query="SELECT 1",
        ))
        registry.register(PredicateDefinition(
            predicate_id="ann_pred",
            name="Annotation",
            description="Annotation predicate",
            category="annotation",
            eval_query="SELECT 1",
        ))

        size_preds = registry.list_predicates(category="size")
        assert len(size_preds) == 1
        assert size_preds[0].predicate_id == "size_pred"


class TestGlobalRegistry:
    """Tests for global predicate registry."""

    def setup_method(self):
        """Reset registry before each test."""
        reset_registry()

    def test_get_registry_initializes_defaults(self):
        """get_registry() should register default predicates."""
        registry = get_registry()
        assert registry.exists("giant")
        assert registry.exists("tiny")
        assert registry.exists("unannotated")
        assert registry.exists("hypothetical")

    def test_get_registry_returns_same_instance(self):
        """get_registry() should return singleton."""
        r1 = get_registry()
        r2 = get_registry()
        assert r1 is r2


class TestPredicateEvaluator:
    """Tests for predicate evaluation."""

    def test_evaluate_giant(self, store):
        """evaluate_predicate() should find giant proteins."""
        result = evaluate_predicate("giant", store)
        # prot_004 (2333aa) and prot_007 (5333aa)
        assert "prot_004" in result
        assert "prot_007" in result
        assert len(result) == 2

    def test_evaluate_tiny(self, store):
        """evaluate_predicate() should find tiny proteins."""
        result = evaluate_predicate("tiny", store)
        # prot_005 (40aa)
        assert "prot_005" in result

    def test_evaluate_unannotated(self, store):
        """evaluate_predicate() should find unannotated proteins."""
        result = evaluate_predicate("unannotated", store)
        # prot_004, 005, 007, 010
        assert len(result) == 4
        assert "prot_004" in result
        assert "prot_005" in result
        assert "prot_007" in result
        assert "prot_010" in result

    def test_evaluate_hypothetical(self, store):
        """evaluate_predicate() should find hypothetical proteins."""
        result = evaluate_predicate("hypothetical", store)
        # prot_003 has "hypothetical" in annotation
        assert "prot_003" in result

    def test_evaluate_with_subset(self, store):
        """evaluate_predicate() should filter to protein subset."""
        result = evaluate_predicate("giant", store, protein_ids=["prot_004", "prot_005"])
        assert "prot_004" in result
        assert "prot_005" not in result

    def test_compute_for_protein(self, store):
        """compute_predicates_for_protein() should return list of predicates."""
        # prot_004 is giant and unannotated
        preds = compute_predicates_for_protein("prot_004", store)
        assert "giant" in preds
        assert "unannotated" in preds

    def test_compute_all_predicates(self, store):
        """compute_all_predicates() should return dict for all proteins."""
        result = compute_all_predicates(store)
        # Should have entries for proteins with any predicate
        assert "prot_004" in result  # giant, unannotated
        assert "giant" in result["prot_004"]
