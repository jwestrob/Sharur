"""Model-parsing invariants for the co-location engine.

The engine reimplements MacSyFinder, so every attribute it fails to read is a
silent divergence from the reference: a missed per-gene spacing override or a
dropped exchangeable does not raise, it just produces a different system. These
tests parse the installed model corpus and assert the attributes survive.

Skipped when the MacSyFinder model corpus is not installed.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from sharur.colocation import _parse_model_xml

MODELS_ROOT = Path.home() / ".macsyfinder" / "models"
TXSSCAN = MODELS_ROOT / "TXSScan" / "definitions"
DEFENSE = MODELS_ROOT / "defense-finder-models" / "definitions"

pytestmark = pytest.mark.skipif(
    not TXSSCAN.is_dir(), reason="MacSyFinder model corpus not installed"
)


def _models(root: Path, family: str):
    for xml in sorted(root.rglob("*.xml")):
        model = _parse_model_xml(xml, family)
        if model is not None:
            yield xml, model


class TestModelParsing:
    def test_every_txsscan_model_parses(self):
        parsed = list(_models(TXSSCAN, "TXSScan"))
        assert parsed, "no TXSScan models parsed"
        for xml, model in parsed:
            assert model.genes, f"{xml.name} parsed with no genes"
            assert model.inter_gene_max_space >= 1
            assert model.min_genes_required >= 1

    def test_quorums_are_not_silently_defaulted(self):
        """min_mandatory must never exceed min_genes; a swap inverts the filter."""
        for xml, model in _models(TXSSCAN, "TXSScan"):
            assert model.min_mandatory_genes_required <= model.min_genes_required, xml.name

    def test_per_gene_spacing_override_is_read(self):
        """A gene may widen its own spacing beyond the model default.

        Missing this silently splits clusters the reference would keep whole.
        """
        conj = TXSSCAN / "bacteria" / "diderm" / "CONJ.xml"
        if not conj.exists():
            pytest.skip("CONJ model not present")
        model = _parse_model_xml(conj, "TXSScan")
        assert model.inter_gene_max_space == 30
        overrides = {
            g.name: g.inter_gene_max_space
            for g in model.genes
            if getattr(g, "inter_gene_max_space", None) is not None
        }
        assert overrides, "per-gene inter_gene_max_space override was dropped"
        assert max(overrides.values()) > model.inter_gene_max_space

    def test_multi_loci_is_read_where_declared(self):
        """multi_loci means several loci on ONE replicon -- never across contigs."""
        flagged = [m.name for _x, m in _models(TXSSCAN, "TXSScan") if m.multi_loci]
        assert flagged, "no TXSScan model parsed with multi_loci=True"

    def test_defensefinder_corpus_is_single_locus(self):
        """DefenseFinder declares no multi_loci and no loner anywhere.

        This is why it is structurally immune to cross-contig assembly, and why
        it is the correct baseline to validate the engine against first.
        """
        if not DEFENSE.is_dir():
            pytest.skip("defense-finder models not installed")
        for xml, model in _models(DEFENSE, "defense-finder-models"):
            assert not model.multi_loci, f"{xml.name} unexpectedly multi_loci"

    def test_exchangeables_do_not_collapse_into_the_parent(self):
        """An exchangeable substitutes for its parent and must remain resolvable."""
        conj = TXSSCAN / "bacteria" / "diderm" / "CONJ.xml"
        if not conj.exists():
            pytest.skip("CONJ model not present")
        model = _parse_model_xml(conj, "TXSScan")
        names = {g.name for g in model.genes}
        exch: set[str] = set()
        for g in model.genes:
            exch |= set(getattr(g, "exchangeables", ()) or ())
        assert "T4SS_t4cp1" in names
        assert exch, "exchangeables were dropped entirely"
