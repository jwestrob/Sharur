"""Genome ring: contigs, genes by strand and category, linked feature marks."""

from sharur.browser.genome_ring import render


def test_ring_draws_contigs_genes_and_marks():
    data = {
        "contigs": [("c1", 6000), ("c2", 3000)],
        "genes": [("p1", "c1", 1, 900, "+", 300), ("p2", "c1", 1000, 1900, "-1", 300), ("p3", "c2", 1, 9000, "+", 3000)],
        "categories": {"p1": "metabolism"},
        "marks": [{"kind": "giant", "contig": "c2", "start": 1, "end": 9000, "label": "p3 · 3,000 aa", "url": "/protein/p3"},
                  {"kind": "defense", "contig": "nowhere", "start": 1, "end": 2, "label": "x", "url": "/call/x"}],
    }
    svg, legend, marks = render(data, title="Genome", subtitle="2 contigs")
    assert svg.count('class="ring-contig') == 2 and svg.count('class="ring-gene"') == 3
    assert 'href="/protein/p3"' in svg and "/call/x" not in svg       # marks on unknown contigs are skipped
    assert "Metabolism" in legend and "Protein ≥ 3,000 aa" in marks and "Defense" not in marks
