# First look at a dataset

These commands assume an ingested dataset at `data/my_dataset/sharur.duckdb`. All of them read the database without changing it.

## What is in this dataset?

```bash
sharur describe --db data/my_dataset/sharur.duckdb
```

`describe` lists the annotation sources and how many proteins each covers, which curated callers ran (defense systems, secretion systems, hydrogenase classification and so on), the genome metadata available, and whether the stored predicates match the installed predicate maps. Start here: it tells you which questions the dataset can answer.

## One protein at a glance

```bash
sharur card PROTEIN_ID --db data/my_dataset/sharur.duckdb
```

A card shows a protein's genome and contig, its length and position (including how close it sits to a contig end; see [Contig edges](../concepts/contig-edges.md)), its annotation hits, and its predicates with the evidence behind each one, in one screen.

## Why does a protein carry a label?

```bash
sharur why PROTEIN_ID sam_binding --db data/my_dataset/sharur.duckdb
```

`why` prints each path from an annotation hit to the predicate: the hit, the map entry and its evidence (for example a GO term, a KEGG EC number or a Swiss-Prot consensus count), and any is-a expansion. [How predicates are built](../predicate_construction.md) explains each evidence type.

## Pathways

```bash
sharur setup-kegg     # once per machine
sharur modules --db data/my_dataset/sharur.duckdb --bin GENOME_ID --min-completeness 0.75
```

Module completeness follows KEGG's module definitions step by step. Genomes assembled from metagenomes are often incomplete, so a missing step means "not detected in this genome", which is weaker than "absent". When found genes of an incomplete module cluster at a contig end, the output says so.

## Searching

```bash
sharur search --help
sharur neighborhood --help
sharur architecture "Big_* {5,} . VWA" --db data/my_dataset/sharur.duckdb
```

The same operations are available from Python:

```python
from sharur.operators import Sharur

b = Sharur("data/my_dataset/sharur.duckdb", read_only=True)
hits = b.search_by_predicates(has=["nife_group3"])
b.get_neighborhood(hits.records[0]["protein_id"], window=10)
b.find_similar(hits.records[0]["protein_id"], k=20)
b.search_architecture("TPR_* {10,}")          # ordered-domain patterns
```

[Domain-architecture search](../guides/architecture-search.md) covers the pattern language.

## Before you report numbers

Read [Interpreting annotations](../biological_interpretation.md). It lists the known annotation traps, keyed on what you observe, and the rules for naming systems and pathways.
