# Contig edges

Metagenome-assembled genomes arrive in pieces. A gene near the end of a contig may have lost its neighbors, or part of itself, to an assembly break. Sharur records where each gene sits so that this is visible when you read a protein, a neighborhood or a pathway.

## What is recorded

| Field | Meaning |
|---|---|
| `genes_to_start`, `genes_to_end` | Genes between this one and each contig end, ordered by coordinate |
| `bp_to_start`, `bp_to_end` | Base pairs to each end. `bp_to_end` needs the assembly's contig length |
| `truncated_start`, `truncated_end` | Prodigal's `partial=` flag: the gene runs off that end of the contig (or into an assembly gap) |
| `edge_status` | `truncated`, `contig_edge` (within 3 genes and 5 kb of an end), `interior`, `circular`, or `no_coordinates` for protein-only records |

`sharur card` prints the contig position, `sharur neighborhood` notes when a contig end falls inside the window and which genes run off it, and `sharur.contig_context.edge_context()` returns the fields for any set of proteins.

## Module completeness

When two or more found genes of an incomplete KEGG module sit together (within 5 genes) and that cluster reaches a contig end, `sharur modules` lists them under `found_at_contig_edge` and notes that the missing genes could lie beyond the contig end.

Read this as context for a missing step. In fragmented MAGs, clusters at contig ends are common among complete modules too: in one test on 120 MAGs, 49% of complete modules had such a cluster, against 57% of modules missing one step. A cluster at an edge says a break *could* explain a gap; a missing step on a long, intact contig is the stronger observation.

## Contig lengths in older databases

Databases built before schema 8 store each contig's length as the end of its last gene, so the distance from the last gene to the contig end is unknown, and Prodigal's truncation flags are absent. New ingests read both. For an existing dataset:

```bash
sharur backfill-contig-context --db data/my_dataset/sharur.duckdb \
  --assemblies path/to/assemblies --proteins data/my_dataset/stage03_prodigal
sharur seal --db data/my_dataset/sharur.duckdb --force
```

Assemblies are matched to bins by file name (`BIN_ID.fna`, `.fa` or `.fasta`, optionally gzipped). The command writes the canonical database, so run it in a maintenance window. It refuses an assembly whose contigs are shorter than the genes on them.
