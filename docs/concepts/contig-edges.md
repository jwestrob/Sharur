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

## Circular and linear contigs

A gene on a circular contig has no contig end, so its `edge_status` is `circular`. `contigs.is_circular`
comes from the first of these that marks the contig when Stage 07 or `sharur backfill-contig-context`
reads the assembly:

1. The assembly header: `circular=yes|no` (metaMDBG, Flye style), `suggestCircular=yes|no` (Canu,
   Trycycler), `[topology=circular|linear]` (NCBI style).
2. A contig-name suffix convention the dataset records. Each belongs to one assembler's or project's
   naming, so it applies only to datasets that record it:

   | Convention | Circular | Linear |
   |---|---|---|
   | `upper_c_rc_l` (renamed contigs, uppercase suffixes) | `..._C`, `..._RC` | `..._L` |
   | `metamdbg` | `ctg<N>_<cov>x_c`, `..._rc` | `ctg<N>_<cov>x_l` |
   | `hifiasm` (hifiasm-meta) | `[sp]tg<N>c` | `[sp]tg<N>l` |

3. A curated `circular` token in the name (`Genome_2_circular_FINAL`; `noncircular` reads as linear).

With no convention recorded, a name ending in `_C` or `_L` says nothing about topology. Record the
conventions when building (`sharur-ingest --contig-name-topology upper_c_rc_l`; a `--force` rebuild keeps what
the replaced database recorded) or on an existing database:

```bash
sharur contig-name-topology --db data/my_dataset/sharur.duckdb             # report matches per convention
sharur contig-name-topology --db data/my_dataset/sharur.duckdb --set upper_c_rc_l  # record
sharur seal --db data/my_dataset/sharur.duckdb --force
```

The setting lives in the database's `dataset_meta` table (key `contig_name_topology`). Recording it
leaves stored `is_circular` values as they are; the report lists, per convention, how many contig IDs it
marks and how many stored values differ, and `sharur backfill-contig-context` re-derives them from the
assemblies under the recorded conventions.

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
