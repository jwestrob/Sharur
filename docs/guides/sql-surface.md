# SQL surface

Sharur's neighborhood, slot and presence logic is available as SQL macros, so an analysis
written in SQL gets the same strand handling, feature matching and evidence labels as the
Python operators. The macros live on the connection only (DuckDB `TEMP` macros): they install
on a read-only connection, leave the database file and its seal untouched, and appear in no
table listing.

## Opening a connection

**Python.** Every `DuckDBStore` connection carries the surface, so `Sharur(..., read_only=True)`
has it:

```python
from sharur.operators import Sharur

b = Sharur("data/DATASET/sharur.duckdb", read_only=True, duckdb_threads=4, duckdb_memory_limit="8GB")
b.store.execute_df("FROM slot_tally('pfam:PF00069', 'pfam:PF00072', split := 'K00001')")
```

`sharur.sql.connect(path)` returns a plain read-only DuckDB connection with the surface
installed. `sql_prelude=False` on `DuckDBStore`, or `SHARUR_SQL_PRELUDE=0`, leaves it out.

**Command line.**

```bash
sharur sql query "FROM neighbors('PF00069', k := 2) LIMIT 20" --db data/DATASET/sharur.duckdb
sharur sql objects --db data/DATASET/sharur.duckdb      # what is installed, with signatures
```

`sharur sql query` opens the database read-only and prints markdown, TSV or JSON.

**Any DuckDB client** (the DuckDB CLI, an MCP DuckDB server, a notebook): attach the database,
make it the default catalog, then run the prelude file.

```sql
ATTACH 'data/DATASET/sharur.duckdb' AS ds (READ_ONLY);
USE ds;
.read sharur/sql/prelude.sql
```

A client without `.read` takes the file's text as one multi-statement query;
`sharur sql prelude --out prelude.sql` writes the same script with a registry macro appended.
`FROM sharur_sql()` then lists every object with its signature. DuckDB 1.3 and later bind a table
macro's tables when the macro is created, so on a database without the protein-family or
remote-homology layer those macros cannot be created; `sharur sql prelude --db
data/DATASET/sharur.duckdb --out prelude.sql` writes the script for that database without them,
and Python connections leave them out on their own (`sharur sql objects` lists them as not
installed).

## Feature specs

Every macro that takes features accepts a comma-separated string or a list:

| Spec | Matches |
|---|---|
| `PF00069`, `K00001`, `Pkinase` | accession or name in any annotation source |
| `pfam:Pkinase`, `kofam:K00001` | that source only, by accession or name |
| `pfam:*` | any hit of that source |
| `ko:K00001` | KofamScan or KEGG |
| `predicate:hydrogenase` | `protein_predicates` |
| `defense:CBASS`, `secretion:TYPE` | a curated `defense_systems` / `secretion_systems` call, by type or subtype |
| `crispr:I-B` | a `crispr_cas_systems` call with status `crispr_cas` or `cas` |
| `cdd:cd00009`, `cdd:*` | CDD hits (remote homology); an unqualified spec leaves them out |
| a protein ID | that protein (anchor arguments) |

Each match carries `evidence`: `observed` for per-domain hits, `caller` for curated system
calls and their annotation projections, `predicate` for predicates, `resembles` for
remote-homology rows ([Remote homology](remote-homology.md)). A named system rests on
`caller` rows; `observed` rows are domains consistent with a name; `resembles` rows record a
similarity to a profile and name nothing.

## Objects

| Object | Returns |
|---|---|
| `genes()` | the `proteins` table without its sequence column, plus `strand_sign` |
| `genomes()` | `bins` with `domain`, `phylum`, `class`, `"order"`, `family`, `genus`, `species` |
| `feature_proteins(features)` | one row per matching protein: `feature`, `protein_id`, `bin_id`, `contig_id`, `gene_index`, `strand`, `source`, `accession`, `name`, `evidence` |
| `feature_presence(features, min_completeness := 0, max_contamination := 100)` | every genome × feature: `present`, `n_proteins`, `sources`, `evidence`, completeness, contamination, GTDB ranks |
| `neighbors(anchor, k := 5)` | genes within `k` positions of each anchor, `offset` in the anchor's orientation (negative = upstream), `same_strand`, distance to each contig end, and each gene's `pfam`, `ko`, `ko_name`, `systems`, `other`, `resembles` (CDD hits) |
| `neighbor_features(anchor, k := 5, src := 'pfam')` | per oriented offset and strand: each feature's number of genomes and anchors |
| `slot(anchor, side := 'upstream', partner := NULL, max_span := 10, min_completeness := 0, max_contamination := 100)` | per anchor protein: the gene beside it (`side` = `upstream` or `downstream`), or every gene between it and its nearest `partner` within `max_span` genes (`side` = `upstream`, `downstream` or `either`); occupant strand, length, annotations, `at_contig_edge`, genome completeness and ranks |
| `slot_tally(anchor, occupant, side, partner, rank := 'order', split := NULL, max_span, min_completeness, max_contamination)` | per taxon (and per `split` feature presence): genomes with an anchor, with an `occupant` match on the same strand or either strand, with an occupant lacking any Pfam hit, at a contig edge, without a partner in range |
| `protein_summary(ids)`, `anchor_proteins(anchor)` | the building blocks above |
| `family_proteins(families)`, `family_presence(families, min_completeness := 0, max_contamination := 100)`, `families_of(ids)` | protein families in the `feature_proteins` / `feature_presence` shapes, and each protein's family in every clustering run; need the family layer ([Protein families](protein-families.md)) |
| `resembles_of(ids, min_probability := 90, min_cols := 30)`, `resembles_proteins(targets, ...)`, `resembles_presence(targets, ..., min_completeness := 0, max_contamination := 100)` | what proteins resemble (their families' Pfam and PDB hits at the rule, their own CDD hits), proteins resembling targets in the `feature_proteins` shape, and the genome × target grid; `evidence` = `resembles`; need the remote-homology layer ([Remote homology](remote-homology.md)) |
| `family_annotations_via(target_run)` | family-tier hits of other clustering runs on the families of `target_run` that share the searched family's representative |
| `strand_sign(s)`, `id_list(x)`, `rank_of(taxonomy, rank)` | scalar helpers |
| `write_fasta(path, ids)` | Python connections only: writes FASTA, returns `{path, n_sequences, n_requested, n_missing}` |

Gene order is the order of called genes on a contig (`gene_index`, then start). An anchor with
an unrecognized strand is read in `+` orientation. Strand codes `+`/`-` and `1`/`-1` both work.

## Examples

The examples run on the four-genome fixture of `tests/test_sql_surface.py`, whose feature
names are invented: an anchor family (KO `K90001`, Pfam `AnchorFam`), a neighbor family
(Pfam `NeighborFam`, `NeighborFam2`), a partner (`K90002`) and a genome-level split feature
(`K90003`). Two genomes sit in order `OrdA`, two in `OrdB`:

```
g1 (OrdA):  p0+  p1+ NeighborFam  p2+ anchor  p3-  p4+ partner  p5+
g2 (OrdA):  q0- partner  q1- anchor  q2- NeighborFam2  q3+  q4-
g3 (OrdB):  r0+ anchor  r1+ split
g4 (OrdB):  s0- (no annotation)  s1+ anchor  s2+ split ;  second contig: s3+ (CBASS system)  s4+
```

### What sits upstream of the anchor, by order

```sql
SELECT taxon, genome_has_split, genomes_with_anchor, genomes_occupant_same_strand,
       genomes_occupant_no_pfam, genomes_at_contig_edge
FROM slot_tally('K90001,pfam:AnchorFam', 'pfam:NeighborFam,pfam:NeighborFam2',
                side := 'upstream', rank := 'order', split := 'K90003');
```

| taxon | genome_has_split | genomes_with_anchor | genomes_occupant_same_strand | genomes_occupant_no_pfam | genomes_at_contig_edge |
|---|---|---|---|---|---|
| OrdA | False | 2 | 2 | 0 | 0 |
| OrdB | True | 2 | 0 | 1 | 1 |

Upstream is read in each anchor's orientation: for the minus-strand anchor `q1` it is `q2`.
In `OrdB`, `r0` starts its contig and `s1`'s upstream gene lies on the other strand with no
Pfam hit. Add `min_completeness := 90, max_contamination := 5` to restrict to high-quality
genomes.

`slot()` returns the rows behind the tally, so further columns are a `GROUP BY` away:

```sql
WITH s AS (FROM slot('K90001,pfam:AnchorFam')),
     split AS (SELECT DISTINCT bin_id FROM feature_proteins('K90003'))
SELECT "order", bin_id IN (SELECT bin_id FROM split) AS has_split,
       count(DISTINCT bin_id) AS genomes_with_anchor,
       count(DISTINCT bin_id) FILTER (WHERE occupant_same_strand AND list_contains(pfam, 'NeighborFam')) AS neighborfam_upstream,
       count(DISTINCT bin_id) FILTER (WHERE occupant_id IS NOT NULL AND len(pfam) = 0) AS upstream_has_no_pfam,
       count(DISTINCT bin_id) FILTER (WHERE at_contig_edge) AS at_contig_edge
FROM s GROUP BY ALL ORDER BY ALL;
```

| order | has_split | genomes_with_anchor | neighborfam_upstream | upstream_has_no_pfam | at_contig_edge |
|---|---|---|---|---|---|
| OrdA | False | 2 | 1 | 0 | 0 |
| OrdB | True | 2 | 0 | 1 | 1 |

### What sits between the anchor and its partner

```sql
SELECT anchor_id, partner_id, partner_offset, n_between, slot_position, occupant_id, occupant_same_strand
FROM slot('K90001', partner := 'K90002', side := 'either');
```

| anchor_id | partner_id | partner_offset | n_between | slot_position | occupant_id | occupant_same_strand |
|---|---|---|---|---|---|---|
| p2 | p4 | 2 | 1 | 1 | p3 | False |
| q1 | q0 | 1 | 0 |  |  |  |
| s1 |  |  | 0 |  |  |  |

One row per gene in the slot; an adjacent partner gives an empty slot, and an anchor without a
partner within `max_span` genes gives a row without one.

### Genome-level presence in one shape

```sql
SELECT bin_id, feature, present, n_proteins, evidence, "order"
FROM feature_presence('K90003,defense:CBASS') ORDER BY bin_id, feature;
```

| bin_id | feature | present | n_proteins | evidence | order |
|---|---|---|---|---|---|
| g1 | K90003 | False | 0 |  | OrdA |
| g1 | defense:CBASS | False | 0 |  | OrdA |
| g2 | K90003 | False | 0 |  | OrdA |
| g2 | defense:CBASS | False | 0 |  | OrdA |
| g3 | K90003 | True | 1 | ["observed"] | OrdB |
| g3 | defense:CBASS | False | 0 |  | OrdB |
| g4 | K90003 | True | 1 | ["observed"] | OrdB |
| g4 | defense:CBASS | True | 1 | ["caller"] | OrdB |

Every genome appears for every feature, so presence/absence tables, alternation counts and
stratified tests start from one query; completeness and ranks ride along for stratification.

### Context of one protein

```sql
SELECT "offset", protein_id, same_strand, pfam, ko, ko_name
FROM neighbors('p2', k := 2);
```

| offset | protein_id | same_strand | pfam | ko | ko_name |
|---|---|---|---|---|---|
| -2 | p0 | True | [] | [] |  |
| -1 | p1 | True | ["NeighborFam", "NeighborFam_C"] | ["K90004"] |  |
| 0 | p2 | True | [] | ["K90001"] | anchor family enzyme |
| 1 | p3 | False | [] | [] |  |
| 2 | p4 | True | [] | ["K90002"] |  |

Pfam names follow domain order. `genes_upstream_to_edge` and `genes_downstream_to_edge` give
the anchor's distance to each contig end; a window that reaches an end returns fewer genes
(see [Contig edges](../concepts/contig-edges.md)).

## Sequences stay out of results

No macro returns a sequence: `genes()` is the protein table to query, and the docs and
examples use it. FASTA leaves the database as a file. In a Python connection:

```sql
SELECT write_fasta('/abs/path/upstream_neighbors.faa',
                   (SELECT list(DISTINCT occupant_id) FROM slot('K90001,pfam:AnchorFam')
                    WHERE occupant_same_strand
                      AND occupant_id IN (SELECT protein_id FROM feature_proteins('pfam:NeighborFam,pfam:NeighborFam2'))));
-- {'path': '/abs/path/upstream_neighbors.faa', 'n_sequences': 2, 'n_requested': 2, 'n_missing': 0}
```

Pass the IDs as one list; each call writes the whole file. `sharur.sql.write_fasta(store, ids,
path)` does the same from Python. In a generic client, `COPY` writes the file and reports only
a row count:

```sql
COPY (SELECT line FROM (
        SELECT protein_id, unnest(['>' || protein_id, rtrim(sequence, '*')]) AS line,
               unnest([1, 2]) AS k
        FROM proteins WHERE protein_id IN (SELECT protein_id FROM feature_proteins('K90001')))
      ORDER BY protein_id, k)
TO '/abs/path/anchors.faa' (FORMAT csv, HEADER false, QUOTE '', DELIMITER '\t');
```

`sharur sql query` refuses to print a result with a sequence column or a run of 30 or more
residue letters.

## In verification queries

`sharur reverify` and the review verifier install the surface on their read-only replay
connections, so a finding's verification query can call the macros:

```json
{"claim": "Both OrdA genomes carry a neighbor-family gene upstream of the anchor on the same strand",
 "query": "SELECT genomes_occupant_same_strand FROM slot_tally('K90001,pfam:AnchorFam', 'pfam:NeighborFam,pfam:NeighborFam2') WHERE taxon = 'OrdA'",
 "expected": 2}
```

## Speed

On a dataset of about 3.6 million proteins and 6 million annotation rows (4 threads,
read-only): `feature_proteins` and `feature_presence` 0.1–0.3 s, `neighbors` of one protein
0.2 s, `neighbors` of about 1,300 anchors 0.4 s, `slot` and `slot_tally` over all genomes
0.7–1.0 s, a partner-mode `slot` 0.8 s.

## Adding a macro

`sharur/sql/prelude.sql` holds one statement per `-- @object NAME KIND: DOC` marker. Use
`CREATE OR REPLACE TEMP MACRO` only (views and tables would need write access and would enter
`information_schema` and dataset seals), keep sequences out of the result, and add a test to
`tests/test_sql_surface.py`. Mark a CTE `MATERIALIZED` when the macro reads it more than once.
