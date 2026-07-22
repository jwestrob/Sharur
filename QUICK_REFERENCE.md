# Sharur Quick Reference Card

A query/operator cheatsheet for agents working against a Sharur DuckDB. For
narrative guidance see `CLAUDE.md` (routing table) and `docs/`.

## Prefer operators over raw SQL

CLAUDE.md's core rule: reach for the `Sharur` facade before hand-writing SQL.
The facade applies the `COUNT(DISTINCT protein_id)` discipline, predicate
mappings, and neighborhood logic for you.

```python
from sharur.operators import Sharur

b = Sharur("data/my_dataset/sharur.duckdb", read_only=True)

# Predicate search -> records is the public programmatic payload
result = b.search_by_predicates(has=["unannotated", "giant"]); proteins = result.records
b.search_by_predicates(has=["nife_group3", "bidirectional_hydrogenase"])

# Genomic neighborhood (all_annotations pulls every source, not just PFAM)
b.get_neighborhood(protein_id, window=10)
b.get_neighborhood(protein_id, window=5, all_annotations=True)

# Embedding similarity / structure / export
b.find_similar(protein_id, k=20)
b.predict_structure(protein_id)
b.search_foldseek_for_protein(protein_id)
b.export_fasta(protein_ids, "out.faa")
```

Drop to raw SQL (below) only for ad-hoc shapes the facade doesn't cover.
For shell pipelines, the main read commands accept
`--format markdown|json|jsonl|tsv`.

## Operator preflight

```bash
# Typed, non-mutating capability brief
sharur preflight --db data/my_dataset/sharur.duckdb
sharur preflight --db data/my_dataset/sharur.duckdb --format json --skip-tools

# Prepare persistent similarity sidecars for a legacy embedding H5
sharur build-vector-index --db data/my_dataset/sharur.duckdb

# Seal a stable dataset and verify it later
sharur seal --db data/my_dataset/sharur.duckdb
sharur verify-seal data/my_dataset/dataset.seal.json
# Archival mode: stream full SHA-256 over large canonical artifacts
sharur seal --db data/my_dataset/sharur.duckdb --full --force
```

Capability states are `available`, `unavailable`, `stale`, and `failed`. The brief discovers
live `annotations`-table sources and structured caller resources separately; it does not
assume a closed annotation list or treat every annotation source as a system caller.
Seal verification returns non-zero on canonical drift; use `--format json` in automation.

## Schema (live as of schema v4)

### proteins
`protein_id, contig_id, bin_id, gene_index, start, end_coord, strand, sequence, sequence_length, gc_content`
- Coordinates: `start` / `end_coord` (NOT `end`). `gene_index` is the ordinal on the contig.

### contigs
`contig_id, bin_id, length, gc_content, is_circular, taxonomy`

### bins
`bin_id, completeness, contamination, taxonomy, n_contigs, total_length`
- Taxonomy lives on **both** `contigs` and `bins` — join to whichever you already have.

### annotations
`annotation_id, protein_id, source, accession, name, description, evalue, score, start_aa, end_aa`
- `accession` = e.g. `PF00142`; `name` = short label; `score` (NOT `bitscore`).
- Annotation sources are dataset-dependent. Examples include `pfam`, `kofam`,
  `hyddb`, and validated-caller products. Inspect the live set with
  `SELECT DISTINCT source FROM annotations ORDER BY source`; also inspect the
  live curated system tables rather than assuming annotations are the complete caller surface.
- **Always `COUNT(DISTINCT protein_id)`** — repeat domains inflate `COUNT(*)`.

### Predicates (V2 backend + V1 compat)
- `protein_predicates` — compat table: `protein_id, predicates (VARCHAR[]), updated_at`. One row per protein; `predicates` is an array like `['giant','unannotated','kegg:K01591', ...]`.
- `semantic_atoms` — `protein_id, atom_id, facet, relation, source_accession, source_db, evidence_evalue, evidence_score`.
- `semantic_terms` — `protein_id, term_id, term_kind, facet, relation, source_db, source_accession`.
- `semantic_state` — one row per protein (resolved state).
- `system_proteins` — `system_id, protein_id, system_source, position, profile_name, score`.

### Loci & systems
- `loci` — `locus_id, locus_type, contig_id, start, end_coord, confidence, metadata (JSON)`. Live `locus_type` values: `crispr`, `prophage`, `island` (dataset-dependent). **BGCs are NOT here.**
- `bgc_loci` — same shape as `loci`; this is where biosynthetic gene clusters (GECCO) go.
- `locus_proteins` — `locus_id, protein_id, position` (membership join table).
- `defense_systems` — MacSyFinder-validated: `system_id, genome_id, system_type, system_subtype, activity, genes_count, protein_ids, profile_names, sys_beg, sys_end, created_at`. **Only this table holds genuine defense system calls** — never report raw `defensefinder` HMM hits as systems (80%+ FP).
- `secretion_systems` — TXSScan-validated, same columns minus `activity`.

### ELSA synteny sidecar

`DATASET/synteny.duckdb` is the normalized agent query surface. Its
`elsa_runs`, `elsa_blocks`, `elsa_anchor_pairs`, `elsa_clusters`,
`elsa_cluster_loci`, and `elsa_cluster_members` tables are keyed by `run_id`.
Singleton clusters and anchor/context membership are explicit.

```python
b = Sharur("data/DATASET/sharur.duckdb", read_only=True)
b.capabilities().get("elsa_synteny")
b.synteny_for_protein("PROTEIN_ID")
b.synteny_anchor_blocks("PROTEIN_ID")
b.get_synteny_cluster("cluster:SOURCE_ID", run_id="elsa-RUN_ID")
```

ELSA store position indices and Sharur `proteins.gene_index` values are
separate coordinate contracts. Exact relational membership replaces legacy
CSV string scans. Dataset-seal mismatch fails closed; use
`allow_stale_synteny=True` only for a reviewed, explicitly historical run.

## SQL patterns

### Domain search with taxonomy
```sql
SELECT DISTINCT p.protein_id, p.bin_id
FROM proteins p
JOIN annotations a ON p.protein_id = a.protein_id
JOIN bins b ON p.bin_id = b.bin_id
WHERE a.source = 'pfam' AND a.accession = :domain
  AND b.taxonomy LIKE :tax_pattern;     -- e.g. '%Nanobdellota%'
```

### Predicate membership (compat array)
```sql
-- proteins tagged both giant and unannotated
SELECT protein_id FROM protein_predicates
WHERE list_contains(predicates, 'giant')
  AND list_contains(predicates, 'unannotated');
```

### Strand-aware neighbors
```python
def upstream_window(protein, n_bp):
    if protein.strand == "+":
        return (protein.start - n_bp, protein.start)
    else:  # "-" strand: upstream is HIGHER coordinates
        return (protein.end_coord, protein.end_coord + n_bp)
```

### Spatial window (linear contig)
```sql
SELECT p2.*
FROM proteins p1
JOIN proteins p2 ON p1.contig_id = p2.contig_id
WHERE p1.protein_id = :anchor
  AND p2.start    >= p1.start    - :window
  AND p2.end_coord <= p1.end_coord + :window
ORDER BY p2.start;
```
(MAG contigs are effectively all linear — `is_circular` is unset across these
datasets. Skip wrap-around logic unless a contig actually flags circular.)

### Operon expansion (adjacent same-strand genes within a gap)
```sql
WITH RECURSIVE operon AS (
    SELECT protein_id, contig_id, start, end_coord, strand, 0 AS depth
    FROM proteins WHERE protein_id = :seed
    UNION ALL
    SELECT p.protein_id, p.contig_id, p.start, p.end_coord, p.strand, o.depth + 1
    FROM proteins p
    JOIN operon o ON p.contig_id = o.contig_id AND p.strand = o.strand
    WHERE (p.start - o.end_coord BETWEEN 0 AND :gap_threshold
           OR o.start - p.end_coord BETWEEN 0 AND :gap_threshold)
      AND p.protein_id != o.protein_id
      AND o.depth < 50  -- safety limit
)
SELECT DISTINCT * FROM operon;
```

### Validated defense systems (never raw HMM hits)
```sql
SELECT genome_id, system_type, system_subtype, genes_count
FROM defense_systems
ORDER BY genome_id;
```

## Common errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| "No results" for an obvious query | Taxonomy string mismatch | Check GTDB format; use `LIKE '%Name%'` |
| Inflated counts | `COUNT(*)` over repeat domains | `COUNT(DISTINCT protein_id)` |
| Wrong neighbors | Ignored strand | Strand-aware upstream/downstream |
| `column "end" does not exist` | Wrong column name | Use `end_coord` |
| BGC query returns nothing | Queried `loci` | BGCs live in `bgc_loci` |
| Defense FPs in report | Raw `defensefinder` hits | Use `defense_systems` table only |
