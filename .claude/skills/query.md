# Quick Query Skill

Fast, ad-hoc queries against the Sharur database. No ceremony, just answers.

**CRITICAL: You are a leaf agent. DO NOT spawn sub-agents or use the Task tool.**

**CONCURRENCY: DuckDB does not support concurrent writes. Only ONE agent should access a database at a time. The coordinator must run DB-accessing skills sequentially, not in parallel.**

> **Mandatory:** Follow the shared validation protocols in `_validation_protocols.md`.
> Verify accession names before reporting. Use COUNT(DISTINCT protein_id) for protein
> counts. Apply Context-First protocol for annotations averaging >10 hits/genome.

---

## Usage

```
/q giant proteins on scaffold_9
/q proteins with cas_domain but annotated as transposase
/q top 20 PFAM domains by frequency
/q unannotated proteins larger than 1000aa
/q how many genomes have hydrogenase?
/q show me the neighborhood around protein_123
```

---

## Prompt

You are running quick queries against a Sharur metagenomic database. Answer the user's question directly and concisely.

**Database**: Use the Sharur instance at `data/{DATASET}/sharur.duckdb` (infer from context or ask).

### Query Patterns

**Predicate searches:**
```python
from sharur.operators import Sharur
b = Sharur("data/DATASET/sharur.duckdb")

# Combine predicates with AND logic
b.search_by_predicates(has=["giant", "unannotated"])
b.search_by_predicates(has=["cas_domain"], lacks=["crispr_associated"])
b.search_by_predicates(has=["hydrogenase", "multi_domain"])
```

**SQL queries for statistics:**
```python
# Domain frequency
b.store.execute("""
    SELECT name, COUNT(*) as n
    FROM annotations
    WHERE source = 'pfam'
    GROUP BY name
    ORDER BY n DESC
    LIMIT 20
""")

# Predicate counts
b.store.execute("""
    SELECT pred, COUNT(*) as n
    FROM protein_predicates, LATERAL unnest(predicates) AS t(pred)
    GROUP BY pred
    ORDER BY n DESC
    LIMIT 30
""")

# Per-genome counts
b.store.execute("""
    SELECT bin_id, COUNT(*) as n_proteins,
           SUM(CASE WHEN 'unannotated' = ANY(pp.predicates) THEN 1 ELSE 0 END) as n_unannotated
    FROM proteins p
    JOIN protein_predicates pp ON p.protein_id = pp.protein_id
    GROUP BY bin_id
""")
```

**Neighborhood lookups:**
```python
# Quick neighborhood view
print(b.get_neighborhood("protein_id", window=10))

# With visualization
b.visualize_neighborhood("protein_id", window=12)
```

**Protein details:**
```python
# Full protein info
print(b.get_protein("protein_id", verbosity=2))

# Just annotations
b.store.execute("""
    SELECT source, name, accession, evalue
    FROM annotations
    WHERE protein_id = 'protein_id'
    ORDER BY evalue
""")
```

### Response Format

- **Be concise** - answer the question directly
- **Show counts** - include numbers when relevant
- **List examples** - show 3-5 representative hits
- **Suggest follow-ups** - if the query reveals something interesting

### Example Responses

**Query:** "how many genomes have group 4 hydrogenase?"

**Response:**
```
23 of 29 genomes (79%) have group 4 hydrogenase (mbh_hydrogenase predicate).

Top genomes by count:
- GCA_003598175.1: 12 proteins
- GCA_020430105.1: 8 proteins
- GCA_016934335.1: 7 proteins

This is expected for syntrophs - group 4 hydrogenases are energy-conserving
and H2-evolving, critical for syntrophic metabolism.
```

**Query:** "giant unannotated proteins"

**Response:**
```
Found 47 giant (>1000 aa) unannotated proteins across 18 genomes.

Largest:
- prot_12345 (3,842 aa) on GCA_003598175.1 - has DUF2341 domain
- prot_67890 (2,156 aa) on GCA_020430105.1 - no domain hits
- prot_11111 (1,892 aa) on GCA_016934335.1 - has TPR repeats

Consider running /characterize on these for structure-based annotation.
```
