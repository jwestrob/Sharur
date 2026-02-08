# Bennu

> A data plane for agent-driven metagenomic discovery

Bennu makes large metagenomic datasets navigable by AI agents. It combines a DuckDB relational store, LanceDB vector store, and a functional predicate system into an operator framework that agents (Claude Code, Codex, etc.) use to search, characterize, and compare proteins across hundreds of genomes.

## What it does

Given a set of metagenome-assembled genomes (MAGs), Bennu:

1. **Ingests** proteins, annotations (PFAM, KEGG, HydDB, VOGdb, CAZy, DefenseFinder), CRISPR arrays, biosynthetic gene clusters, and ESM2 embeddings into a unified database
2. **Computes predicates** -- functional tags derived from annotation combinations (e.g., `nife_group3`, `crispr_associated`, `giant_unannotated`) that make semantic search possible
3. **Exposes operators** that agents call to explore the data: search by predicate, navigate genomic neighborhoods, find similar proteins by embedding, detect loci, export results

Agents bring the reasoning; Bennu brings the data access.

## Quick start

```bash
git clone https://github.com/jwestrob/Sharur.git
cd Sharur
pip install -e ".[dev]"
```

### Ingest a dataset

```bash
python scripts/ingest.py \
  --input-dir /path/to/genomes \
  --data-dir data/my_dataset \
  --output data/my_dataset/bennu.duckdb \
  --force
```

This runs Prodigal, Astra (PFAM/KEGG/HydDB), GECCO, dbCAN, minced, ESM2 embeddings, and builds the DuckDB knowledge base. See [`QUICKSTART.md`](QUICKSTART.md) for manual step-by-step instructions.

### Use the operators

```python
from bennu.operators import Bennu

b = Bennu("data/my_dataset/bennu.duckdb")

# Predicate search
hydrogenases = b.search("nife_group3 AND bidirectional_hydrogenase")
giants = b.search("giant AND unannotated")
defense = b.search("crispr_associated OR restriction_modification")

# Genomic neighborhood (with all annotation sources)
b.get_neighborhood(protein_id, window=10, all_annotations=True)

# Embedding similarity
similar = b.find_similar(protein_id, k=20)

# Structure prediction + remote homology
b.predict_structure(protein_id)
hits = b.search_foldseek_for_protein(protein_id)

# Export
b.export_fasta(protein_ids, "output.faa")
```

### Use with Claude Code

Bennu ships with skill specs in `.claude/skills/` that give Claude Code structured workflows for metagenomic analysis:

```
/survey    # Systematic comprehensive survey of a dataset
/explore   # Curiosity-driven hypothesis testing
/defense   # Defense system inventory
/metabolism # Metabolic pathway reconstruction
/literature # Literature search for functional claims
/characterize # Deep-dive on unknown proteins
```

## Architecture

```
┌─────────────────────────────────────────────────────────┐
│                   Agent (Claude Code, etc.)              │
│  Skills • Workflows • Multi-turn reasoning              │
└────────────────────────┬────────────────────────────────┘
                         │
┌────────────────────────┴────────────────────────────────┐
│                    Operator Layer                        │
│  search • navigate • similarity • export • structure    │
│  predicates • visualization • introspection             │
└────────────────────────┬────────────────────────────────┘
                         │
┌────────────────────────┴────────────────────────────────┐
│                     Data Layer                           │
│  DuckDB (proteins, annotations, loci, predicates)       │
│  LanceDB (ESM2 embeddings, similarity search)           │
└─────────────────────────────────────────────────────────┘
```

## Project structure

```
├── bennu/                 # Core package
│   ├── core/              # Data models, session state, types
│   ├── storage/           # DuckDB store, vector store, schema, migrations
│   ├── operators/         # Search, navigation, similarity, export, visualization
│   ├── predicates/        # Functional predicate system + PFAM/KEGG/CAZy/VOG mappings
│   └── reports/           # PDF report generation
├── src/ingest/            # Ingestion pipeline (stages 00-07)
├── scripts/               # Reusable CLI utilities
├── tests/                 # Unit and integration tests
├── .claude/skills/        # Claude Code skill specifications
├── CLAUDE.md              # Agent knowledge base (protocols, patterns, tools)
└── pyproject.toml
```

## Predicate system

The predicate system is what makes Bennu more than a database wrapper. Annotations are mapped to functional predicates via curated rules:

- **PFAM**: 2000+ domain-to-predicate mappings + regex patterns
- **KEGG**: KO-to-predicate mappings for metabolic functions
- **CAZy**: Carbohydrate-active enzyme families
- **VOGdb**: Viral orthologous groups
- **Computed**: `giant` (>1000 aa), `unannotated` (no hits), `membrane_protein` (TM helices)

This lets agents ask functional questions ("find electron-bifurcating hydrogenases") instead of remembering accession numbers.

## Ingest pipeline

| Stage | Tool | Output |
|-------|------|--------|
| 00 | Prepare inputs | Organized genome/protein files |
| 01 | QUAST | Assembly quality metrics |
| 02 | DFAST (optional) | QC and taxonomic classification |
| 03 | Prodigal | Gene calling (.faa, .gff) |
| 04 | Astra | PFAM, KEGG, HydDB, DefenseFinder annotations |
| 04 (opt-in) | Astra + VOGdb | Viral orthologous groups (`--databases ... VOGdb`) |
| 05a | GECCO | Biosynthetic gene clusters |
| 05b | dbCAN | CAZyme annotations |
| 05c | minced | CRISPR arrays |
| 06 | ESM2 | Protein embeddings (320-dim) |
| 07 | Builder | DuckDB knowledge base + predicates |

Stage 07 also runs **hydrogenase subgroup classification** when HydDB annotations are present: DIAMOND search against the HydDB reference database assigns NiFe/FeFe subgroups (e.g., Group 1a, Group 4e). All classified hits receive subgroup predicates (`nife_group1`, `fefe_groupB`, etc.). Hits lacking PFAM corroboration (including all Group 4 NiFe) are tagged `hyddb_needs_curation` for agent-level neighborhood validation. See `scripts/classify_hydrogenases.py`.

## Development

```bash
pip install -e ".[dev]"
pytest tests/ --override-ini addopts=""
```

## Key documents

| Document | Purpose |
|----------|---------|
| [`CLAUDE.md`](CLAUDE.md) | Agent knowledge base -- tools, patterns, protocols |
| [`QUICKSTART.md`](QUICKSTART.md) | Step-by-step dataset ingestion guide |
| [`QUICK_REFERENCE.md`](QUICK_REFERENCE.md) | SQL patterns and operator cheatsheet |
| [`DATA_ORGANIZATION.md`](DATA_ORGANIZATION.md) | Data directory conventions |

## Citation

```bibtex
@software{bennu2025,
  author = {West-Roberts, Jacob},
  title = {Bennu: Agent-driven metagenomic discovery},
  year = {2025},
  url = {https://github.com/jwestrob/Sharur}
}
```

## License

MIT
