<p align="center">
  <img src="assets/logo.png" alt="Sharur logo" width="480">
</p>

# Sharur

> *The weapon which loved the lord, obedient to its master — it spun around the horizon of heaven to find out what was happening, and joyfully it reported the message.*
>
> — [*Ninurta's exploits*](https://etcsl.orinst.ox.ac.uk/cgi-bin/etcsl.cgi?text=t.1.6.2&charenc=j#) (Lugal-e), Sumerian, c. 2100 BCE

In Sumerian mythology, Sharur is the sentient mace of the god Ninurta. It flies ahead to scout the unknown, gathers intelligence, and reports back what it finds.

> A data plane for agent-driven metagenomic discovery

Sharur makes large metagenomic datasets navigable by AI agents. It combines a DuckDB relational store, FAISS vector store, and a functional predicate system into an operator framework that agents (Claude Code, Codex, etc.) use to search, characterize, and compare proteins across hundreds of genomes.

## What it does

Given a set of metagenome-assembled genomes (MAGs), Sharur:

1. **Ingests** proteins, annotations (PFAM, KEGG, HydDB, VOGdb, CAZy, DefenseFinder), CRISPR arrays, biosynthetic gene clusters, and ESM2 embeddings into a unified database
2. **Computes predicates** -- functional tags derived from annotation combinations (e.g., `nife_group3`, `crispr_associated`, `giant_unannotated`) that make semantic search possible
3. **Exposes operators** that agents call to explore the data: search by predicate, navigate genomic neighborhoods, find similar proteins by embedding, detect loci, export results

Agents bring the reasoning; Sharur brings the data access.

## Quick start

```bash
git clone https://github.com/jwestrob/Sharur.git
cd Sharur
pip install -e ".[dev]"
```

### Ingest a dataset

```bash
sharur-ingest \
  --input-dir /path/to/genome_fastas \
  --data-dir data/my_dataset \
  --output data/my_dataset/sharur.duckdb \
  --force
```

`sharur-ingest` is the primary command-line interface for the standard ingest pipeline. Use it by default for new datasets. See [`QUICKSTART.md`](QUICKSTART.md) for the operator-facing guide and [`src/ingest/README.md`](src/ingest/README.md) for the manual stage-by-stage reference. [`scripts/ingest.py`](scripts/ingest.py) remains as a repo-local compatibility shim for the same implementation. If `sharur-ingest` is not on `PATH`, refresh the editable install with `pip install -e ".[dev]"`.

### Use the operators

```python
from sharur.operators import Sharur

b = Sharur("data/my_dataset/sharur.duckdb")

# Predicate search
hydrogenases = b.search_by_predicates(has=["nife_group3", "bidirectional_hydrogenase"])
giants = b.search_by_predicates(has=["giant", "unannotated"])
defense = b.search_by_predicates(has=["crispr_associated"])

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

Sharur ships with skill specs in `.claude/skills/` that give Claude Code structured workflows for metagenomic analysis:

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
│  FAISS (ESM2 embeddings, similarity search)           │
└─────────────────────────────────────────────────────────┘
```

## Project structure

```
├── sharur/                # Core package
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

The predicate system is what makes Sharur more than a database wrapper. Annotations are mapped to functional predicates via curated rules:

- **PFAM**: 2000+ domain-to-predicate mappings + regex patterns
- **KEGG**: KO-to-predicate mappings for metabolic functions
- **CAZy**: Carbohydrate-active enzyme families
- **VOGdb**: Viral orthologous groups
- **Computed**: `giant` (>1000 aa), `unannotated` (no hits), `membrane_protein` (TM helices)

This lets agents ask functional questions ("find electron-bifurcating hydrogenases") instead of remembering accession numbers.

The predicate system is a project very much in progress; there are tens of thousands of accessions that need to be described. Around 2-3,000 have been tagged so far. Pull requests to add more predicate tags to our system are more than welcome; they require careful attention and review. Please notify if you find any errors in classification among the predicate system. 

## Ingest pipeline

| Stage | Tool | Output |
|-------|------|--------|
| 00 | Prepare inputs | Standard. Validate and organize genome FASTAs |
| 01 | QUAST | Optional assembly QC metrics |
| 02 | DFAST | Optional QC and taxonomy |
| 03 | Prodigal | Standard gene calling (`.faa`, `.genes.fna`) |
| 04 | Astra | Standard PFAM, KOFAM, HydDB, DefenseFinder, dbCAN annotation |
| 04 (opt-in) | Astra + extra DBs | Optional VOGdb, TXSScan, CANT-HYD via repeated `-d` flags |
| 05a | GECCO | Optional biosynthetic gene clusters |
| 05b | Legacy dbCAN | Deprecated; standard CAZyme ingest comes from Stage 04 + Stage 07 |
| 05c | minced | Standard CRISPR array detection |
| 07 | Builder | Standard DuckDB knowledge base + predicates |
| 06 | ESM2 | Standard post-pipeline protein embeddings (required for ELSA) |

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
| [`QUICKSTART.md`](QUICKSTART.md) | Primary `sharur-ingest` workflow for new datasets |
| [`src/ingest/README.md`](src/ingest/README.md) | Manual stage-by-stage ingest reference |
| [`QUICK_REFERENCE.md`](QUICK_REFERENCE.md) | SQL patterns and operator cheatsheet |
| [`DATA_ORGANIZATION.md`](DATA_ORGANIZATION.md) | Data directory conventions |

## Citation

```bibtex
@software{sharur2025,
  author = {West-Roberts, Jacob},
  title = {Sharur: Agent-driven metagenomic discovery},
  year = {2025},
  url = {https://github.com/jwestrob/Sharur}
}
```

## License

MIT
