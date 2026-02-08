# 🪶 Bennu

> Agent-driven metagenomic discovery system

Bennu is an LLM-powered interface for exploring metagenomic datasets. Ask natural language questions and get structured, reproducible answers.

## Features

- **Natural language queries**: "Find hydrogenases in Archaea", "What's weird in Bin_023?"
- **Multi-turn exploration**: Build on previous results, save working sets, track hypotheses
- **Multiple search modes**: Annotations, taxonomy, embedding similarity, spatial proximity
- **Reproducible**: Full provenance tracking, export to notebooks

## Quick Start

### Install

```bash
git clone https://github.com/jacobwestroberts/bennu
cd bennu
pip install -e ".[dev]"
```

### Run an agent (DuckDB)

Python:
```python
from bennu import BennuAgent, ExplorationSession

session = ExplorationSession(db_path="data/bennu.duckdb")
agent = BennuAgent(session)  # heuristic routing if no LM configured
print(agent.process("Find proteins with PF00142"))
```

CLI (one-shot ask):
```bash
bennu ask --db data/bennu.duckdb "Find proteins with PF00142"
# or equivalently
python -m bennu.cli ask "Find proteins with PF00142" --db data/bennu.duckdb
```
Requires `OPENAI_API_KEY` in the environment (optional `BENNU_LM_MODEL`, default `gpt-5-mini-2025-08-07`).

### Build the knowledge base (ingest)

Pipeline runner (tools by default):
```bash
python scripts/ingest.py --input-dir dummy_dataset --data-dir data --output data/bennu.duckdb --force
```
*DFAST is optional; CRISPR uses `minced` if available. Embeddings and vector store are built by default.*

Stage 07 builder (manual):
```bash
python -m src.ingest.07_build_knowledge_base --data-dir data --output data/bennu.duckdb --force
```

Expected stage dirs: `stage02_dfast_qc` (optional), `stage03_prodigal`, `stage04_astra`, `stage05a_gecco`, `stage05b_dbcan`, `stage05c_crispr` (may be empty), `stage06_embeddings`.

## Installation

```bash
# From source
git clone https://github.com/jacobwestroberts/bennu
cd bennu
pip install -e ".[dev]"
```

## Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                    Exploration Session                       │
│  Working sets • Focus stack • Hypotheses • Provenance       │
└─────────────────────────────────────────────────────────────┘
                              │
┌─────────────────────────────────────────────────────────────┐
│                      DSPy Agent                              │
│  Router → Parameter Extraction → Execution → Synthesis      │
└─────────────────────────────────────────────────────────────┘
                              │
┌─────────────────────────────────────────────────────────────┐
│                      Tool Layer                              │
│  find_proteins • get_context • detect_loci • find_similar   │
│  find_anomalies • compare_across • manage_sets • export     │
└─────────────────────────────────────────────────────────────┘
                              │
┌─────────────────────────────────────────────────────────────┐
│                      Data Layer                              │
│  DuckDB (relational) + LanceDB (vectors)                   │
└─────────────────────────────────────────────────────────────┘
```

## Key Documents

| Document | Purpose |
|----------|---------|
| [`CLAUDE.md`](CLAUDE.md) | **Agent knowledge base** - tools, patterns, protocols |
| [`QUICK_REFERENCE.md`](QUICK_REFERENCE.md) | Critical patterns and SQL snippets |
| [`DATA_ORGANIZATION.md`](DATA_ORGANIZATION.md) | Data directory structure |

## Development

```bash
# Install dev dependencies
pip install -e ".[dev]"

# Run tests
pytest

# Format code
black bennu/ tests/
ruff check --fix bennu/ tests/

# Type check
mypy bennu/
```

## Project Structure

```
bennu/
├── bennu/
│   ├── core/           # Data models, session state
│   ├── storage/        # DuckDB store, vector store, schema
│   ├── tools/          # DSPy agent tools
│   ├── agent/          # DSPy signatures and orchestrator
│   ├── operators/      # Bennu operators (search, navigate, visualize)
│   ├── predicates/     # Functional predicate system
│   └── reports/        # PDF report generation
├── tests/
├── scripts/            # CLI utilities
├── examples/           # Notebooks
└── CLAUDE.md           # Agent knowledge base
```

## Tools

| Tool | Description |
|------|-------------|
| `find_proteins` | Multi-modal search by domain, function, taxonomy, similarity |
| `get_genomic_context` | Neighborhood around a protein with ASCII visualization |
| `detect_loci` | Find prophages, BGCs, CRISPR arrays, operons |
| `find_similar` | Embedding similarity at protein or locus level |
| `find_anomalies` | Statistical outliers on pre-computed metrics |
| `compare_across` | Cross-genome feature comparisons |
| `manage_sets` | Create and manipulate working sets |
| `export` | Export to FASTA, GFF, TSV, JSON |

## Data Requirements

Bennu expects:
- **Proteins** with genomic coordinates and annotations
- **ESM2 embeddings** (320-dim) for similarity search

See [scripts/load_metagenome.py](scripts/load_metagenome.py) for data loading utilities.

## Citation

If you use Bennu in your research, please cite:

```bibtex
@software{bennu2024,
  author = {West-Roberts, Jacob},
  title = {Bennu: Agent-driven metagenomic discovery},
  year = {2024},
  url = {https://github.com/jacobwestroberts/bennu}
}
```

## License

MIT

---

*Named after the ancient Egyptian deity associated with creation and rebirth, and the OSIRIS-REx asteroid.*
