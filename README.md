<p align="center">
  <img src="assets/logo.png" alt="Sharur logo" width="480">
</p>

# Sharur

> *The weapon which loved the lord, obedient to its master — it spun around the horizon of heaven to find out what was happening, and joyfully it reported the message.*
>
> — [*Ninurta's exploits*](https://etcsl.orinst.ox.ac.uk/cgi-bin/etcsl.cgi?text=t.1.6.2&charenc=j#) (Lugal-e), Sumerian, c. 2100 BCE

In Sumerian mythology, Sharur is the sentient mace of the god Ninurta. It flies ahead to scout the unknown, gathers intelligence, and reports back what it finds.

> A data plane for agent-driven metagenomic discovery

Sharur makes large metagenomic datasets navigable by AI agents. It combines a DuckDB relational store, a FAISS vector store, and an evidence-backed functional predicate system into typed operators that agents use to search, characterize, and compare proteins across thousands of genomes.

Sharur is model-agnostic. Its interface is a CLI and a Python API whose outputs are bounded, typed, and sequence-free, so any agent that can run shell commands or Python can use it: **Claude Code**, **Codex**, or **open-weight models** running in an agent harness of your choice. Agents bring the reasoning; Sharur brings the data access and the guardrails.

## What it does

Given a set of metagenome-assembled genomes (MAGs), Sharur:

1. **Ingests** proteins, annotations (Pfam, KOfam, HydDB, VOGdb, CAZy, DefenseFinder, TXSScan), CRISPR arrays, biosynthetic gene clusters, and ESM2 embeddings into one database.
2. **Computes predicates**: functional tags such as `nife_group3`, `sam_binding` or `crispr_associated`. Every mapping carries recorded evidence (GO, Pfam/KEGG text, ENZYME, KEGG BRITE and modules, HydDB, reviewed-protein consensus), and system-level claims come only from validated system callers. See [How predicates are built](docs/predicate_construction.md).
3. **Exposes operators** for agents: predicate search, genomic neighborhoods, protein cards, predicate explanations, KEGG module completeness, embedding similarity, structure search, and export.
4. **Records provenance**: dataset seals, a capability preflight, and a stamp of the predicate maps behind every generation.

## Working with agents

| Agent | How it uses Sharur |
|---|---|
| **Claude Code** | Reads `CLAUDE.md` and the skills in `.claude/skills/` (`/survey`, `/explore`, `/characterize`, `/hydrogenase`, …); calls the CLI and Python API. |
| **Codex** | Reads `AGENTS.md` (the same instructions as `CLAUDE.md`); calls the CLI and Python API. |
| **Open-weight models** | Any agent harness that runs shell commands or Python can drive the same CLI and API. The bounded outputs of `sharur describe`, `sharur card`, and `sharur why` are designed for small context budgets. |

For large campaigns, Sharur's worker executors drive headless `claude -p` and `codex exec` sessions against coordinated Atlas tasks (see [Scaling up](#scaling-up)).

Rules every agent follows live in `CLAUDE.md`: report domain observations separately from named functions, name a system only when a curated caller called it, attach a verification query to every number, and keep sequences out of model-visible text.

## Quick start

```bash
git clone https://github.com/jwestrob/Sharur.git
cd Sharur
pip install -e ".[all,dev]"
sharur doctor          # which external tools, reference databases and keys are present
sharur setup-kegg      # build the KEGG predicate map locally (KEGG REST, academic use)
```

`pip install -e "."` is the lean core; focused extras include `parquet`, `vectors`, `embeddings`, `ops`, `visualization`, `structure`, `reports`, and `notebooks`. External tools (Prodigal, DIAMOND, HMMER, Astra, Foldseek, …) and reference databases install separately; see [`INSTALL.md`](INSTALL.md).

KEGG data is subject to [KEGG's terms](https://www.kegg.jp/kegg/legal.html): `sharur setup-kegg` builds the KO → predicate map on your machine, and non-academic users can build from a licensed KEGG copy with `--inputs`. Sharur ships only its own rules; see [`DATA_LICENSES.md`](DATA_LICENSES.md).

### Ingest a dataset

```bash
sharur-ingest \
  --input-dir /path/to/genome_fastas \
  --data-dir data/my_dataset \
  --output data/my_dataset/sharur.duckdb \
  --profile auto
```

`sharur-ingest` runs a ledger-backed, resumable stage DAG and builds DuckDB before launching embeddings. Optional stages are opt-in (`--with-quast`, `--with-dfast`, `--with-gecco`; the dbCAN three-tool consensus with `--enable-cazymes`). Choose `--profile local`, `mps`, or `slurm` when `auto` does not fit, and inspect the plan with `--dry-run`. See [`QUICKSTART.md`](QUICKSTART.md) and [`src/ingest/README.md`](src/ingest/README.md).

### First look

```bash
sharur describe --db data/my_dataset/sharur.duckdb        # sources, curated callers, genome metadata, predicate state
sharur card PROTEIN_ID --db data/my_dataset/sharur.duckdb # one protein: context, hits, evidence-backed predicates
sharur why PROTEIN_ID nad_binding --db data/my_dataset/sharur.duckdb   # the evidence behind one predicate
sharur modules --db data/my_dataset/sharur.duckdb --bin GENOME_ID --min-completeness 0.75
sharur preflight --db data/my_dataset/sharur.duckdb --format json      # typed capability brief
```

Seal a finished dataset and verify it later:

```bash
sharur seal --db data/my_dataset/sharur.duckdb
sharur verify-seal data/my_dataset/dataset.seal.json
```

### Use the operators

```python
from sharur.operators import Sharur

b = Sharur("data/my_dataset/sharur.duckdb", read_only=True)

b.describe()                                          # what the dataset holds
hydrogenases = b.search_by_predicates(has=["nife_group3", "bidirectional_hydrogenase"])
giants = b.search_by_predicates(has=["giant", "unannotated"])

b.card(protein_id)                                    # bounded protein summary
b.why(protein_id, "sam_binding")                      # evidence paths for one predicate
b.get_neighborhood(protein_id, window=10, all_annotations=True)
b.modules(bins=[genome_id], min_completeness=0.75)    # KEGG module completeness
b.locus_modules(protein_id, window=10)                # module steps co-encoded nearby

case = b.inspect(caller_or_protein_id, entity_type="system", upstream_orfs=4, downstream_orfs=8)
similar = b.find_similar(protein_id, k=20)            # embedding similarity
b.export_fasta(protein_ids, "output.faa")
```

Operator results expose `result.records`, `result.raw`, `result.status`, and `result.to_json()`. CLI listing commands accept `--format markdown|json|jsonl|tsv`. `sharur inspect` and `sharur compare-context` expose typed cases and foreground/background context tests; see [`docs/cases_and_evidence.md`](docs/cases_and_evidence.md).

## Core and extensions

**Core**: ingest, the predicate system, the operator API and CLI (`describe`, `card`, `why`, `modules`, search, neighborhoods), preflight and seals. This is everything a single analyst or agent needs.

**Extensions**, installed and used as needed:

- **Embeddings and similarity**: ESM2 embeddings and FAISS indexes (`embeddings`, `vectors` extras).
- **Synteny**: ELSA embedding-based conserved gene blocks, exposed through a run-scoped `synteny.duckdb` sidecar.
- **Structure**: ESM3 structure prediction and Foldseek remote homology (`structure` extra).
- **Campaign scale**: `sharur-ops`, `sharur-query`, `sharur-atlas`, `sharur-review` (below).
- **Reports and figures**: PDF reports and publication-quality neighborhood and domain figures.

## Scaling up

Coordinated campaigns over large databases route through three services: `sharur-ops` owns tasks, leases and findings; `sharur-query` serves one sealed, read-only DuckDB through bounded typed operators and a shared cache; `sharur-review` reduces candidate records through a hierarchical review DAG.

```bash
export SHARUR_OPS_URL=http://ops-host:8811
sharur-query --db data/my_dataset/sharur.duckdb --direct --host 0.0.0.0

sharur-atlas plan --db data/my_dataset/sharur.duckdb --output-dir data/my_dataset/atlas --packet-bytes 524288
sharur-atlas packet-census --plan-dir data/my_dataset/atlas
sharur-atlas enqueue --plan-dir data/my_dataset/atlas --ops-url http://ops-host:8811 --query-url http://query-host:8812

sharur-review reduce --ops-url http://ops-host:8811 --campaign-id CAMPAIGN_ID
sharur-review route --ops-url http://ops-host:8811 --campaign-id CAMPAIGN_ID --watch
```

Atlas gives each genome one task; packets combine whole consecutive contigs from that genome, and workers checkpoint cursors and write a coverage manifest, typed candidates and one unit disposition. See [`docs/query_service.md`](docs/query_service.md), [`docs/agent_ops_spec.md`](docs/agent_ops_spec.md) and [`docs/review_workflow.md`](docs/review_workflow.md).

```
┌────────────────────────────────────────────────────────────┐
│ Agents: skills • workflows • multi-turn reasoning          │
└───────────────┬────────────────────────┬───────────────────┘
                │ coordination           │ typed queries
                v                        v
┌───────────────────────────┐  ┌─────────────────────────────┐
│ sharur-ops                │  │ sharur-query                │
│ one SQLite control owner  │  │ one DuckDB owner and cache  │
└───────────────┬───────────┘  └──────────────┬──────────────┘
                v                             v
┌───────────────────────────┐  ┌─────────────────────────────┐
│ sharur_ops.db             │  │ Typed operator layer        │
│ tasks • leases • findings │  │ search • navigate • V2      │
└───────────────────────────┘  └──────────────┬──────────────┘
                                              v
                               ┌─────────────────────────────┐
                               │ Sealed local DuckDB replica │
                               │ FAISS / ELSA sidecars       │
                               └─────────────────────────────┘
```

## Predicate system

Annotations map to predicates through evidence-backed tables:

- **Pfam**: a generated family → predicate table in which every pair cites GO, the family's name/description, an ENZYME name, or reviewed-protein consensus.
- **KEGG**: a KO → predicate table built locally by `sharur setup-kegg` from KEGG's EC numbers, BRITE placements and modules, HydDB labels for hydrogenase KOs, and reviewed-protein consensus.
- **CAZy**: family predicates from CAZy class definitions and reviewed-protein consensus.
- **VOGdb**: functional categories and vetted consensus-description matches.
- **Computed**: size classes, annotation status, and transmembrane topology.

Predicates are either component-level (what one gene or domain shows) or system-level (a validated multi-gene system, emitted only by system callers). The V2 backend stores each claim as a typed atom with a facet, a relation (`implies` / `supports` / `flags`) and evidence metadata; composite rules build higher-order conclusions. [How predicates are built](docs/predicate_construction.md) lists every source, rule and threshold, and `sharur why` shows the chain for any protein. Maintainers rebuild all maps from a dated snapshot with `make snapshots predicate-maps recount`.

Proposals for new mappings are welcome: add them to the proposal files and rebuild; a pair ships once its evidence verifies.

## Ingest pipeline

| Stage | Tool | Output |
|-------|------|--------|
| 00 | Prepare inputs | Validate and organize genome FASTAs |
| 01 | QUAST | Optional assembly QC metrics |
| 02 | DFAST | Optional QC and taxonomy |
| 03 | Prodigal | Gene calling (`.faa`, `.genes.fna`) |
| 04 | Astra | Pfam, KOfam, HydDB, DefenseFinder, dbCAN annotation |
| 04 (opt-in) | Astra + extra DBs | VOGdb, TXSScan, CANT-HYD via repeated `-d` flags |
| 05a | GECCO | Optional biosynthetic gene clusters |
| 05c | minced | CRISPR array detection |
| 07 | Builder | DuckDB knowledge base and predicates; optional dbCAN consensus with `--enable-cazymes` |
| 06 | ESM2 | Protein embeddings (required for ELSA) |

Stage 07 also assigns hydrogenase subgroups: each HydDB-hit protein receives the subgroup of its nearest HydDB reference, recorded in `hydrogenase_classifications` with a catalytic-domain check and KOfam support. `scripts/classify_hydrogenases.py` refreshes an existing database through a validated staged copy.

## Project structure

```
├── sharur/                # Core package
│   ├── operators/         # Search, navigation, cards, similarity, export, visualization
│   ├── predicates/        # Vocabulary, Pfam/KEGG/CAZy/VOG mappings and evidence, provenance
│   ├── predicates_v2/     # Semantic-atom backend (atoms, composites, persistence)
│   ├── hydrogenase/       # HydDB subgroup assignment and KOfam support
│   ├── storage/           # DuckDB store, vector store, schema, migrations
│   ├── ingest/            # Packaged ingest stages
│   ├── ops/ query/ review/ workers/   # Campaign-scale services and worker executors
│   └── modules.py         # KEGG module completeness
├── config/predicates_v2/  # V2 facets, relations, composites
├── src/ingest/            # Stage-by-stage ingest reference
├── scripts/               # Map builders, callers, loaders, report renderers
├── docs/                  # Reference guides (routed from CLAUDE.md)
├── .claude/skills/        # Claude Code skills
├── CLAUDE.md / AGENTS.md  # Agent instructions (one file)
└── tests/
```

## Development

```bash
pip install -e ".[dev,vectors,ops,reports]"
make test                  # pytest
make lint
```

## Key documents

| Document | Purpose |
|----------|---------|
| [`CLAUDE.md`](CLAUDE.md) | Agent rules and routing table to `docs/` |
| [`QUICKSTART.md`](QUICKSTART.md) | `sharur-ingest` workflow for new datasets |
| [`docs/predicate_construction.md`](docs/predicate_construction.md) | How predicates are built |
| [`docs/biological_interpretation.md`](docs/biological_interpretation.md) | Annotation provenance and claim discipline |
| [`docs/predicates_v2.md`](docs/predicates_v2.md) | V2 semantic-atom predicate system |
| [`docs/findings_spec.md`](docs/findings_spec.md) | Structured, verifiable findings |
| [`docs/tools_reference.md`](docs/tools_reference.md) | Astra, ELSA, ESM3, Foldseek |
| [`QUICK_REFERENCE.md`](QUICK_REFERENCE.md) | SQL patterns and operator cheatsheet |
| [`DATA_LICENSES.md`](DATA_LICENSES.md) | Licenses of shipped data and locally built KEGG data |

## Citation

See [`CITATION.cff`](CITATION.cff) and [`CITATIONS.md`](CITATIONS.md) for the tools and databases Sharur builds on.

```bibtex
@software{sharur,
  author = {West-Roberts, Jacob},
  title = {Sharur: agent-driven exploration of metagenomic datasets},
  year = {2026},
  url = {https://github.com/jwestrob/Sharur}
}
```

## License

MIT for the code; see [`DATA_LICENSES.md`](DATA_LICENSES.md) for shipped data.
