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
2. **Computes predicates** -- functional tags derived from annotation combinations (e.g., `nife_group3`, `crispr_associated`, `giant_unannotated`) that make semantic search possible. The current backend is **Predicate V2**, a typed semantic-atom system (see below)
3. **Exposes operators** that agents call to explore the data: search by predicate, navigate genomic neighborhoods, find similar proteins by embedding, detect loci, export results
4. **Maps synteny** with ELSA -- embedding-based conserved gene-block
   discovery over a per-dataset index, exposed through a normalized,
   run-scoped `synteny.duckdb` sidecar and exact Sharur operators

Agents bring the reasoning; Sharur brings the data access.

## Quick start

```bash
git clone https://github.com/jwestrob/Sharur.git
cd Sharur
pip install -e ".[all,dev]"
```

### Verify your toolchain

```bash
sharur doctor        # report which external tools, reference DBs, and API keys are present
sharur --version
```

The command above installs every optional Python surface for a development checkout.
`pip install -e "."` is the lean query/ingest-orchestration core; focused extras include
`parquet`, `vectors`, `embeddings`, `ops`, `visualization`, `structure`, `reports`, and
`notebooks`.
External bioinformatics tools (Prodigal, DIAMOND, HMMER, Astra, Foldseek, …) and reference
databases install separately — see [`INSTALL.md`](INSTALL.md). `sharur doctor` shows exactly
what's present and what's still missing.

### Ingest a dataset

```bash
sharur-ingest \
  --input-dir /path/to/genome_fastas \
  --data-dir data/my_dataset \
  --output data/my_dataset/sharur.duckdb \
  --profile auto
```

`sharur-ingest` is the primary command-line interface for the standard ingest pipeline. Use it by default for new datasets. See [`QUICKSTART.md`](QUICKSTART.md) for the operator-facing guide and [`src/ingest/README.md`](src/ingest/README.md) for the manual stage-by-stage reference. [`scripts/ingest.py`](scripts/ingest.py) remains as a repo-local compatibility shim for the same implementation. If `sharur-ingest` is not on `PATH`, refresh the editable install with `pip install -e ".[embeddings]"`.

The default plan runs the core stages and builds DuckDB before launching embeddings. Optional
QC/BGC stages are opt-in with `--with-quast`, `--with-dfast`, or `--with-gecco`;
the deprecated helper is available only as `--with-legacy-dbcan`. The slower dbCAN
three-tool consensus classifier is a separate Stage 07 opt-in,
`--enable-cazymes`. Ingest is a ledger-backed dependency DAG and resumes matching successful
stages by default. Embedding inference and the CPU FAISS build are separate attempts, so an
index failure never forces re-embedding or imports FAISS into the Torch process. Choose
`--profile local`, `--profile mps`, or `--profile slurm` explicitly when `auto` is not
appropriate. Inspect the exact plan with `--dry-run`.

After ingest, get one machine-readable view of what the dataset and runtime can actually do:

```bash
sharur preflight --db data/my_dataset/sharur.duckdb --format json
```

This inspects live sources/caller resources, semantic coverage, persistent similarity
sidecars, the run ledger, resource profiles, and tool availability without mutating data.

Once dataset writes are complete, record and later verify its canonical state:

```bash
sharur seal --db data/my_dataset/sharur.duckdb
sharur verify-seal data/my_dataset/dataset.seal.json
```

The default seal is disk-light and samples large canonical artifacts. Use `--full` when
making an archival seal that should stream full SHA-256 over DuckDB, embeddings, and active
index sidecars.

For a coordinated campaign over a large database, serve one sealed immutable
database through the bounded query data plane. Direct mode is appropriate
when a replica would occupy the same storage tier:

```bash
export SHARUR_OPS_URL=http://ops-host:8811
sharur-query --db data/my_dataset/sharur.duckdb --direct --host 0.0.0.0
```

Agents reuse their Sharur Ops credentials and call typed operators through
`SharurQuery`. One service owns the DuckDB cache, threads, memory, and spill
budget. See [`docs/query_service.md`](docs/query_service.md).

Exhaustive genome reading uses deterministic Atlas units over the same data
plane:

```bash
sharur migrate --db data/my_dataset/sharur.duckdb
sharur seal --db data/my_dataset/sharur.duckdb --force

sharur-atlas plan --db data/my_dataset/sharur.duckdb \
  --output-dir data/my_dataset/atlas \
  --packet-contigs 128 --packet-proteins 500 --packet-bytes 524288
sharur-atlas packet-census --plan-dir data/my_dataset/atlas
sharur-atlas verify-packet-census --plan-dir data/my_dataset/atlas --deep
sharur-atlas enqueue --plan-dir data/my_dataset/atlas \
  --ops-url http://ops-host:8811 --query-url http://query-host:8812
```

Each logical task owns one genome. Every model packet contains data from that
single bin, combines whole consecutive contigs when they fit, and splits only
oversized contigs. Enqueue requires the zero-model-call invocation/payload
census. Workers checkpoint packet cursors and write a coverage manifest plus
typed candidates and one reconciled unit disposition. Reduce and route those
records through the hierarchical review DAG:

```bash
sharur-review reduce --ops-url http://ops-host:8811 \
  --campaign-id CAMPAIGN_ID
sharur-review route --ops-url http://ops-host:8811 \
  --campaign-id CAMPAIGN_ID --watch
```

See [`docs/review_workflow.md`](docs/review_workflow.md).

### Use the operators

```python
from sharur.operators import Sharur

b = Sharur("data/my_dataset/sharur.duckdb", read_only=True)

# Predicate search
hydrogenases = b.search_by_predicates(has=["nife_group3", "bidirectional_hydrogenase"])
giants = b.search_by_predicates(has=["giant", "unannotated"])
defense = b.search_by_predicates(has=["crispr_associated"])

# Genomic neighborhood (with all annotation sources)
b.get_neighborhood(protein_id, window=10, all_annotations=True)

# First-class, provenance-separated case with asymmetric ORF context
case = b.inspect(
    caller_or_protein_id,
    entity_type="system",
    upstream_orfs=4,
    downstream_orfs=8,
)
comparison = case.compare_context(
    features=["pfam:PF00589", "other_called_system"],
    window=8,
    deduplicate_by="replicon",
)

# Embedding similarity
similar = b.find_similar(protein_id, k=20)

# Structure prediction + remote homology
b.predict_structure(protein_id)
hits = b.search_foldseek_for_protein(protein_id)

# Export
b.export_fasta(protein_ids, "output.faa")
```

Operator results expose `result.records`, `result.raw`, `result.status`, and
`result.to_json()` for programmatic use. The CLI commands `overview`, `genomes`,
`proteins`, `neighborhood`, and `search` accept
`--format markdown|json|jsonl|tsv`.

`sharur inspect` and `sharur compare-context` expose typed cases and exact
foreground/background context tests. Optional assembly/host evidence lives in
a separate sidecar, and composition scanning is available only through the
explicit `compute-composition-evidence` command. See
[`docs/cases_and_evidence.md`](docs/cases_and_evidence.md).

### Use with Claude Code

Sharur ships with skill specs in `.claude/skills/` that give Claude Code structured workflows for metagenomic analysis:

```
# Analysis
/survey       # Systematic comprehensive survey of a dataset
/explore      # Curiosity-driven hypothesis testing
/characterize # Deep-dive on unknown proteins
/compare      # Cross-genome comparative analysis
/atlas        # Exhaustive genome-by-genome reading
/metabolism   # Metabolic pathway reconstruction
/pathway      # KEGG pathway completeness check

# Domain specialists
/defense      # Defense system inventory (MacSyFinder-validated)
/prophage     # Prophage & viral element detection
/hydrogenase  # NiFe/FeFe hydrogenase subgroup validation
/synteny      # ELSA-powered conserved gene-block discovery
/foldseek     # ESM3 structure prediction + Foldseek homology
/visualize    # Publication-quality neighborhood/domain figures

# Reasoning & orchestration
/literature   # Literature search for functional claims
/reviewer_2   # Adversarial manuscript claim verification
/brainstorm   # Cross-domain synthesis + investigation proposals
/coordinator  # Multi-agent analysis orchestration
/query        # Fast ad-hoc database queries
```

## Architecture

```
┌────────────────────────────────────────────────────────────┐
│ Agents: skills • workflows • multi-turn reasoning          │
└───────────────┬────────────────────────┬───────────────────┘
                │ coordination           │ typed queries
                v                        v
┌───────────────────────────┐  ┌─────────────────────────────┐
│ sharur-ops                │  │ sharur-query                │
│ one SQLite control owner  │  │ one DuckDB owner and cache │
└───────────────┬───────────┘  └──────────────┬──────────────┘
                │                             │
                v                             v
┌───────────────────────────┐  ┌─────────────────────────────┐
│ sharur_ops.db             │  │ Typed operator layer        │
│ tasks • leases • findings │  │ search • navigate • V2     │
└───────────────────────────┘  └──────────────┬──────────────┘
                                              │
                                              v
                               ┌─────────────────────────────┐
                               │ Sealed local DuckDB replica │
                               │ FAISS / ELSA sidecars       │
                               └─────────────────────────────┘
```

## Project structure

```
├── sharur/                # Core package
│   ├── core/              # Data models, session state, types
│   ├── storage/           # DuckDB store, vector store, schema, migrations
│   ├── operators/         # Search, navigation, similarity, export, visualization
│   ├── predicates/        # V1 predicate system + PFAM/KEGG/CAZy/VOG mappings
│   ├── predicates_v2/     # V2 semantic-atom backend (atoms, composites, review queue)
│   ├── ops/               # Multi-agent coordination server + client
│   ├── query/             # Bounded shared DuckDB query server + client
│   └── reports/           # PDF report generation
├── config/predicates_v2/  # V2 YAML config: facets, relations, composites
├── src/ingest/            # Ingestion pipeline (stages 00-07)
├── scripts/               # Reusable CLI utilities
├── tests/                 # Unit and integration tests
├── docs/                  # On-demand reference guides (routed from CLAUDE.md)
├── .claude/skills/        # Claude Code skill specifications
├── CLAUDE.md              # Agent knowledge base (protocols, patterns, tools)
└── pyproject.toml
```

## Predicate system

The predicate system is what makes Sharur more than a database wrapper. Annotations are mapped to functional predicates via curated rules:

- **PFAM**: domain-to-predicate mappings (~1,900 rows in `sharur/predicates/mappings/data/pfam_predicates.tsv`) + regex patterns
- **KEGG**: KO-to-predicate mappings for metabolic functions
- **CAZy**: Carbohydrate-active enzyme families
- **VOGdb**: Viral orthologous groups
- **Computed**: `giant` (>1000 aa), `unannotated` (no hits), `membrane_protein` (TM helices)

This lets agents ask functional questions ("find electron-bifurcating hydrogenases") instead of remembering accession numbers.

### Predicate V2 (current backend)

New Stage 07 builds use **Predicate V2**, a typed semantic-atom system. Instead of flat booleans, each atom is a structured claim carrying **facet** (activity / role / architecture / localization / topology / size_class / quality_flag), **relation** (how strongly the evidence implies/supports/flags the claim), and evidence metadata. Composite YAML rules combine atoms into higher-order conclusions, and unmapped accessions are surfaced in a frequency-ranked review queue. For backward compatibility, Stage 07 still materializes the legacy flat `protein_predicates` table from V2 output, so `search_by_predicates`, reports, and older scripts work unchanged. See [`docs/predicates_v2.md`](docs/predicates_v2.md) and `config/predicates_v2/`.

The predicate system is a project very much in progress; there are tens of thousands of accessions that need to be described. Pull requests to add more predicate tags are more than welcome; they require careful attention and review. Please notify if you find any errors in classification.

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
| 07 | Builder | Standard DuckDB knowledge base + predicates; optional dbCAN consensus with `--enable-cazymes` |
| 06 | ESM2 | Standard post-pipeline protein embeddings (required for ELSA) |

Stage 07 also runs **hydrogenase subgroup classification** when HydDB annotations are present: DIAMOND search against the HydDB reference database assigns NiFe/FeFe subgroups (e.g., Group 1a, Group 4e). All classified hits receive subgroup predicates (`nife_group1`, `fefe_groupB`, etc.). Hits lacking PFAM corroboration (including all Group 4 NiFe) are tagged `hyddb_needs_curation` for agent-level neighborhood validation. See `scripts/classify_hydrogenases.py`.

## Development

```bash
pip install -e ".[dev,vectors,ops,reports]"
pytest tests/ --override-ini addopts=""
```

## Key documents

| Document | Purpose |
|----------|---------|
| [`CLAUDE.md`](CLAUDE.md) | Agent knowledge base -- core rules + routing table to `docs/` |
| [`QUICKSTART.md`](QUICKSTART.md) | Primary `sharur-ingest` workflow for new datasets |
| [`src/ingest/README.md`](src/ingest/README.md) | Manual stage-by-stage ingest reference |
| [`QUICK_REFERENCE.md`](QUICK_REFERENCE.md) | SQL patterns and operator cheatsheet |
| [`DATA_ORGANIZATION.md`](DATA_ORGANIZATION.md) | Data directory conventions |
| [`docs/predicates_v2.md`](docs/predicates_v2.md) | Predicate V2 semantic-atom system |
| [`docs/tools_reference.md`](docs/tools_reference.md) | Astra, ELSA, ESM3, Foldseek, V2 atoms |
| [`docs/findings_spec.md`](docs/findings_spec.md) | Canonical schema for structured, verifiable findings |
| [`docs/biological_interpretation.md`](docs/biological_interpretation.md) | Annotation provenance and claim discipline |
| [`docs/analysis_workflow.md`](docs/analysis_workflow.md) | Full 5-phase analysis pipeline |
| [`docs/manuscript_guide.md`](docs/manuscript_guide.md) | Manuscript and report compilation |
| [`agent_ops_spec.md`](agent_ops_spec.md) | Multi-agent coordination layer (`sharur/ops/`) |

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
