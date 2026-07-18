# Tools Reference

**Load this doc when:** Using Astra, ESM2/ESM3, Foldseek, ELSA, CAZyme/HydDB pipelines, V2 predicates, or hypothesis tracking.

## Astra Annotation Pipeline

Astra manages pre-installed HMM databases for annotation searches.

**Location:** `~/astra/` (source), installed via pyenv shim

**Installed databases** (`~/.config/Astra/`): PFAM, KOFAM, VOGdb, HydDB, DefenseFinder, CRISPRCasFinder, CANT-HYD, dbCAN, padloc, TXSScan

**NEVER call `astra search` directly.** Use the stage 04 script, which knows the correct flags for each database:

```bash
python src/ingest/04_astra_scan.py \
    -i data/DATASET/stage03_prodigal \
    -o data/DATASET/stage04_astra \
    -d KOFAM -t 10 --force
```

`--force` only removes the specific database subdirectory being re-run, not other results. Multiple databases: `-d PFAM -d KOFAM -d HydDB`.

**Why not raw `astra search`?** Each database needs different flags:
- PFAM/HydDB/DefenseFinder: `--cut_ga`
- KOFAM: `--cut_ga --cascade` (adaptive per-profile thresholds)
- DefenseFinder/TXSScan: `--write_macsyfinder` (for system validation)
- VOGdb/CANT-HYD: no cutoffs (filtered at load time)

The stage 04 script handles all of this automatically.

**Notes:**
- `--prot_in` expects a **directory** containing `.faa` files, not a single file
- Output: consolidated hits at `<outdir>/<database_lowercase>_results/<DATABASE>_hits_df.tsv` (per-genome intermediates land in `tmp_results/*_results.tsv`)
- For single files: `mkdir source/ && cp proteins.faa source/`

### Secretion System Identification (TXSScan)
**Validation:** Use the authoritative co-location engine below. Do not reconstruct systems
from raw profile rows.
**Requires:** TXSScan HMMs via Astra (`-d TXSScan` in stage 04), TXSScan models in `~/.macsyfinder/models/`
**Output:** replicon-provenanced `secretion_systems`, normalized `system_proteins`, and
`txsscan_system` annotations in DuckDB. Inspect the live caller table for the exact systems
emitted by the currently installed model definitions.

## Defense System Validation (Co-location Engine)

**Module:** `sharur/colocation.py`
**Purpose:** Calls systems from raw DefenseFinder/TXSScan profile hits using the installed
MacSyFinder XML definitions and ordered-replicon semantics, in process against DuckDB.

`scripts/validate_defense_systems.py` remains only as a compatibility CLI and delegates to
this module; it no longer runs an independent concatenated-FASTA workflow.

### Usage

```python
import logging
logging.basicConfig(level=logging.INFO, format='%(message)s')

from sharur.colocation import validate_systems, integrate_defense_results, integrate_secretion_results

db_path = "data/DATASET/sharur.duckdb"
affected = set()

# DefenseFinder
systems_df, genes_df = validate_systems(db_path, source="defensefinder")
affected |= integrate_defense_results(db_path, systems_df, genes_df)

# TXSScan (if loaded)
systems_df, genes_df = validate_systems(db_path, source="txsscan")
affected |= integrate_secretion_results(db_path, systems_df, genes_df)

# Integration deliberately invalidates derived semantics for old/new members.
from sharur.operators import Sharur
b = Sharur(db_path)
b.generate_v2(
    protein_ids=sorted(affected),
    update_legacy_predicates=True,
    predict_topology=False,
)
```

### What it does
1. Parses the currently installed model definitions and profile-name metadata.
2. Quarantines ambiguous internal HMM names that cannot be mapped to one model component.
3. Selects one deterministic global best profile per `(bin, contig, gene_index)`.
4. Clusters and evaluates each composite `(bin_id, contig_id)` replicon independently using
   model- and gene-specific gaps, quorum, forbidden, multi-locus, loner, and multi-system
   rules.
5. Uses an exact deterministic conflict solver with MacSyFinder-compatible scoring and
   sharing rules.
6. Transactionally replaces the structured table, normalized membership, and validated
   annotation surface, including when the new result is empty.

### Caller validation

On the full DPANN translated TXSScan hit set, MacSyFinder 2.1.4 / TXSScan 1.1.3 produced
904 candidates and selected 687 systems across all 19 installed models. Sharur selected
exactly the same 687 instances (zero disagreements), with identical output under three
Python hash seeds.

### Requirements
- MacSyFinder model XMLs at `~/.macsyfinder/models/defense-finder-models/` (installed by defense-finder)
- Raw HMM hits loaded in DuckDB `annotations` table with `source='defensefinder'`
- Astra `--write_macsyfinder` is NOT needed (engine reads from DuckDB, not hmmsearch files)

Schema-v5 migration quarantines pre-v5 named DefenseFinder/TXSScan calls and invalidates
their affected semantic cache rows. Rerun Stage 07 (or the caller plus a V2 subset refresh)
before treating a migrated legacy database's system surface as available.

## Hydrogenase Classification
**Script:** `scripts/classify_hydrogenases.py`
**Requires:** HydDB HMMs via Astra, DIAMOND database (`data/reference/hyddb/HydDB_all.dmnd`)
**Output:** Subgroup-level classification (NiFe Group 1-4, FeFe A-C)
**Note:** Hits lacking PFAM corroboration are tagged `hyddb_needs_curation` for agent neighborhood curation.
**Full protocol:** `.claude/skills/hydrogenase.md`

## CAZyme Classification (dbCAN 3-tool consensus)
**Script:** `scripts/classify_cazymes.py`
**Requires:** `data/dbcan_db/` with CAZy.dmnd, dbCAN.hmm, dbCAN-sub.hmm
**Method:** DIAMOND (1e-18) + dbCAN.hmm (1e-15) + dbCAN-sub.hmm (1e-15), consensus ≥2 tools.
**Pipeline:** Optional in Stage 07 before final V2 predicate generation. Enable it with
`sharur-ingest --enable-cazymes` or standalone Stage 07 `--enable-cazymes`.

## ESM2 Embeddings
**Script:** `src/ingest/06_esm2_embeddings.py`
**Usage:** `python src/ingest/06_esm2_embeddings.py data/DATASET/stage03_prodigal data/DATASET/embeddings/`
**Model:** `facebook/esm2_t6_8M_UR50D` (320-dim, auto-detects MPS/CUDA/CPU)
**Output:** canonical `embeddings/protein_embeddings.h5`, `embedding_manifest.json`,
an atomic `protein_embeddings.index.json`, a generation-scoped `.faiss` index, and a
generation-scoped `.ids.sqlite` stable row-to-protein map. Stage 06 builds the sidecars by
default in a Torch-free subprocess; use `--skip-index` only when deliberately deferring the
CPU build. The standard ingest DAG deliberately does defer it and records the follow-on CPU
build as `06i`, allowing index retry without another model run.

FASTA input and H5 output are streamed in batches. Residue pooling excludes padding and
special tokens, invalid/zero embeddings fail closed, the H5 is published atomically, and the
embedding manifest counts sequences truncated at the model limit.

Sharur discovers artifact paths at session creation but does not open H5 or FAISS until the
first similarity call, so ordinary DuckDB operators stay lightweight. A valid generation is
opened read-only with mmap. Missing or stale sidecars are rebuilt on first use, or explicitly:

```bash
sharur build-vector-index --db data/DATASET/sharur.duckdb
```

The H5 remains canonical. The index manifest is the atomic commit point, and its source
signature prevents reuse after the H5 changes.
**Status:** Standard pipeline stage — run after stage 07. Required for ELSA synteny.

## ELSA Synteny Discovery

**Location:** External package (`pip install -e /path/to/ELSA`)
**Purpose:** Embedding-based synteny detection across genomes.
**Requires:** `faiss-cpu`, `h5py`, `duckdb`. Set `KMP_DUPLICATE_LIB_OK=TRUE` on macOS.

**Data contract — ELSA reads Sharur outputs directly:**
1. `sharur.duckdb` — reads `proteins` table
2. `embeddings/protein_embeddings.h5` — standard embedding format
3. `sharur.duckdb` `annotations` table (via `--annotations-db`)

**Basic CLI:**
```bash
elsa synteny \
    --db data/DATASET/sharur.duckdb \
    --embeddings data/DATASET/embeddings/protein_embeddings.h5 \
    --annotations-db data/DATASET/sharur.duckdb \
    --store data/DATASET/synteny/store \
    -o data/DATASET/synteny/
```

**For complete query patterns, gene ID formats, and citing rules, see `.claude/skills/synteny.md`.**

## Embedding Visualization
**Script:** `scripts/visualize_embeddings.py`
**Usage:** `python scripts/visualize_embeddings.py --db data/DATASET/sharur.duckdb --output figures/umap.html --color-by genome`
**Color modes:** `genome`, `predicate`, `annotation`

## Local Foldseek
**Binary:** Auto-detected via `which foldseek`
**Database path:** `~/.foldseek/{db_name}/{db_name}`
**Behavior:** Tries local binary first, falls back to web API.

## Structure Prediction & Foldseek

```python
b = Sharur("data/DATASET/sharur.duckdb", read_only=True)

# Database proteins
result = b.predict_structure("protein_id", output_path="structures/protein.pdb")  # requires ESM_API_KEY

# Arbitrary sequences
from sharur.operators.structure import predict_structure_from_sequence
result = predict_structure_from_sequence(sequence="MKVL...", output_path="output.pdb", name="my_protein")

# Batch
result = b.batch_predict_structures(protein_ids=["id1", "id2"], output_dir="structures/", max_length=1024)

# Foldseek search
hits = b.search_foldseek("structures/protein.pdb", databases=["afdb50", "pdb100"], top_k=10)
hits = b.search_foldseek_for_protein("protein_id")  # convenience
```

**TM-score:** >0.5 = similar fold; >0.7 = high confidence homology.

## Visualization Operators

```python
b.visualize_neighborhood(protein_id, window=12, output_path="output.png")
```

**Multi-source locus diagrams:**
```bash
python scripts/plot_locus_multisource.py \
    --db data/dataset/sharur.duckdb \
    --protein PROTEIN_ID \
    --window 12 \
    --output figures/locus.png
```

Features: Multi-source annotation priority (Foldseek > DefenseFinder > PADLOC > PFAM/KEGG/VOGdb), gene numbers below track, absolute genome coordinates.

## Analysis Manifest System

Each dataset may have a `manifest.json` derived continuity/status cache. Canonical findings
and live database/artifact state remain authoritative. `b.resume()` refreshes the manifest
view in memory from those sources; call `b.manifest.save()` explicitly to persist the
reconciled cache. Legacy manifest shapes are normalized with a warning, and malformed
manifests are not overwritten automatically. Other APIs:
`b.manifest.log_session(phase, note)`. Migration:
`python scripts/migrate_to_manifest.py data/my_dataset/`.

## V2 Semantic Atoms

Full docs: `docs/predicates_v2.md`

```python
b = Sharur("data/YOUR_DATASET/sharur.duckdb", read_only=True)
# New Stage 07 builds already materialize V2. Use generate_v2() for manual
# refreshes or subsets.
b.generate_v2(output_review_queue="review_queue.tsv")

state = b.get_semantic_state("protein_id")
# state.activities, state.roles, state.architecture, state.composite_predicates

b.search_by_facet("activity", atom_ids=["hydrogenase"], relation="implies")
b.search_by_atoms(has=["giant", "unannotated"])
b.search_atoms(atom_id="hydrogenase", relation="implies", source_db="hyddb")
atoms = b.get_atoms("protein_id")  # each has facet, relation, source_db, source_accession
explanation = b.explain("protein_id")  # includes composite_explanations witnesses
b.list_composites()  # YAML-declared rules (config/predicates_v2/composites.yaml)
```

## Hypothesis Tracking & Provenance

```python
b = Sharur("data/YOUR_DATASET/sharur.duckdb", read_only=True)
h = b.propose_hypothesis("Group 4 NiFe hydrogenases are energy-conserving")
b.add_evidence(h.hypothesis_id, "NiFe Group 4 survey", "12/41 genomes", True, 0.8)
print(b.hypothesis_summary())
# Provenance DAG: b.render_provenance(title="Analysis DAG", output_path="figures/provenance.mermaid")
```

Hypotheses persist in `exploration/hypotheses.json` across sessions.
