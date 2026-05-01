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
- Output: `<outdir>/<database>_results/tmp_results/bulk_results.tsv`
- For single files: `mkdir source/ && cp proteins.faa source/`

### Secretion System Identification (TXSScan)
**Validation:** Use the co-location engine (see section above), NOT `scripts/validate_secretion_systems.py`.
**Requires:** TXSScan HMMs via Astra (`-d TXSScan` in stage 04), TXSScan models in `~/.macsyfinder/models/`
**Output:** `secretion_systems` table + `txsscan_system` annotations in DuckDB.
**Systems detected:** T1SS–T9SS, T4P, Tad, Flagellum, MSH, ComM (280 HMM profiles).

## Defense System Validation (Co-location Engine)

**Module:** `sharur/colocation.py`
**Purpose:** Validates raw HMM hits (DefenseFinder/TXSScan) into genuine multi-gene systems using MacSyFinder's XML model definitions, but runs in-process with DuckDB — **~1-2 seconds** vs 20-30 minutes for the MacSyFinder subprocess.

**ALWAYS use this instead of `scripts/validate_defense_systems.py`** (which calls MacSyFinder single-threaded). The old script was written for the susan_genomes dataset before this engine existed.

### Usage

```python
import logging
logging.basicConfig(level=logging.INFO, format='%(message)s')

from sharur.colocation import validate_systems, integrate_defense_results, integrate_secretion_results

db_path = "data/DATASET/sharur.duckdb"

# DefenseFinder
systems_df, genes_df = validate_systems(db_path, source="defensefinder")
integrate_defense_results(db_path, systems_df, genes_df)

# TXSScan (if loaded)
systems_df, genes_df = validate_systems(db_path, source="txsscan")
integrate_secretion_results(db_path, systems_df, genes_df)
```

### What it does
1. Parses 556 DefenseFinder XML model definitions from `~/.macsyfinder/models/`
2. Loads all HMM hits from DuckDB `annotations` table (single indexed query)
3. Clusters hits on contigs within each model's `inter_gene_max_space`
4. Validates clusters against quorum rules (mandatory/accessory/forbidden genes)
5. Resolves conflicts (greedy score-ranked non-overlapping selection)
6. Writes validated systems to `defense_systems` table + `defensefinder_system` annotations

### Key numbers (omni_production benchmark)
- 312,144 HMM hits → 18,430 systems (42,954 genes) in 1.3 seconds
- 110,210 raw proteins → 22,158 system-validated proteins (79.9% FP reduction)
- Typical FP rates: 80-85% of raw HMM hits are not part of genuine multi-gene systems

### Requirements
- MacSyFinder model XMLs at `~/.macsyfinder/models/defense-finder-models/` (installed by defense-finder)
- Raw HMM hits loaded in DuckDB `annotations` table with `source='defensefinder'`
- Astra `--write_macsyfinder` is NOT needed (engine reads from DuckDB, not hmmsearch files)

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
**Pipeline:** Integrated into stage 07 before final V2 predicate generation.

## ESM2 Embeddings
**Script:** `src/ingest/06_esm2_embeddings.py`
**Usage:** `python src/ingest/06_esm2_embeddings.py data/DATASET/stage03_prodigal data/DATASET/embeddings/`
**Model:** `facebook/esm2_t6_8M_UR50D` (320-dim, auto-detects MPS/CUDA/CPU)
**Output:** `embeddings/protein_embeddings.h5` + `embedding_manifest.json` (FAISS index built at runtime from H5)
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
b = Sharur("data/DATASET/sharur.duckdb")

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

Each dataset has a `manifest.json` for session continuity. Key APIs: `b.resume()` (status overview), `b.manifest.log_session(phase, note)`, `b.manifest.save()`. Migration: `python scripts/migrate_to_manifest.py data/my_dataset/`

## V2 Semantic Atoms

Full docs: `docs/predicates_v2.md`

```python
b = Sharur("data/YOUR_DATASET/sharur.duckdb")
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
b = Sharur("data/YOUR_DATASET/sharur.duckdb")
h = b.propose_hypothesis("Group 4 NiFe hydrogenases are energy-conserving")
b.add_evidence(h.hypothesis_id, "NiFe Group 4 survey", "12/41 genomes", True, 0.8)
print(b.hypothesis_summary())
# Provenance DAG: b.render_provenance(title="Analysis DAG", output_path="figures/provenance.mermaid")
```

Hypotheses persist in `exploration/hypotheses.json` across sessions.
