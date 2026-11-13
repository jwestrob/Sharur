# Installation Guide

## Prerequisites

- **conda** or **mamba** (recommended: [miniforge](https://github.com/conda-forge/miniforge))
- **git** and **curl**
- **Linux x86-64** for the full ingest pipeline. [Aksha](https://github.com/jwestrob/aksha), the
  Stage 04 HMM search, ships its compiled runtime for Linux x86-64 with Python 3.12. On macOS
  (Apple Silicon or Intel) the same environment installs everything else, which covers querying,
  browsing and analyzing a dataset built on Linux.
- **Disk:** about 14 GB for the environment (19 GB while package caches fill) and 23 GB for the
  Stage 04 HMM databases (KOFAM 15 GB, PFAM 7.3 GB). Aksha stores
  databases under `~/.config/Astra`; on a machine with a small home quota, point Aksha's `db_path`
  elsewhere before installing them (see Aksha's
  [database setup guide](https://github.com/jwestrob/aksha/blob/main/docs/database-setup.md)).
- **GPU:** optional, for Stage 06 embeddings; see
  [Without a GPU](QUICKSTART.md#without-a-gpu).

## Quick Start

```bash
# Clone the repository
git clone https://github.com/jwestrob/Sharur.git
cd Sharur

# Create and activate the conda environment: Python, compiled libraries, bioinformatics tools,
# Sharur (editable, with every optional extra), Aksha (Linux x86-64) and MacSyFinder
conda env create -f environment.yml
conda activate sharur

# Stage 04 HMM databases (about 2.5 GB to download; KOFAM is the largest)
for db in PFAM KOFAM HydDB DefenseFinder dbCAN; do aksha initialize --hmms "$db"; done

# DefenseFinder system models for Stage 07's defense-system caller (~/.macsyfinder/models)
macsydata install --user --org mdmparis defense-finder-models

# HydDB reference for Stage 07's hydrogenase subgroup calls (data/reference/hyddb)
mkdir -p data/reference/hyddb
curl -L -o data/reference/hyddb/HydDB_all_hydrogenases.faa \
  https://raw.githubusercontent.com/GreeningLab/HydDB/main/fastas/HydDB_all_hydrogenases.faa
diamond makedb --in data/reference/hyddb/HydDB_all_hydrogenases.faa --db data/reference/hyddb/HydDB_all

# KEGG predicate map, built on this machine from KEGG REST (academic use; see below)
sharur setup-kegg

# Every core requirement present?
sharur doctor --strict
```

`environment.yml` installs Sharur in editable mode with all optional Python surfaces plus
development tools (`pip install -e ".[all,dev]"`), Aksha from PyPI on Linux x86-64, and
MacSyFinder, whose `macsydata` command fetches the DefenseFinder models. The DefenseFinder CLI
itself is optional (see [DefenseFinder co-location models](#defensefinder-co-location-models)).

The base install contains the DuckDB query/predicate surface and ingest orchestration without
model frameworks, vector libraries, web serving, plotting, or notebooks. Add only the
capabilities needed:

| Extra | Adds |
|---|---|
| `parquet` | Optional PyArrow interoperability |
| `vectors` | H5 inspection and persistent FAISS similarity |
| `embeddings` | `vectors` dependencies plus Torch/Transformers embedding inference |
| `ops` | FastAPI/uvicorn/Jinja2: multi-agent coordination service, query service and the `sharur browse` website |
| `visualization` | Matplotlib, Plotly, and neighborhood rendering |
| `structure` | Biopython and ESM structure prediction |
| `reports` | PDF report backends |
| `notebooks` | Jupyter and ipywidgets |
| `all` | Every optional runtime surface |
| `docs` | MkDocs Material for building the documentation site |
| `dev` | Tests, HTTP test client, linting, typing, and wheel building |

### Minimal install (no conda)

If you only need the Python library without bioinformatics CLI tools:

```bash
pip install -e "."                              # lean core
pip install -e ".[embeddings]"                  # standard ingest including Stage 06
pip install -e ".[ops]"                         # coordination server, query service, dataset browser
pip install -e ".[dev,vectors,ops,reports]"     # regression-test surfaces
pip install -e ".[all,dev]"                     # full development environment
```

## External Tools

### Aksha HMM search

[Aksha](https://github.com/jwestrob/aksha) runs Stage 04's HMM searches and manages the HMM
databases. It is the successor to Astra and accepts the same `search` flags; Stage 04 falls
back to an `astra` executable when `aksha` is absent. `environment.yml` installs it from PyPI on
Linux x86-64; in another Python 3.12 environment:

```bash
python -m pip install aksha
# or from source:
git clone https://github.com/jwestrob/aksha && pip install -e aksha

aksha --help
```

Aksha's compiled runtime (`aksha-runtime`) ships for Linux x86-64 with Python 3.12; see Aksha's
`release/INSTALL.md` for other platforms. `sharur doctor` reports whether the runtime imports and
which of the five standard databases are installed.

**Note:** `--prot_in` expects a **directory** containing `.faa` files, not a single file.

#### Installing HMM databases

```bash
# List available and installed databases
aksha initialize --show_available
aksha initialize --show_installed

# Install the standard pipeline databases
aksha initialize --hmms PFAM
aksha initialize --hmms KOFAM
aksha initialize --hmms HydDB
aksha initialize --hmms DefenseFinder
aksha initialize --hmms dbCAN
```

Database files are stored at `~/.config/Astra/` (the `db_path` Aksha inherits from Astra); the
five standard databases take about 23 GB there (KOFAM 15 GB, PFAM 7.3 GB). Installing them
downloads about 2.4 GB and took 1.5 h in a test, most of it KOFAM's 1.6 GB from genome.jp.
Aksha's registry of installed databases lives in its own config directory
(`~/.config/Aksha/hmm_databases.json` on Linux). The standard pipeline (Stage 04) defaults to
PFAM, KOFAM, HydDB, DefenseFinder, and dbCAN.

**Moving from Astra.** Databases installed through Astra (or KOFAM installed by Aksha 0.2.1) are
already on disk; register them with Aksha instead of downloading them again:

```bash
aksha initialize --show_available     # creates Aksha's registry on first run
sharur adopt-astra-hmms --dry-run     # list what will be registered
sharur adopt-astra-hmms               # mark them installed (backup: hmm_databases.json.bak)
```

### Protein families and remote homology (optional Stages 09-11)

`environment.yml` installs MMseqs2 for Stage 09 (`--with-families`) and, through dbCAN, BLAST+
(`rpsblast`) for Stage 10. Stage 11 also needs HH-suite and MAFFT:

```bash
conda install -c conda-forge -c bioconda hhsuite mafft
```

Both stages read databases you name in a configuration file (CDD for Stage 10; Pfam and PDB in
HH-suite format for Stage 11); [Remote homology](docs/guides/remote-homology.md) lists them, and
`sharur remote-homology check --config CONFIG` reports what is missing.

### DefenseFinder co-location models

The co-location validation engine (Stage 07) reads MacSyFinder model definitions in addition to
the DefenseFinder and TXSScan HMM profiles. MacSyFinder (installed by `environment.yml`, or
`pip install macsyfinder==2.1.4`) downloads them:

```bash
macsydata install --user --org mdmparis defense-finder-models
```

This installs the XML system definitions to `~/.macsyfinder/models/defense-finder-models/`,
which Stage 07 uses to validate raw HMM hits into multi-gene defense systems.

The DefenseFinder CLI (`defense-finder update` installs the same models) is optional; Sharur
uses it only in `scripts/run_defensefinder_systems.py`. Install it after the environment, with
`--no-deps` because `mdmparis-defense-finder` pins `click==8.0.3` and typer needs
`click>=8.2.1` (the CLI works with the newer click):

```bash
pip install --no-deps mdmparis-defense-finder
```

### BasicTeX (macOS only, for PDF rendering)

Required for rendering manuscripts with pandoc + xelatex:

```bash
brew install --cask basictex

# Add to PATH (add to your shell profile)
export PATH="/Library/TeX/texbin:$PATH"

# Verify
xelatex --version
```

### ESM API key (for quick-look structure prediction)

The `predict_structure` operators call the ESM3 API, which needs an API key from
EvolutionaryScale. Protenix v2 and ESMFold2 run locally on a GPU node and need no key:

```bash
# Set in your shell profile
export ESM_API_KEY="your-key-here"
```

## Reference Databases

### KEGG predicate map (built locally)

KEGG data is subject to KEGG's terms, so Sharur ships rules and builds the KO → predicate map on
your machine:

```bash
sharur setup-kegg                        # KEGG REST (academic use), about 3 minutes
sharur setup-kegg --inputs /path/to/kegg # build from a licensed KEGG copy
```

The map lands in `$SHARUR_KEGG_DIR` when set, else in `data/reference/kegg` of a source checkout,
else in `~/.sharur/kegg`; Sharur reads it from the same places. See [`DATA_LICENSES.md`](DATA_LICENSES.md).

### VOGdb annotation table

VOG descriptions and functional categories come from VOGdb's `vog.annotations.tsv`:

```bash
mkdir -p data/reference/vogdb
curl -L -o data/reference/vogdb/vog.annotations.tsv.gz \
  https://fileshare.lisc.univie.ac.at/vog/latest/vog.annotations.tsv.gz
```

Sharur looks for `vog.annotations.tsv` (or `.gz`) in `$SHARUR_VOG_ANNOTATIONS`, `data/reference/vogdb/`,
`~/.sharur/vogdb/` and `~/.config/Astra/VOGdb/`.

### Annotation name maps

Stage 07 enriches annotation rows from `pfam_id_desc.tsv` and KOfam's `ko_list` when present. It
looks in `$SHARUR_REFERENCE_DIR`, the `reference` directory beside the dataset directory (for
`data/my_dataset`, `data/reference`), `./data/reference` and `~/.sharur/reference`; the KOFAM install's
`~/.config/Astra/ko_list` also serves. Without `pfam_id_desc.tsv`, Stage 07 reads Pfam accessions,
names and descriptions from the `Pfam-A.hmm` that `aksha initialize --hmms PFAM` installed, the
file Stage 04 searched.

### CoverM (optional)

Per-sample coverage (`sharur-ingest --reads`, Stage 08) uses [CoverM](https://github.com/wwood/CoverM):

```bash
conda install -c bioconda coverm
```

### HydDB DIAMOND database

Used by `scripts/classify_hydrogenases.py` for hydrogenase subgroup classification:

```bash
mkdir -p data/reference/hyddb
curl -L -o data/reference/hyddb/HydDB_all_hydrogenases.faa \
  https://raw.githubusercontent.com/GreeningLab/HydDB/main/fastas/HydDB_all_hydrogenases.faa
diamond makedb --in data/reference/hyddb/HydDB_all_hydrogenases.faa \
  --db data/reference/hyddb/HydDB_all
```

This produces `data/reference/hyddb/HydDB_all.dmnd`, which Stage 07 uses automatically when HydDB annotations are present.

### Foldseek databases

Local Foldseek databases are stored at `~/.foldseek/`:

```bash
# Download databases (pdb100, afdb50, afdb-swissprot)
foldseek databases PDB100 ~/.foldseek/pdb100/pdb100 /tmp/foldseek
foldseek databases AlphaFoldDB50 ~/.foldseek/afdb50/afdb50 /tmp/foldseek
```

## Verification

The quickest check is **`sharur doctor`**, which checks the installed ingest entry point,
external tools, reference directories, and API keys and flags whatever is missing:

```bash
sharur doctor            # one-shot health check (tools, reference DBs, API keys)
sharur doctor --strict   # exit non-zero if a core requirement is missing (for scripts)
```

Once a dataset exists, `doctor` is complemented by the non-mutating, dataset-aware brief:

```bash
sharur preflight --db data/DATASET/sharur.duckdb
sharur preflight --db data/DATASET/sharur.duckdb --format json --strict
```

This distinguishes unavailable, stale, and failed dataset/runtime capabilities, giving each its
own status.

To look around a dataset in a web browser:

```bash
sharur browse --db data/DATASET/sharur.duckdb   # http://localhost:8800/
```

For granular manual checks:

```bash
# Core Python import
python -c "from sharur.operators import Sharur; print('Core import: OK')"

# pyTMHMM (requires numpy <2.0)
python -c "from sharur.predicates.topology import is_available; print(f'pyTMHMM: {is_available()}')"

# Bioinformatics CLI tools
prodigal -v
diamond version
hmmsearch -h | head -2
foldseek version

# Run tests
python -m pytest tests/ --override-ini addopts="" -x -q
```

## Troubleshooting

### pyTMHMM: "numpy.core" import errors

pyTMHMM is compiled against NumPy 1.x and crashes with NumPy 2.x. It is optional — topology predicates (transmembrane, signal peptide) are silently skipped when unavailable. If you need topology predictions:

```bash
pip install "numpy>=1.24.0,<2.0"
pip install --no-build-isolation pyTMHMM
```

Note: this downgrades numpy, which may affect other packages.

### click version conflict (defense-finder vs typer)

`mdmparis-defense-finder` pins `click==8.0.3` but `typer` requires `click>=8.2.1`. The `environment.yml` resolves this by installing defense-finder with `--no-deps`. If you see `TypeError: type 'Choice' is not subscriptable`, your click was downgraded:

```bash
pip install "click>=8.2.1"
```

### torch: MPS on Apple Silicon

PyTorch supports Apple Silicon GPUs via the MPS backend. If you see MPS-related warnings, you can disable it:

```bash
export PYTORCH_MPS_HIGH_WATERMARK_RATIO=0.0
```

Or force CPU-only in Python:

```python
import torch
device = torch.device("cpu")
```

### DuckDB version conflicts

Sharur runs on DuckDB 1.x (`pyproject.toml` pins `duckdb<2`). DuckDB 2.0 stops Stage 07 with
`ResourceInUseException: connection has an unfinished result from execute()` and leaves the SQL
surface without its macros; `sharur doctor` names the installed version. Install a supported
release with `pip install -U "duckdb<2"`.

If you see `duckdb.InvalidInputException` about database versions, your database was created with
a different DuckDB version. Either:

- Install the DuckDB 1.x release that wrote it: `pip install -U "duckdb<2"`
- Re-ingest the dataset (see `QUICKSTART.md`)
