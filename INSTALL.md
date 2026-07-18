# Installation Guide

## Prerequisites

- **conda** or **mamba** (recommended: [miniforge](https://github.com/conda-forge/miniforge))
- **git**
- **macOS** (Apple Silicon or Intel) or **Linux** (x86_64)

## Quick Start

```bash
# Clone the repository
git clone https://github.com/jwestrob/Sharur.git
cd Sharur

# Create the conda environment (installs Python, compiled libs, bioinformatics tools)
conda env create -f environment.yml

# Activate the environment
conda activate sharur
```

The `environment.yml` installs Sharur in editable mode with all dependencies plus dev tools via `pip install -e ".[dev]"`. It also installs `macsyfinder` and `mdmparis-defense-finder` (the latter with `--no-deps` to avoid a click version conflict — see Troubleshooting).

After environment creation, download the DefenseFinder model definitions:

```bash
defense-finder update
```

All runtime dependencies (torch, transformers, matplotlib, biopython, plotly, reportlab, jupyter, etc.) are installed by default with a bare `pip install -e "."`. The only optional extra is `[dev]` for testing and linting tools (pytest, ruff, mypy).

### Minimal install (no conda)

If you only need the Python library without bioinformatics CLI tools:

```bash
pip install -e "."        # all runtime deps
pip install -e ".[dev]"   # + pytest, ruff, mypy
```

## External Tools (not in conda)

### Astra annotation pipeline

Astra manages HMM databases for functional annotation.

```bash
# Astra is installed from source (not on conda/PyPI)
cd ~/astra
pip install -e .

# Verify
astra --help
```

**Note:** `--prot_in` expects a **directory** containing `.faa` files, not a single file.

#### Installing Astra HMM databases

```bash
# List available and installed databases
astra initialize --show_available
astra initialize --show_installed

# Install the standard pipeline databases
astra initialize --hmms PFAM
astra initialize --hmms KOFAM
astra initialize --hmms HydDB
astra initialize --hmms DefenseFinder
astra initialize --hmms dbCAN
```

Databases are stored at `~/.config/Astra/`. The standard pipeline (Stage 04) defaults to PFAM, KOFAM, HydDB, DefenseFinder, and dbCAN.

### DefenseFinder co-location models

The co-location validation engine (Stage 07) requires MacSyFinder model definitions in addition to the Astra HMM profiles. If you used `environment.yml`, `mdmparis-defense-finder` is already installed — just download the models:

```bash
defense-finder update
```

This installs XML system definitions to `~/.macsyfinder/models/defense-finder-models/`, which Stage 07 uses to validate raw HMM hits into genuine multi-gene defense systems.

If installing manually (without `environment.yml`):

```bash
pip install macsyfinder==2.1.4
pip install --no-deps mdmparis-defense-finder
defense-finder update
```

The `--no-deps` flag is required because `mdmparis-defense-finder` pins `click==8.0.3`, which conflicts with typer's `click>=8.2.1` requirement. The defense-finder CLI works fine with newer click versions.

### BasicTeX (macOS only, for PDF rendering)

Required for rendering manuscripts with pandoc + xelatex:

```bash
brew install --cask basictex

# Add to PATH (add to your shell profile)
export PATH="/Library/TeX/texbin:$PATH"

# Verify
xelatex --version
```

### ESM API key (for structure prediction)

ESM3 structure prediction requires an API key from EvolutionaryScale:

```bash
# Set in your shell profile
export ESM_API_KEY="your-key-here"
```

## Reference Databases

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

This distinguishes unavailable, stale, and failed dataset/runtime capabilities instead of
reducing the result to a single install check.

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

If you see `duckdb.InvalidInputException` about database versions, your database was created with a different DuckDB version. Either:

- Upgrade DuckDB: `pip install --upgrade duckdb`
- Re-ingest the dataset (see `QUICKSTART.md`)
