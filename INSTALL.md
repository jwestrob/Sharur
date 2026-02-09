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

The `environment.yml` installs Sharur in editable mode with all dependencies plus dev tools via `pip install -e ".[dev]"`.

All runtime dependencies (torch, transformers, matplotlib, biopython, plotly, reportlab, jupyter, etc.) are installed by default with a bare `pip install -e "."`. The only optional extra is `[dev]` for testing and linting tools (pytest, ruff, mypy).

### Minimal install (no conda)

If you only need the Python library without bioinformatics CLI tools:

```bash
pip install -e "."        # all runtime deps
pip install -e ".[dev]"   # + pytest, ruff, mypy
```

## External Tools (not in conda)

### Astra annotation pipeline

Astra manages pre-installed HMM databases (PFAM, KOFAM, VOGdb, HydDB, DefenseFinder, etc.).

```bash
# Astra is installed from source (not on conda/PyPI)
cd ~/astra
pip install -e .

# Verify
astra --help

# Example: run DefenseFinder against a protein directory
astra search --installed_hmms DefenseFinder --threads 12 \
    --prot_in source/ --outdir annotations/ --cut_ga
```

**Note:** `--prot_in` expects a **directory** containing `.faa` files, not a single file.

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

Used by `scripts/classify_hydrogenases.py` for hydrogenase classification:

```
data/reference/hyddb/HydDB_all.dmnd
```

### Foldseek databases

Local Foldseek databases are stored at `~/.foldseek/`:

```bash
# Download databases (pdb100, afdb50, afdb-swissprot)
foldseek databases PDB100 ~/.foldseek/pdb100/pdb100 /tmp/foldseek
foldseek databases AlphaFoldDB50 ~/.foldseek/afdb50/afdb50 /tmp/foldseek
```

### PFAM / KOFAM / other HMM databases

Managed by Astra. After installing Astra, databases are stored at `~/.config/Astra/`.

## Verification

Run these commands to confirm everything is installed correctly:

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

pyTMHMM is compiled against NumPy 1.x and crashes with NumPy 2.x. The `environment.yml` pins `numpy<2.0` to prevent this. If you installed without conda:

```bash
pip install "numpy>=1.24.0,<2.0"
pip install --no-build-isolation pyTMHMM
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
