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

The `environment.yml` installs Bennu in editable mode with all optional extras via `pip install -e ".[all]"`.

### Minimal install (no conda)

If you only need the core Python library without bioinformatics CLI tools:

```bash
pip install -e "."
```

Or with specific extras:

```bash
pip install -e ".[viz,structure]"
```

## Optional Extras

| Extra | Packages | Use case |
|-------|----------|----------|
| `topology` | pyTMHMM | Transmembrane helix prediction |
| `ml` | torch, transformers, h5py, scipy | ESM2 embeddings, ML workflows |
| `viz` | matplotlib, dna-features-viewer, plotly | Locus diagrams, UMAP plots |
| `structure` | biopython, esm | ESM3 structure prediction |
| `agent` | dspy-ai | BennuAgent orchestrator (DSPy) |
| `reports` | reportlab | PDF report generation |
| `ingest` | biopython | GFF/FASTA parsing for ingestion |
| `dev` | pytest, ruff, mypy, etc. | Development and testing |
| `notebooks` | jupyter, ipywidgets, matplotlib | Jupyter notebook workflows |
| `all` | All of the above | Full installation |

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
# Core Python import (should work WITHOUT dspy-ai)
python -c "from bennu.operators import Bennu; print('Core import: OK')"

# pyTMHMM (requires numpy <2.0)
python -c "from bennu.predicates.topology import is_available; print(f'pyTMHMM: {is_available()}')"

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

### `import bennu` fails with dspy-ai error

This should be fixed — `dspy-ai` is no longer a core dependency. If you still see it:

1. Verify you're on the latest code (`git pull`)
2. Reinstall: `pip install -e "."`
3. Check that `bennu/__init__.py` uses `__getattr__` for `BennuAgent`

### DuckDB version conflicts

If you see `duckdb.InvalidInputException` about database versions, your database was created with a different DuckDB version. Either:

- Upgrade DuckDB: `pip install --upgrade duckdb`
- Re-ingest the dataset (see `QUICKSTART.md`)
