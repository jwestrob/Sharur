# Deploying Sharur

Sharur can be consumed three ways, from lightest to most turnkey. Whichever you
choose, **run `sharur doctor` first** — it reports which external tools, reference
databases, and API keys are available.

> **Partial-container / partial-environment note.** Two tools are installed from
> git checkouts and are currently left as `TODO(user)` pins: **Astra** (stage-04
> HMM annotation) and **ELSA** (synteny). Until they are pinned (`pixi.toml`,
> `Dockerfile`, `build-container.yml`) they must be installed manually per
> `INSTALL.md`. Everything else — the library, the bioconda tools, embeddings,
> structure search — works without them.

## 1. Library (programmatic / notebook use)

For importing `sharur` in Python without the bioinformatics binaries:

```bash
pip install -e .            # or: pixi install
python -c "from sharur.operators import Sharur; print('ok')"
```

This is enough for querying an existing `sharur.duckdb`, predicate search, and the
Python API. It does **not** provide prodigal/diamond/foldseek/etc.

## 2. Full analysis environment (pixi)

The reproducible way to get the full toolchain locally:

```bash
pixi install               # locked conda (bioconda) + PyPI environment
pixi run doctor            # verify tools / reference DBs
```

Then provision reference databases out-of-band (see `INSTALL.md`):

```bash
astra initialize --hmms PFAM KOFAM HydDB DefenseFinder dbCAN
defense-finder update
foldseek databases Alphafold/UniProt50 ~/.foldseek/afdb50 /tmp/foldseek-tmp
```

Set `ESM_API_KEY` for ESM3 structure prediction (optional).

## 3. Container (homelab: Docker Compose)

The image bundles the pixi environment (library + bioconda tools). Reference
databases are **not** baked in — they live in named volumes populated once.

```bash
docker compose build

# One-time reference DBs (profile-gated so it never runs implicitly).
# FOLDSEEK_DBS picks which structure DBs to fetch (comma-separated).
FOLDSEEK_DBS="Alphafold/UniProt50" docker compose --profile init run --rm db-init

# Coordination server on :8811 (SQLite persisted to the ops-db volume).
docker compose up sharur-ops

# One-shot CLI runs (datasets bind-mounted from ./data).
docker compose run --rm sharur-cli overview --db data/<dataset>/sharur.duckdb
docker compose run --rm sharur-cli doctor
```

### GPU image

```bash
docker build --target runtime-gpu -t sharur:gpu .
docker run --gpus all --rm sharur:gpu sharur doctor
```

### Environment variables

| Variable | Purpose | Default |
|----------|---------|---------|
| `ESM_API_KEY` | ESM3 structure-prediction API (EvolutionaryScale) | unset (disabled) |
| `SHARUR_OPS_DB_PATH` | SQLite path for the ops coordination server | `/data/ops/sharur_ops.db` (image) |
| `FOLDSEEK_DBS` | Which Foldseek DBs the `db-init` step downloads | `Alphafold/UniProt50` |

## 4. HPC (Apptainer / Singularity)

Any OCI image Sharur publishes is directly runnable on shared clusters where you
can't run a Docker daemon — Apptainer converts it to a SIF on pull:

```bash
apptainer pull sharur.sif docker://<registry>/sharur:<tag>
apptainer exec sharur.sif sharur doctor          # smoke check
apptainer exec --nv sharur.sif sharur ...        # --nv for GPU nodes
```

Bind your data and reference-DB directories with `--bind`, e.g.
`apptainer exec --bind /scratch/data:/app/data sharur.sif sharur overview --db data/<dataset>/sharur.duckdb`.

> The publish registry is a TODO — GHCR (`ghcr.io/jwestrob/sharur`) is the starting
> point, kept configurable in `build-container.yml`.

## Building & publishing images (CI)

Image builds run only via the **Build and Push Container** GitHub Action
(`workflow_dispatch`), which takes `registry`, `tag`, and `push` inputs (push
defaults to off). See `.github/workflows/build-container.yml`. All Sharur CI is
manually triggered — nothing runs automatically on push/PR yet.

## Future distribution channels (not yet enabled)

- **PyPI** (`pip install sharur`) — enabled by the release scaffolding once opted
  into; covers the Python library only (binaries still need conda/container).
- **Bioconda** — the community standard, but blocked until Astra & ELSA are
  packaged upstream (bioconda accepts only conda-forge/bioconda deps) and Sharur
  is on PyPI with tagged releases.
