# syntax=docker/dockerfile:1
#
# Multi-stage, pixi-based image for Sharur.
#
#   CPU (default):  docker build --target runtime     -t sharur:cpu .
#   GPU:            docker build --target runtime-gpu  -t sharur:gpu .
#
# Requires pixi.lock to exist (generate with `pixi install`). Reference databases
# are NOT baked in — they are large and provisioned into volumes at runtime (see
# docker-compose.yml `db-init`). PDF manuscript rendering (BasicTeX/xelatex, a
# macOS-only path) is intentionally out of scope for this Linux image.
#
# NOTE (partial container): the git-source tools Astra (stage-04 annotation) and
# ELSA (synteny) are NOT installed — their pins are TODO(user). The image runs the
# library + bioconda tools; stage 04 / synteny require the commented steps below.

ARG PIXI_VERSION=0.34.0

# --------------------------------------------------------------------------- #
# Builder (CPU): resolve and install the locked `default` pixi environment.
# --------------------------------------------------------------------------- #
FROM ghcr.io/prefix-dev/pixi:${PIXI_VERSION} AS builder
WORKDIR /app

# Copy only what the environment solve + editable install needs first.
COPY pixi.toml pixi.lock pyproject.toml VERSION ./
COPY sharur/ ./sharur/
COPY src/ ./src/
COPY scripts/ ./scripts/

RUN pixi install --locked --environment default

# Produce an activation script so the conda env (LD_LIBRARY_PATH, PATH, …) is
# fully activated at runtime.
RUN pixi shell-hook --environment default -s bash > /shell-hook.sh \
    && echo 'exec "$@"' >> /shell-hook.sh

# --- TODO(user): pin + install the git-source tools once available ---
# ARG ASTRA_REPO=https://github.com/jwestrob/astra.git
# ARG ASTRA_REF=<TODO(user): commit SHA>
# RUN git clone "$ASTRA_REPO" /opt/astra \
#     && git -C /opt/astra checkout "$ASTRA_REF" \
#     && pixi run -e default pip install --no-deps -e /opt/astra
# ARG ELSA_REPO=<TODO(user): ELSA repo URL — not recorded in the codebase>
# ARG ELSA_REF=<TODO(user): commit SHA>
# RUN git clone "$ELSA_REPO" /opt/elsa \
#     && git -C /opt/elsa checkout "$ELSA_REF" \
#     && pixi run -e default pip install --no-deps -e /opt/elsa

# --------------------------------------------------------------------------- #
# Runtime (CPU): slim image carrying the resolved environment + source.
# --------------------------------------------------------------------------- #
FROM ubuntu:24.04 AS runtime
WORKDIR /app

# Env path must match the builder (/app/.pixi/envs/default) for relocation.
COPY --from=builder /app/.pixi/envs/default /app/.pixi/envs/default
COPY --from=builder /shell-hook.sh /shell-hook.sh
COPY pyproject.toml VERSION ./
COPY sharur/ ./sharur/
COPY src/ ./src/
COPY scripts/ ./scripts/

ENV PATH=/app/.pixi/envs/default/bin:$PATH \
    SHARUR_OPS_DB_PATH=/data/ops/sharur_ops.db

EXPOSE 8811

# Non-strict doctor (always exits 0) — reports health without failing the
# partial image on the absent Astra/ELSA. Switch to `sharur doctor --strict`
# once those tools are baked in.
HEALTHCHECK --interval=30s --timeout=15s --start-period=20s --retries=3 \
    CMD sharur doctor || exit 1

ENTRYPOINT ["/bin/bash", "/shell-hook.sh"]
CMD ["sharur", "--help"]

# --------------------------------------------------------------------------- #
# Builder (GPU): the `gpu` pixi environment (CUDA-enabled PyTorch, linux-64).
# --------------------------------------------------------------------------- #
FROM ghcr.io/prefix-dev/pixi:${PIXI_VERSION} AS builder-gpu
WORKDIR /app
COPY pixi.toml pixi.lock pyproject.toml VERSION ./
COPY sharur/ ./sharur/
COPY src/ ./src/
COPY scripts/ ./scripts/
RUN pixi install --locked --environment gpu
RUN pixi shell-hook --environment gpu -s bash > /shell-hook.sh \
    && echo 'exec "$@"' >> /shell-hook.sh

# --------------------------------------------------------------------------- #
# Runtime (GPU): CUDA runtime base + the resolved `gpu` environment.
# --------------------------------------------------------------------------- #
FROM nvidia/cuda:12.4.1-runtime-ubuntu24.04 AS runtime-gpu
WORKDIR /app
COPY --from=builder-gpu /app/.pixi/envs/gpu /app/.pixi/envs/gpu
COPY --from=builder-gpu /shell-hook.sh /shell-hook.sh
COPY pyproject.toml VERSION ./
COPY sharur/ ./sharur/
COPY src/ ./src/
COPY scripts/ ./scripts/
ENV PATH=/app/.pixi/envs/gpu/bin:$PATH \
    SHARUR_OPS_DB_PATH=/data/ops/sharur_ops.db
EXPOSE 8811
HEALTHCHECK --interval=30s --timeout=15s --start-period=20s --retries=3 \
    CMD sharur doctor || exit 1
ENTRYPOINT ["/bin/bash", "/shell-hook.sh"]
CMD ["sharur", "--help"]
