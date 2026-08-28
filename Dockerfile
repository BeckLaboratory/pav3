# syntax=docker/dockerfile:1
###############################################################################
# PAV3 container image
#
# The image behaves exactly like the `pav3` command (it is the ENTRYPOINT).
#
# Build (run from the repo root; a built wheel must already exist in dist/):
#   docker build -t pav3:latest .
#
# Run (mount the working directory and run as the calling user):
#   docker run --rm -v "$PWD:$PWD" -w "$PWD" --user "$(id -u):$(id -g)" \
#       pav3:latest batch
#
# Convert to Singularity/Apptainer:
#   apptainer build pav3.sif docker-daemon://pav3:latest
#
# Published images (no login required):
#   docker pull ghcr.io/becklaboratory/pav3:latest
#   apptainer pull oras://ghcr.io/becklaboratory/pav3-sif:latest
#
# The published .sif is signed; check it with "apptainer verify pav3_<ver>.sif".
###############################################################################

# --- Pinned versions ---------------------------------------------------------
# Change a version here and every dependency is re-resolved at build time.
ARG PYTHON_VERSION=3.12
ARG MINIMAP2_VERSION=2.31
# UCSC tools are versioned by dated directory on hgdownload, e.g.
#   https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v472/
# Browse https://hgdownload.soe.ucsc.edu/admin/exe/ for available versions.
ARG UCSC_VERSION=479


###############################################################################
# Stage 1 - fetch pinned external binaries.
#
# Both tools ship as self-contained, statically-built x86_64 binaries, so there
# is no compiler or htslib to manage. wget/ca-certificates stay in this
# throwaway stage and never reach the runtime image.
###############################################################################
FROM python:${PYTHON_VERSION}-slim AS pav3_tools

ARG MINIMAP2_VERSION
ARG UCSC_VERSION

RUN set -eux; \
    apt-get update; \
    apt-get install -y --no-install-recommends ca-certificates wget bzip2; \
    rm -rf /var/lib/apt/lists/*

WORKDIR /tools/bin

# minimap2: official prebuilt release (binary is self-contained).
RUN set -eux; \
    base="minimap2-${MINIMAP2_VERSION}_x64-linux"; \
    wget -q "https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VERSION}/${base}.tar.bz2"; \
    tar -xjf "${base}.tar.bz2"; \
    install -m 0755 "${base}/minimap2" ./minimap2; \
    rm -rf minimap2-*; \
    ./minimap2 --version

# UCSC bedToBigBed: pinned versioned directory.
RUN set -eux; \
    url="https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v${UCSC_VERSION}/bedToBigBed"; \
    wget -q "$url" -O bedToBigBed; \
    chmod 0755 bedToBigBed


###############################################################################
# Stage 2 - runtime image.
###############################################################################
FROM python:${PYTHON_VERSION}-slim AS pav3

# Set by files/docker/build (the recommended entry point). When empty, the
# newest wheel in dist/ is installed instead.
ARG PAV_VERSION

LABEL org.opencontainers.image.title="PAV3" \
      org.opencontainers.image.description="Assembly-based variant caller for phased genomes" \
      org.opencontainers.image.source="https://github.com/BeckLaboratory/pav3" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.version="${PAV_VERSION}"

# External binaries from the tools stage.
COPY --from=pav3_tools /tools/bin/ /opt/pav/bin/
ENV PATH="/opt/pav/bin:${PATH}"

# Install the locally-built PAV3 wheel. pip resolves every Python dependency
# (including agglovar) from PyPI. When PAV_VERSION is given that exact wheel is
# installed; otherwise the newest wheel in dist/ is selected by version sort.
#
# The bind mount reads dist/ straight from the build context, so the wheels are
# never baked into an image layer.
RUN --mount=type=bind,source=dist,target=/dist \
    set -eux; \
    if [ -n "${PAV_VERSION}" ]; then \
        wheel="/dist/pav3-${PAV_VERSION}-py3-none-any.whl"; \
    else \
        wheel="$(ls /dist/pav3-*.whl | sort -V | tail -n1)"; \
    fi; \
    echo "Installing ${wheel}"; \
    pip install --no-cache-dir "$wheel"; \
    pav3 --version

# Container-friendly writable cache/config locations. Lets the image run as an
# arbitrary UID (Singularity, `--user`) and on a read-only root filesystem.
ENV HOME=/tmp \
    XDG_CACHE_HOME=/tmp/.cache \
    MPLCONFIGDIR=/tmp/.cache/matplotlib

ENTRYPOINT ["pav3"]
CMD ["--help"]
