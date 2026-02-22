# =============================================================================
#  CellJanus — Docker image
#  Dual-Perspective Host–Microbe Deconvolution Pipeline
#
#  Build:
#    docker build -t celljanus:0.1.1 .
#
#  Run examples:
#    # Check tool availability
#    docker run --rm celljanus:0.1.1 celljanus check
#
#    # Run full pipeline (mount data directory)
#    docker run --rm -v /path/to/data:/data celljanus:0.1.1 \
#        celljanus run \
#            --read1 /data/sample_R1.fastq.gz \
#            --read2 /data/sample_R2.fastq.gz \
#            --host-index /data/refs/bowtie2_index/GRCh38_noalt_as \
#            --kraken2-db /data/refs/standard_8 \
#            --output-dir /data/results \
#            --threads 8
# =============================================================================

# ---------- Stage 1: Miniforge + bioconda tools ----------
FROM condaforge/miniforge3:24.11.3-2 AS base

LABEL maintainer="Zhaoqing Wang" \
    description="CellJanus: Dual-Perspective Host–Microbe Deconvolution" \
    version="0.1.1" \
    url="https://github.com/zhaoqing-wang/CellJanus"

# Avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# System dependencies (minimal — most tools come from conda)
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    procps \
    && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Create the CellJanus conda environment with pinned versions
RUN conda create -y -n celljanus \
    -c bioconda -c conda-forge \
    python=3.11 \
    fastp=1.1.0 \
    bowtie2=2.5.4 \
    samtools=1.23 \
    kraken2=2.17.1 \
    bracken=3.1 \
    && conda clean -afy

# Activate the environment by default
ENV PATH="/opt/conda/envs/celljanus/bin:${PATH}" \
    CONDA_DEFAULT_ENV=celljanus

# ---------- Stage 2: Install CellJanus Python package ----------
WORKDIR /opt/celljanus

# Copy only dependency specification first (Docker layer caching)
COPY pyproject.toml README.md LICENSE ./

# Copy the package source
COPY celljanus/ ./celljanus/
COPY tests/ ./tests/
COPY testdata/ ./testdata/

# Install CellJanus into the conda env
RUN pip install --no-cache-dir . && \
    pip install --no-cache-dir ".[dev]"

# Smoke test — verify all tools are importable and discoverable
RUN python -c "import celljanus; print(f'CellJanus {celljanus.__version__} installed')" && \
    celljanus check && \
    python -m pytest tests/ -v --tb=short && \
    celljanus qc --read1 testdata/reads_SE.fastq.gz --output-dir /tmp/smoke_qc && \
    echo 'Smoke tests passed' && rm -rf /tmp/smoke_qc

# ---------- Runtime ----------
WORKDIR /data

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "celljanus"]
CMD ["celljanus", "--help"]
