# CellJanus

**A Dual-Perspective Tool for Deconvolving Host Single-Cell and Microbial Transcriptomes**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.9+](https://img.shields.io/badge/Python-3.9%2B-blue.svg)](https://www.python.org/)
[![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg)](https://hub.docker.com/)

CellJanus separates host gene expression from microbial abundance in the same
FASTQ data. It orchestrates proven C/C++ tools (fastp, Bowtie2, samtools,
Kraken2, Bracken) through a single Python CLI.

```
FASTQ ──▶ fastp (QC) ──▶ Bowtie2 (align hg38) ──▶ samtools (extract unmapped)
                                                         │
                                 ┌───────────────────────┘
                                 ▼
                          Kraken2 + Bracken ──▶ Abundance tables + Plots
```

---

## Table of Contents

- [Installation (Docker)](#installation)
- [Try It Now (Test Data)](#try-it-now)
- [Quick Start (Real Data)](#quick-start)
- [Alternative Installation](#alternative-installation)
- [Commands](#commands)
- [Project Structure](#project-structure)
- [Citation](#citation)

---

## Installation

Docker is the recommended installation method — it bundles Python, all
bioinformatics tools, and test data in a single image.

### Prerequisites

- [Docker](https://docs.docker.com/get-docker/) (Docker Desktop on Windows/macOS, or `docker.io` on Linux)

### Build the Image

```bash
git clone https://github.com/zhaoqing-wang/CellJanus.git
cd CellJanus
docker build -t celljanus:0.1.1 .
```

The build takes ~5–10 min on a fresh machine. It automatically:
1. Installs Python 3.11 + all bioinformatics tools via bioconda
2. Installs the CellJanus Python package
3. Runs the full test suite (18 tests)
4. Runs a smoke QC on the bundled test data

### Verify

```bash
docker run --rm celljanus:0.1.1 celljanus check
```

Expected output — all core tools show **✔ Found**:

```
┃ Tool          ┃ Status    ┃ Path                                        ┃
│ fastp         │ ✔ Found   │ /opt/conda/envs/celljanus/bin/fastp         │
│ bowtie2       │ ✔ Found   │ /opt/conda/envs/celljanus/bin/bowtie2       │
│ samtools      │ ✔ Found   │ /opt/conda/envs/celljanus/bin/samtools      │
│ kraken2       │ ✔ Found   │ /opt/conda/envs/celljanus/bin/kraken2       │
│ bracken       │ ✔ Found   │ /opt/conda/envs/celljanus/bin/bracken       │
```

> **Note:** STAR appears as "Missing" — this is expected. STAR is an optional
> roadmap feature and is not required for the core pipeline.

---

## Try It Now

The repository ships real-format FASTQ test data (`testdata/`) so you can
test the pipeline immediately — **no downloads needed**.

### Test Data Overview

| File | Reads | Description |
|------|------:|-------------|
| `reads_R1.fastq.gz` | 1,000 | Paired-end R1 |
| `reads_R2.fastq.gz` | 1,000 | Paired-end R2 |
| `reads_SE.fastq.gz` | 1,000 | Single-end |

Read composition: 600 human (HBB, ACTB), 300 microbial (*S. aureus*, *E. coli*,
*K. pneumoniae* 16S rRNA), 100 low-quality reads with adapter read-through.

### Run QC on Test Data (Docker)

```bash
# Paired-end QC
docker run --rm -v $(pwd)/results:/data celljanus:0.1.1 \
    celljanus qc \
        --read1 /opt/celljanus/testdata/reads_R1.fastq.gz \
        --read2 /opt/celljanus/testdata/reads_R2.fastq.gz \
        --output-dir /data/qc_output

# Single-end QC
docker run --rm -v $(pwd)/results:/data celljanus:0.1.1 \
    celljanus qc \
        --read1 /opt/celljanus/testdata/reads_SE.fastq.gz \
        --output-dir /data/qc_se_output
```

Expected result: **90% reads retained** (100 low-quality pairs filtered),
Q20 improves from 88% → 98%.

### Run Tests (Docker)

```bash
docker run --rm celljanus:0.1.1 \
    python -m pytest /opt/celljanus/tests/ -v
```

Expected: **18 passed**.

---

## Quick Start

For real data, you need reference databases (hg38 genome + Kraken2 DB).

### 1. Download Reference Data

```bash
mkdir -p refs

# Human genome hg38 + pre-built Bowtie2 index (~5 GB)
docker run --rm -v $(pwd)/refs:/data celljanus:0.1.1 \
    celljanus download hg38 --output-dir /data

# Kraken2 database (~8 GB, fits in 32 GB RAM)
docker run --rm -v $(pwd)/refs:/data celljanus:0.1.1 \
    celljanus download kraken2 --output-dir /data --db-name standard_8
```

### 2. Run Full Pipeline

```bash
docker run --rm \
    -v /path/to/your/data:/data \
    -v $(pwd)/refs:/refs \
    celljanus:0.1.1 \
    celljanus run \
        --read1 /data/sample_R1.fastq.gz \
        --read2 /data/sample_R2.fastq.gz \
        --host-index /refs/bowtie2_index/GRCh38_noalt_as \
        --kraken2-db /refs/standard_8 \
        --output-dir /data/results \
        --threads 8
```

### 3. Output Structure

```
results/
├── 01_qc/
│   ├── sample_R1_qc.fastq.gz        # Trimmed reads
│   ├── sample_R1_fastp.json          # QC metrics (JSON)
│   └── sample_R1_fastp.html          # QC report (interactive HTML)
├── 02_alignment/
│   ├── host_aligned.sorted.bam       # Host reads (for scRNA-seq)
│   ├── host_aligned.sorted.bam.bai
│   └── host_align_stats.txt
├── 03_unmapped/
│   └── unmapped_R1.fastq.gz          # Non-host → microbial candidates
├── 04_classification/
│   ├── kraken2_report.txt
│   ├── kraken2_output.txt
│   └── bracken_S.txt                 # Species-level abundance
├── 05_visualisation/
│   └── plots/
│       ├── abundance_bar.png
│       ├── abundance_pie.png
│       ├── abundance_heatmap.png
│       └── pipeline_dashboard.png
└── celljanus.log
```

---

## Alternative Installation

### Conda (Linux / macOS / WSL2)

```bash
conda create -n CellJanus -c bioconda -c conda-forge \
    python=3.11 fastp=1.1.0 bowtie2=2.5.4 samtools=1.23 \
    kraken2=2.17.1 bracken=3.1
conda activate CellJanus

git clone https://github.com/zhaoqing-wang/CellJanus.git
cd CellJanus
pip install .

# Verify
celljanus check

# Run QC on test data
celljanus qc --read1 testdata/reads_R1.fastq.gz \
    --read2 testdata/reads_R2.fastq.gz --output-dir ./qc_output
```

### Windows (WSL2 recommended)

bioconda does not support Windows natively. Use WSL2 with the Conda
instructions above, or use Docker.

<details>
<summary>Windows native (partial — QC & alignment steps only)</summary>

```bash
conda create -n CellJanus python=3.11
conda activate CellJanus
pip install .
```

Download fastp and Bowtie2 binaries manually and add to PATH.
Kraken2/Bracken require Linux — use WSL2 or Docker for the full pipeline.

</details>

---

## Commands

| Command | Purpose |
|---------|---------|
| `celljanus run` | Full pipeline (QC → Align → Extract → Classify → Visualize) |
| `celljanus qc` | Quality control only |
| `celljanus align` | Host alignment only |
| `celljanus extract` | Extract unmapped reads from BAM |
| `celljanus classify` | Microbial classification only |
| `celljanus visualize` | Generate plots from Bracken output |
| `celljanus download` | Download reference data (hg38, Kraken2 DB) |
| `celljanus check` | Verify external tool availability |

### `celljanus run` Options

| Option | Default | Description |
|--------|---------|-------------|
| `--read1` | *required* | R1 FASTQ file |
| `--read2` | — | R2 FASTQ (paired-end) |
| `--host-index` | *required* | Bowtie2 index prefix |
| `--kraken2-db` | *required* | Kraken2 database directory |
| `--output-dir` | `celljanus_output` | Output directory |
| `--threads` | auto | Number of threads |
| `--max-memory` | 80% RAM | Max memory (GB) |
| `--min-quality` | 15 | Phred quality threshold |
| `--min-length` | 36 | Min read length after trim |
| `--confidence` | 0.05 | Kraken2 confidence |
| `--bracken-level` | S | Taxonomic level (D/P/C/O/F/G/S) |
| `--skip-qc` | — | Skip QC step |
| `--skip-classify` | — | Skip classification |
| `--skip-visualize` | — | Skip plots |

---

## Python API

```python
from pathlib import Path
from celljanus.config import CellJanusConfig
from celljanus.pipeline import run_pipeline

cfg = CellJanusConfig(
    output_dir=Path("./results"),
    host_index=Path("./refs/bowtie2_index/GRCh38_noalt_as"),
    kraken2_db=Path("./refs/standard_8"),
    threads=8,
)

result = run_pipeline(
    Path("sample_R1.fastq.gz"),
    read2=Path("sample_R2.fastq.gz"),
    cfg=cfg,
)
print(result.bracken_df.head())
```

---

## Project Structure

```
CellJanus/
├── celljanus/               # Python package
│   ├── __init__.py
│   ├── __main__.py          # python -m celljanus
│   ├── cli.py               # Click CLI
│   ├── config.py            # Configuration + tool discovery
│   ├── utils.py             # Subprocess, FASTQ I/O, logging
│   ├── download.py          # Reference data downloader
│   ├── qc.py                # fastp wrapper
│   ├── align.py             # Bowtie2 alignment
│   ├── extract.py           # samtools unmapped extraction
│   ├── classify.py          # Kraken2 + Bracken
│   ├── visualize.py         # matplotlib / seaborn plots
│   └── pipeline.py          # Pipeline orchestrator
├── testdata/                # Shipped test FASTQ files
│   ├── reads_R1.fastq.gz   # PE R1 (1,000 reads, ~127 KB)
│   ├── reads_R2.fastq.gz   # PE R2
│   ├── reads_SE.fastq.gz   # Single-end
│   └── README.md
├── tests/                   # Test suite (18 tests)
│   ├── conftest.py
│   ├── generate_test_data.py
│   ├── generate_testdata.py # Generates testdata/ files
│   └── test_celljanus.py
├── Dockerfile               # Fully reproducible image
├── docker-compose.yml
├── pyproject.toml
├── CHANGELOG.md
├── LICENSE
└── README.md
```

---

## Performance

| Component | Memory | Note |
|-----------|--------|------|
| fastp | < 1 GB | Streaming |
| Bowtie2 + hg38 | ~3.5 GB | Memory-mapped |
| Kraken2 (8 GB DB) | ~8 GB | `--memory-mapping` |
| **Peak total** | **~12–14 GB** | Fits 32 GB laptop |

| Data Size | Threads | Time |
|-----------|---------|------|
| 10 GB FASTQ | 8 | ~30 min |
| 50 GB FASTQ | 8 | ~2–4 hours |

---

## Roadmap

- [ ] STAR aligner support
- [ ] 10x Genomics barcode-aware demux
- [ ] Spatial transcriptomics integration
- [ ] Interactive HTML reports (Plotly)
- [ ] Nextflow / Snakemake export
- [ ] MetaPhlAn4 classifier support
- [ ] Multi-sample batch processing

---

## Citation

```
Wang Z (2026). CellJanus: A Dual-Perspective Tool for Deconvolving Host Single-Cell
and Microbial Transcriptomes. Python package version 0.1.1.
https://github.com/zhaoqing-wang/CellJanus
```

## License

MIT License. See [LICENSE](LICENSE) for details.

## Acknowledgements

- [fastp](https://github.com/OpenGene/fastp) — Chen et al., *Bioinformatics*, 2018
- [Bowtie2](https://github.com/BenLangmead/bowtie2) — Langmead & Salzberg, *Nature Methods*, 2012
- [samtools](https://github.com/samtools/samtools) — Li et al., *Bioinformatics*, 2009
- [Kraken2](https://github.com/DerrickWood/kraken2) — Wood et al., *Genome Biology*, 2019
- [Bracken](https://github.com/jenniferlu717/Bracken) — Lu et al., *PeerJ Computer Science*, 2017
- [MetaScope](https://github.com/wejlab/metascope) — Odom et al. (conceptual framework)
