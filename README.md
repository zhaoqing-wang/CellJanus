# CellJanus

**A Dual-Perspective Tool for Deconvolving Host Single-Cell and Microbial Transcriptomes**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.9+](https://img.shields.io/badge/Python-3.9%2B-blue.svg)](https://www.python.org/)
[![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg)](https://hub.docker.com/)

<p align="center">
  <a href="#highlights">Highlights</a> &bull;
  <a href="#pipeline-overview">Pipeline</a> &bull;
  <a href="#installation">Installation</a> &bull;
  <a href="#quick-start">Quick Start</a> &bull;
  <a href="#commands">Commands</a> &bull;
  <a href="#performance">Performance</a> &bull;
  <a href="#citation">Citation</a>
</p>

---

## Highlights

| Feature | Description |
|---------|-------------|
| **Dual-perspective deconvolution** | Simultaneously extracts host gene expression and microbial abundance from the same FASTQ |
| **CellRanger-like CLI** | Simple `celljanus run` command with sensible defaults |
| **Gentle QC** | fastp with permissive settings (Q15 / 36 bp) preserves low-abundance microbial reads |
| **High-performance tools** | Python orchestrates C/C++ backends: Bowtie2, samtools, Kraken2, Bracken |
| **Memory-efficient** | Handles 50 GB FASTQ on a 32 GB laptop — streaming I/O, memory-mapped DBs |
| **Docker-ready** | One command to build a fully reproducible environment |
| **Publication-quality plots** | Bar charts, donut plots, heatmaps of microbial abundance |

---

## Pipeline Overview

```
FASTQ ──▶ fastp (QC) ──▶ Bowtie2 (align hg38) ──▶ samtools (extract unmapped)
                                                         │
                                 ┌───────────────────────┘
                                 ▼
                          Kraken2 + Bracken ──▶ Abundance tables + Plots
```

| Step | Tool | Purpose | Output |
|------|------|---------|--------|
| 1. QC | [fastp](https://github.com/OpenGene/fastp) | Adapter trimming, quality filtering (gentle Q15) | `*_qc.fastq.gz`, HTML/JSON report |
| 2. Align | [Bowtie2](https://github.com/BenLangmead/bowtie2) | Align reads to human genome (hg38) | `host_aligned.sorted.bam` |
| 3. Extract | [samtools](https://github.com/samtools/samtools) | Separate host-mapped and unmapped reads | `unmapped_R1.fastq.gz` |
| 4. Classify | [Kraken2](https://github.com/DerrickWood/kraken2) + [Bracken](https://github.com/jenniferlu717/Bracken) | Taxonomic classification + abundance re-estimation | `bracken_S.txt` |
| 5. Visualize | matplotlib + seaborn | Abundance bar charts, donut plots, heatmaps | PNG/PDF/SVG |

---

## Installation

### Option A: Docker (Recommended — Any Platform)

Docker provides a fully self-contained environment with all tools pre-installed.
No conda or manual binary downloads needed.

```bash
# Build the image
git clone https://github.com/zhaoqing-wang/CellJanus.git
cd CellJanus
docker build -t celljanus:0.1.0 .

# Verify
docker run --rm celljanus:0.1.0 celljanus check
```

To run the pipeline, mount your data directory:

```bash
docker run --rm -v /path/to/data:/data celljanus:0.1.0 \
    celljanus run \
        --read1 /data/sample_R1.fastq.gz \
        --read2 /data/sample_R2.fastq.gz \
        --host-index /data/refs/bowtie2_index/GRCh38_noalt_as \
        --kraken2-db /data/refs/standard_8 \
        --output-dir /data/results \
        --threads 8
```

Or use **docker compose**:

```bash
DATA_DIR=/path/to/data docker compose run --rm celljanus \
    celljanus run --read1 /data/R1.fastq.gz --read2 /data/R2.fastq.gz ...
```

### Option B: Conda — Linux / macOS / WSL2

One-step install via bioconda:

```bash
# Create environment with all external tools
conda create -n CellJanus -c bioconda -c conda-forge \
    python=3.11 \
    fastp=1.1.0 \
    bowtie2=2.5.4 \
    samtools=1.23 \
    kraken2=2.17.1 \
    bracken=3.1

conda activate CellJanus

# Install CellJanus
git clone https://github.com/zhaoqing-wang/CellJanus.git
cd CellJanus
pip install .
```

### Option C: Windows Native (Two-Step)

bioconda packages are Linux/macOS only. On Windows native, install the Python
package first, then add external tool binaries to your PATH manually.

**Step 1 — Python package**

```bash
conda create -n CellJanus python=3.11
conda activate CellJanus
git clone https://github.com/zhaoqing-wang/CellJanus.git
cd CellJanus
pip install .
```

**Step 2 — External tools**

| Tool | Source | Notes |
|------|--------|-------|
| fastp | [GitHub Releases](https://github.com/OpenGene/fastp/releases) | `fastp.exe` available |
| Bowtie2 | [GitHub Releases](https://github.com/BenLangmead/bowtie2/releases) | MinGW `.exe` available |
| samtools | [MSYS2](https://packages.msys2.org/package/mingw-w64-x86_64-samtools) | No official Windows binary |
| Kraken2 | **WSL2 only** | No Windows build |
| Bracken | **WSL2 only** | No Windows build |

> **Tip:** For full pipeline functionality on Windows, use **WSL2** (Option B inside
> your WSL2 distribution) or **Docker** (Option A).

### Verify Installation

```bash
celljanus check
# or: python -m celljanus check
```

<details>
<summary><strong>Troubleshooting</strong></summary>

| Problem | Solution |
|---------|----------|
| `celljanus` command not found | Use `python -m celljanus check`, or add conda `Scripts/` to PATH |
| `pysam` / `htslib` build errors | CellJanus does **not** depend on `pysam` — unrelated |
| `conda install fastp` fails on Windows | bioconda is Linux/macOS only — use manual binary or WSL2 |
| STAR missing in `celljanus check` | STAR is optional (roadmap), not required for the core pipeline |

</details>

---

## Quick Start

### 1. Download Reference Data

```bash
# Human genome hg38 + pre-built Bowtie2 index (~5 GB)
celljanus download hg38 --output-dir ./refs

# Kraken2 database (~8 GB, fits in 32 GB RAM)
celljanus download kraken2 --output-dir ./refs --db-name standard_8
```

### 2. Run Full Pipeline

```bash
celljanus run \
    --read1 sample_R1.fastq.gz \
    --read2 sample_R2.fastq.gz \
    --host-index ./refs/bowtie2_index/GRCh38_noalt_as \
    --kraken2-db ./refs/standard_8 \
    --output-dir ./results \
    --threads 8
```

### 3. Output Structure

```
results/
├── 01_qc/
│   ├── sample_R1_qc.fastq.gz
│   ├── sample_R1_fastp.json
│   └── sample_R1_fastp.html
├── 02_alignment/
│   ├── host_aligned.sorted.bam
│   ├── host_aligned.sorted.bam.bai
│   └── host_align_stats.txt
├── 03_unmapped/
│   └── unmapped_R1.fastq.gz
├── 04_classification/
│   ├── kraken2_report.txt
│   ├── kraken2_output.txt
│   └── bracken_S.txt
├── 05_visualisation/
│   └── plots/
│       ├── abundance_bar.png
│       ├── abundance_pie.png
│       ├── abundance_heatmap.png
│       └── pipeline_dashboard.png
└── celljanus.log
```

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
| `--read2` | — | R2 FASTQ file (paired-end) |
| `--host-index` | *required* | Bowtie2 index prefix for hg38 |
| `--kraken2-db` | *required* | Kraken2 database directory |
| `--output-dir` | `celljanus_output` | Output directory |
| `--threads` | auto (CPUs − 2) | Number of threads |
| `--max-memory` | 80% of RAM | Max memory usage (GB) |
| `--min-quality` | 15 | Phred quality threshold |
| `--min-length` | 36 | Minimum read length after trim |
| `--confidence` | 0.05 | Kraken2 confidence threshold |
| `--bracken-level` | S | Taxonomic level: D/P/C/O/F/G/S |
| `--skip-qc` | — | Skip QC step |
| `--skip-classify` | — | Skip classification step |
| `--skip-visualize` | — | Skip visualisation step |

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

## Performance

### Memory Budget (32 GB laptop)

| Component | Memory | Strategy |
|-----------|--------|----------|
| fastp | < 1 GB | Streaming I/O |
| Bowtie2 + hg38 | ~3.5 GB | Memory-mapped index |
| samtools | ~1–2 GB | Streaming BAM |
| Kraken2 (8 GB DB) | ~8 GB | `--memory-mapping` |
| Bracken | < 500 MB | Post-processing |
| **Peak total** | **~12–14 GB** | Fits comfortably |

### Speed Estimates

| Data Size | Threads | Full Pipeline |
|-----------|---------|---------------|
| 10 GB FASTQ | 8 | ~30 min |
| 50 GB FASTQ | 8 | ~2–4 hours |

---

## Project Structure

```
CellJanus/
├── celljanus/               # Package source
│   ├── __init__.py          # Package metadata
│   ├── __main__.py          # python -m celljanus
│   ├── cli.py               # Click CLI with subcommands
│   ├── config.py            # Configuration dataclass + tool discovery
│   ├── utils.py             # Subprocess runner, FASTQ I/O, logging
│   ├── download.py          # Reference data downloader
│   ├── qc.py                # fastp wrapper
│   ├── align.py             # Bowtie2 alignment
│   ├── extract.py           # samtools unmapped extraction
│   ├── classify.py          # Kraken2 + Bracken classification
│   ├── visualize.py         # matplotlib / seaborn plots
│   └── pipeline.py          # Full pipeline orchestrator
├── tests/                   # Test suite
│   ├── conftest.py          # pytest fixtures
│   ├── generate_test_data.py
│   └── test_celljanus.py    # 18 tests
├── Dockerfile               # Docker image definition
├── docker-compose.yml       # Docker Compose convenience
├── pyproject.toml           # Build config + dependencies
├── CHANGELOG.md             # Release notes
├── LICENSE                  # MIT License
└── README.md                # This file
```

---

## Comparison with MetaScope (R)

| Feature | MetaScope (R) | CellJanus (Python) |
|---------|--------------|-------------------|
| QC | Manual / external | fastp (integrated, gentle Q15) |
| Aligner | Rbowtie2 / Rsubread | Bowtie2 direct (piped to samtools) |
| BAM handling | Rsamtools | samtools CLI (streaming) |
| Microbial ID | Bayesian EM in R | Kraken2 + Bracken (C++) |
| Speed (10 GB) | Hours | ~30 min |
| Memory | High (R in-memory) | Low (streaming + memory-mapping) |
| Deployment | R scripts | CLI + Docker |

---

## Roadmap

- [ ] STAR aligner support for splice-aware host alignment
- [ ] 10x Genomics barcode-aware demultiplexing
- [ ] Spatial transcriptomics coordinate integration
- [ ] Interactive HTML reports (Plotly)
- [ ] Nextflow / Snakemake workflow export
- [ ] MetaPhlAn4 integration as alternative classifier
- [ ] Multi-sample batch processing with merged reports

---

## Citation

If you use CellJanus in your research, please cite:

```
Wang Z (2026). CellJanus: A Dual-Perspective Tool for Deconvolving Host Single-Cell
and Microbial Transcriptomes. Python package version 0.1.0.
https://github.com/zhaoqing-wang/CellJanus
```

---

## License

MIT License. See [LICENSE](LICENSE) for details.

## Acknowledgements

CellJanus builds on the conceptual framework of
[MetaScope](https://github.com/wejlab/metascope) (Odom et al.) and integrates
established bioinformatics tools:

- [fastp](https://github.com/OpenGene/fastp) — Chen et al., *Bioinformatics*, 2018
- [Bowtie2](https://github.com/BenLangmead/bowtie2) — Langmead & Salzberg, *Nature Methods*, 2012
- [samtools](https://github.com/samtools/samtools) — Li et al., *Bioinformatics*, 2009
- [Kraken2](https://github.com/DerrickWood/kraken2) — Wood et al., *Genome Biology*, 2019
- [Bracken](https://github.com/jenniferlu717/Bracken) — Lu et al., *PeerJ Computer Science*, 2017
