# CellJanus

**Dual-perspective deconvolution of host and microbial transcriptomes from FASTQ data.**

[![Version](https://img.shields.io/badge/dynamic/toml?url=https%3A%2F%2Fraw.githubusercontent.com%2Fzhaoqing-wang%2FCellJanus%2Fmain%2Fpyproject.toml&label=Version&query=project.version&color=blue&style=flat-square)](https://github.com/zhaoqing-wang/CellJanus/releases)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.9+](https://img.shields.io/badge/Python-3.9%2B-blue.svg)](https://www.python.org/)

```
FASTQ ─→ fastp (QC) ─→ Bowtie2 (host alignment) ─→ unmapped reads
                             │                            │
                             ▼                            ▼
                    host_aligned.bam           Kraken2 + Bracken
                   (gene expression)        (microbial abundance)
                                                     │
                                                     ▼
                                          plots (PNG/PDF) + CSV tables
```

## Contents

1. [Installation](#installation)
2. [Quick Start](#quick-start)
3. [Real Data](#real-data)
4. [CLI Reference](#cli-reference)
5. [Python API](#python-api)
6. [Output Structure](#output-structure)
7. [Citation](#citation)

---

## Installation

> Requires **Linux / macOS / WSL2**. Bioconda packages are not available on native Windows.

```bash
# Create conda environment with all dependencies
conda create -n CellJanus -c bioconda -c conda-forge \
    python=3.11 fastp bowtie2 samtools kraken2 bracken

# Install CellJanus
conda activate CellJanus
git clone https://github.com/zhaoqing-wang/CellJanus.git
cd CellJanus
pip install .

# Verify installation
celljanus check
```

All tools should show **✔ Found**. STAR is optional (for future RNA-seq alignment support).

<details>
<summary>Docker alternative</summary>

```bash
docker build -t celljanus .
docker run --rm celljanus celljanus check
```
</details>

---

## Quick Start

The repository includes test data and pre-built reference databases — run the full pipeline immediately with **no downloads required**.

```bash
conda activate CellJanus
cd CellJanus

celljanus run \
    --read1 testdata/reads_R1.fastq.gz \
    --read2 testdata/reads_R2.fastq.gz \
    --host-index testdata/refs/host_genome/host \
    --kraken2-db testdata/refs/kraken2_testdb \
    --output-dir test_results \
    --threads 4
```

**Test data**: 1,000 paired-end reads (600 human, 300 microbial, 100 low-quality).

**Results** (~4 seconds):

| Step | Metric |
|------|--------|
| QC | 1,000 → 900 pairs retained (90%), Q20 improved 88% → 98% |
| Host alignment | 66.39% aligned to host genome |
| Classification | 300 reads classified → 3 species detected |
| Top species | *S. aureus* 38.7%, *K. pneumoniae* 31.3%, *E. coli* 30.0% |
| Output | 8 plots (PNG + PDF), 3 CSV tables, QC reports |

Pre-computed results are available in [`test_results/`](test_results/).

### Run Individual Steps

```bash
# QC only
celljanus qc -1 testdata/reads_R1.fastq.gz -2 testdata/reads_R2.fastq.gz -o results/01_qc

# Align to host
celljanus align -1 results/01_qc/reads_R1_qc.fastq.gz \
    -2 results/01_qc/reads_R2_qc.fastq.gz \
    -x testdata/refs/host_genome/host -o results/02_alignment

# Classify microbial reads
celljanus classify -1 results/02_alignment/unmapped_R1.fastq.gz \
    -2 results/02_alignment/unmapped_R2.fastq.gz \
    -d testdata/refs/kraken2_testdb -o results/04_classification

# Generate plots
celljanus visualize -b results/04_classification/bracken_S.txt -o results/05_visualisation
```

---

## Real Data

### 1. Download reference databases

```bash
# Human genome hg38 + Bowtie2 index (~5 GB)
celljanus download hg38 -o ./refs

# Kraken2 standard database (~8 GB)
celljanus download kraken2 -o ./refs --db-name standard_8
```

### 2. Run pipeline

```bash
celljanus run \
    -1 /path/to/sample_R1.fastq.gz \
    -2 /path/to/sample_R2.fastq.gz \
    -x ./refs/bowtie2_index/GRCh38_noalt_as \
    -d ./refs/standard_8 \
    -o ./results \
    --threads 8
```

### Key Options

| Option | Default | Description |
|--------|---------|-------------|
| `-1, --read1` | *required* | R1 FASTQ (or single-end FASTQ) |
| `-2, --read2` | — | R2 FASTQ for paired-end |
| `-x, --host-index` | *required* | Bowtie2 index prefix |
| `-d, --kraken2-db` | *required* | Kraken2 database path |
| `-o, --output-dir` | `celljanus_output` | Output directory |
| `-t, --threads` | auto (CPUs − 2) | Worker threads |
| `--min-quality` | 15 | Phred quality threshold |
| `--confidence` | 0.05 | Kraken2 confidence |
| `--bracken-level` | S | Taxonomic level (D/P/C/O/F/G/S) |
| `--skip-qc` | — | Skip QC step |
| `--skip-classify` | — | Skip classification |
| `--skip-visualize` | — | Skip visualisation |

---

## CLI Reference

| Command | Description |
|---------|-------------|
| `celljanus run` | Full pipeline: QC → Align → Classify → Visualize |
| `celljanus qc` | Quality control (fastp) |
| `celljanus align` | Host alignment + unmapped extraction (Bowtie2) |
| `celljanus extract` | Extract unmapped reads from BAM |
| `celljanus classify` | Taxonomic classification (Kraken2 + Bracken) |
| `celljanus visualize` | Generate abundance plots |
| `celljanus download` | Download reference databases |
| `celljanus check` | Verify external tool installation |

Run `celljanus <command> --help` for full option details.

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

result.bracken_df          # Species abundance (pandas DataFrame)
result.qc_report.summary() # QC statistics
```

---

## Output Structure

```
output_dir/
├── 01_qc/                           # Quality control
│   ├── *_qc.fastq.gz               # Trimmed reads
│   ├── *_fastp.json                 # QC metrics
│   └── *_fastp.html                 # Interactive report
├── 02_alignment/                    # Host alignment
│   ├── host_aligned.sorted.bam      # Full alignment
│   ├── host_mapped.sorted.bam       # Host-only reads
│   ├── unmapped_R{1,2}.fastq.gz     # Non-host reads → classification
│   └── host_align_stats.txt         # Alignment statistics
├── 04_classification/               # Microbial classification
│   ├── kraken2_report.txt           # Taxonomic report
│   ├── kraken2_output.txt           # Per-read assignments
│   └── bracken_S.txt                # Species abundance
├── 05_visualisation/plots/          # Figures (PNG + PDF)
│   ├── abundance_bar.*              # Horizontal bar chart
│   ├── abundance_pie.*              # Donut chart
│   ├── abundance_heatmap.*          # Heatmap (log₁₀ scale)
│   └── pipeline_dashboard.*         # Summary dashboard
├── 06_tables/                       # Machine-readable results
│   ├── pipeline_summary.csv         # Per-step metrics
│   ├── species_abundance.csv        # Species × reads × fraction
│   └── output_manifest.csv          # File inventory with sizes
└── celljanus.log                    # Pipeline log
```

### CSV Tables

**`species_abundance.csv`**:

| name | taxonomy_id | bracken_estimated | fraction_pct |
|------|-------------|------------------:|--------------:|
| Staphylococcus aureus | 1280 | 116 | 38.67 |
| Klebsiella pneumoniae | 573 | 94 | 31.33 |
| Escherichia coli | 562 | 90 | 30.00 |

**`pipeline_summary.csv`**: one row per metric (Step, Metric, Value) covering QC, alignment, and classification statistics.

---

## Performance

| Component | Memory | Note |
|-----------|--------|------|
| fastp | < 1 GB | Streaming I/O |
| Bowtie2 + hg38 | ~3.5 GB | Memory-mapped index |
| Kraken2 (standard DB) | ~8 GB | `--memory-mapping` flag |
| **Peak total** | **~12–14 GB** | Fits a 32 GB laptop |

---

## Citation

```
Wang Z (2026). CellJanus: A Dual-Perspective Tool for Deconvolving Host
Single-Cell and Microbial Transcriptomes. Python package version 0.1.3.
https://github.com/zhaoqing-wang/CellJanus
```

## License

[MIT](LICENSE)
