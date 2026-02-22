# CellJanus

**A Dual-Perspective Tool for Deconvolving Host Single-Cell and Microbial Transcriptomes**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.9+](https://img.shields.io/badge/Python-3.9%2B-blue.svg)](https://www.python.org/)

CellJanus separates host gene expression from microbial abundance in the same
FASTQ data. It orchestrates proven C/C++ tools (fastp, Bowtie2, samtools,
Kraken2, Bracken) through a single Python CLI.

```
FASTQ ──▶ fastp (QC) ──▶ Bowtie2 (align host) ──▶ unmapped reads
                              │                         │
                              ▼                         ▼
                     host_aligned.bam         Kraken2 + Bracken
                  (gene expression BAM)    (microbial abundance)
```

---

## Table of Contents

- [Installation](#installation)
- [Try It Now (Test Data)](#try-it-now)
- [Full Pipeline (Real Data)](#full-pipeline)
- [Commands](#commands)
- [Python API](#python-api)
- [Project Structure](#project-structure)

---

## Installation

CellJanus requires **Linux / macOS / WSL2** — bioconda tools don't support
native Windows. Install into a dedicated conda environment:

```bash
# 1. Create environment with all external tools
conda create -n CellJanus -c bioconda -c conda-forge \
    python=3.11 fastp bowtie2 samtools kraken2 bracken

# 2. Activate and install CellJanus
conda activate CellJanus
git clone https://github.com/zhaoqing-wang/CellJanus.git
cd CellJanus
pip install .

# 3. Verify — all core tools should show "Found"
celljanus check
```

Expected output:

```
┃ Tool          ┃ Status    ┃ Path                                             ┃
│ fastp         │ ✔ Found   │ /path/to/envs/CellJanus/bin/fastp                │
│ bowtie2       │ ✔ Found   │ /path/to/envs/CellJanus/bin/bowtie2              │
│ samtools      │ ✔ Found   │ /path/to/envs/CellJanus/bin/samtools             │
│ kraken2       │ ✔ Found   │ /path/to/envs/CellJanus/bin/kraken2              │
│ bracken       │ ✔ Found   │ /path/to/envs/CellJanus/bin/bracken              │
```

### Docker (Alternative)

```bash
docker build -t celljanus .
docker run --rm celljanus celljanus check
```

---

## Try It Now

The repository ships test data and reference databases so you can run the
**full pipeline** immediately — no downloads needed.

### Test Data

| File | Reads | Content |
|------|------:|---------|
| `testdata/reads_R1.fastq.gz` | 1,000 | Paired-end R1 |
| `testdata/reads_R2.fastq.gz` | 1,000 | Paired-end R2 |
| `testdata/reads_SE.fastq.gz` | 1,000 | Single-end |

Read composition: 600 human (HBB, ACTB gene fragments), 300 microbial
(*S. aureus*, *E. coli*, *K. pneumoniae* 16S rRNA), 100 low-quality
adapter reads.

### Run Full Pipeline on Test Data

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

### Expected Results

The pipeline completes in **~3 seconds** and produces:

```
┃ Step                  ┃ Result                                         ┃
│ QC (fastp)            │ 1,000 → 900 pairs (90% retained)              │
│                       │ Q20: 88% → 98%                                │
│ Alignment (Bowtie2)   │ 66.39% host alignment rate                    │
│ Unmapped extraction   │ ~300 pairs → microbial candidates             │
│ Kraken2               │ 300 classified (33.33%), 600 unclassified     │
│ Bracken               │ 3 species detected:                           │
│                       │   S. aureus    38.7% (116 reads)              │
│                       │   K. pneumoniae 31.3% (94 reads)              │
│                       │   E. coli      30.0% (90 reads)              │
│ Visualisation         │ 4 plots (bar, pie, heatmap, dashboard)        │
```

### Output Files

```
test_results/
├── 01_qc/
│   ├── reads_R1_qc.fastq.gz         # Trimmed R1
│   ├── reads_R2_qc.fastq.gz         # Trimmed R2
│   ├── reads_R1_fastp.json           # QC metrics
│   └── reads_R1_fastp.html           # Interactive QC report
├── 02_alignment/
│   ├── host_aligned.sorted.bam       # All aligned reads
│   ├── host_aligned.sorted.bam.bai   # BAM index
│   ├── host_mapped.sorted.bam        # Host-only (gene expression)
│   ├── host_mapped.sorted.bam.bai
│   ├── host_align_stats.txt          # Bowtie2 statistics
│   ├── unmapped_R1.fastq.gz          # Non-host R1 (microbial)
│   └── unmapped_R2.fastq.gz          # Non-host R2
├── 04_classification/
│   ├── kraken2_report.txt            # Taxonomic report
│   ├── kraken2_output.txt            # Per-read classification
│   └── bracken_S.txt                 # Species-level abundance
├── 05_visualisation/plots/
│   ├── abundance_bar.png             # Top species bar chart
│   ├── abundance_pie.png             # Community donut chart
│   ├── abundance_heatmap.png         # Abundance heatmap
│   └── pipeline_dashboard.png        # Summary dashboard
└── celljanus.log                     # Full pipeline log
```

### Run Individual Steps

```bash
# Step 1: Quality control only
celljanus qc \
    --read1 testdata/reads_R1.fastq.gz \
    --read2 testdata/reads_R2.fastq.gz \
    --output-dir results/01_qc

# Step 2: Align to host genome
celljanus align \
    --read1 results/01_qc/reads_R1_qc.fastq.gz \
    --read2 results/01_qc/reads_R2_qc.fastq.gz \
    --host-index testdata/refs/host_genome/host \
    --output-dir results/02_alignment

# Step 3: Classify microbial reads
celljanus classify \
    --read1 results/02_alignment/unmapped_R1.fastq.gz \
    --read2 results/02_alignment/unmapped_R2.fastq.gz \
    --kraken2-db testdata/refs/kraken2_testdb \
    --output-dir results/04_classification

# Step 4: Generate plots
celljanus visualize \
    --bracken-file results/04_classification/bracken_S.txt \
    --output-dir results/05_visualisation
```

---

## Full Pipeline

For real data, download reference databases first.

### 1. Download References

```bash
mkdir -p refs

# Human genome hg38 + Bowtie2 index (~5 GB)
celljanus download hg38 --output-dir ./refs

# Kraken2 database (~8 GB standard)
celljanus download kraken2 --output-dir ./refs --db-name standard_8
```

### 2. Run Pipeline

```bash
celljanus run \
    --read1 /path/to/sample_R1.fastq.gz \
    --read2 /path/to/sample_R2.fastq.gz \
    --host-index ./refs/bowtie2_index/GRCh38_noalt_as \
    --kraken2-db ./refs/standard_8 \
    --output-dir ./results \
    --threads 8
```

### 3. Key Options

| Option | Default | Description |
|--------|---------|-------------|
| `--read1` | *required* | R1 FASTQ file |
| `--read2` | — | R2 FASTQ (paired-end) |
| `--host-index` | *required* | Bowtie2 index prefix |
| `--kraken2-db` | *required* | Kraken2 database directory |
| `--output-dir` | `celljanus_output` | Output directory |
| `--threads` | auto | Number of threads |
| `--min-quality` | 15 | Phred quality threshold |
| `--confidence` | 0.05 | Kraken2 confidence |
| `--bracken-level` | S | Taxonomic level (D/P/C/O/F/G/S) |
| `--skip-qc` | — | Skip QC step |
| `--skip-classify` | — | Skip classification |
| `--skip-visualize` | — | Skip plots |

---

## Commands

| Command | Purpose |
|---------|---------|
| `celljanus run` | Full pipeline (QC → Align → Extract → Classify → Visualize) |
| `celljanus qc` | Quality control only (fastp) |
| `celljanus align` | Host alignment + unmapped extraction (Bowtie2) |
| `celljanus extract` | Extract unmapped reads from existing BAM |
| `celljanus classify` | Microbial classification (Kraken2 + Bracken) |
| `celljanus visualize` | Generate plots from Bracken output |
| `celljanus download` | Download reference data (hg38, Kraken2 DB) |
| `celljanus check` | Verify external tool availability |

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

# Access results
print(result.bracken_df.head())       # Species abundance table
print(result.qc_report.summary())     # QC statistics
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
│   ├── align.py             # Bowtie2 alignment + unmapped FASTQ
│   ├── extract.py           # samtools host/unmapped extraction
│   ├── classify.py          # Kraken2 + Bracken
│   ├── visualize.py         # matplotlib / seaborn plots
│   └── pipeline.py          # Pipeline orchestrator
├── testdata/                # Test data + reference databases
│   ├── reads_R1.fastq.gz   # PE R1 (1,000 reads)
│   ├── reads_R2.fastq.gz   # PE R2
│   ├── reads_SE.fastq.gz   # Single-end
│   ├── refs/host_genome/    # Tiny Bowtie2 index
│   ├── refs/kraken2_testdb/ # Tiny Kraken2+Bracken DB
│   ├── build_test_refs.py   # Script to regenerate refs/
│   └── README.md
├── tests/                   # Test suite
│   ├── conftest.py
│   ├── test_celljanus.py    # 18 tests
│   ├── generate_test_data.py
│   └── generate_testdata.py
├── Dockerfile
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

---

## Rebuilding Test References

The shipped `testdata/refs/` contains pre-built databases. To regenerate:

```bash
conda activate CellJanus
cd CellJanus/testdata
python build_test_refs.py    # Downloads NCBI taxonomy (~60 MB)
```

This builds a Bowtie2 index from 2 human gene fragments and a minimal
Kraken2+Bracken database from 3 bacterial 16S sequences.

---

## Citation

```
Wang Z (2026). CellJanus: A Dual-Perspective Tool for Deconvolving Host Single-Cell and Microbial Transcriptomes. Python package version 0.1.2.
https://github.com/zhaoqing-wang/CellJanus
```

## License

MIT License. See [LICENSE](LICENSE) for details.
