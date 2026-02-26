<table>
  <tr>
    <td>
      <h1>CellJanus</h1>
      <p><b>Dual-Perspective Deconvolution of Host and Microbial Transcriptomes from FASTQ Data</b></p>
      <p>
        <a href="https://pypi.org/project/celljanus/"><img src="https://img.shields.io/pypi/v/celljanus?color=blue&style=flat-square" alt="PyPI Version" /></a>
        <a href="LICENSE"><img src="https://img.shields.io/badge/License-MIT-green.svg" alt="License: MIT" /></a>
        <a href="https://www.python.org/"><img src="https://img.shields.io/badge/Python-3.9%2B-blue.svg" alt="Python 3.9+" /></a>
        <a href="https://github.com/zhaoqing-wang"><img src="https://img.shields.io/badge/Maintainer-Zhaoqing_Wang-green" alt="GitHub Maintainer" /></a>
      </p>
    </td>
    <td width="200">
      <img src="docs/Sticker.png" alt="CellJanus Logo" width="200" />
    </td>
  </tr>
</table>

## Features

- **Dual-mode analysis**: Bulk RNA-seq and scRNA-seq pipelines
- **Per-cell tracking**: Cell barcode (CB) and UMI extraction for 10x Genomics, Parse Biosciences
- **Automated pipeline**: QC → Host alignment → Microbial classification → Visualization
- **Comprehensive exports**: 6 CSV formats, 8+ publication-ready plots (PNG/PDF)
- **Large dataset support**: Optimized for 10,000+ cells with smart sampling
- **WSL2 optimized**: Cross-filesystem detection and performance warnings

---

## Pipeline Overview

**Bulk RNA-seq:**
```
FASTQ ─→ fastp (QC) ─→ Bowtie2 (host) ─→ unmapped reads ─→ Kraken2+Bracken ─→ plots + CSV
```

**scRNA-seq:**
```
10x FASTQ ─→ Extract CB+UMI ─→ Kraken2 ─→ Per-cell abundance ─→ Cell×Species matrix
```

---

## Contents

1. [Installation](#installation)
2. [Bulk RNA-seq Mode](#bulk-rna-seq-mode)
3. [scRNA-seq Mode](#scrna-seq-mode)
4. [CLI Reference](#cli-reference)
5. [Python API](#python-api)
6. [Output Structure](#output-structure)
7. [Advanced](#advanced)

---

## Installation

### One-Click Setup (Recommended)

Creates a complete environment with CellJanus **and all external tools** (fastp, Bowtie2, samtools, Kraken2, Bracken):

```bash
# Option 1: From GitHub (recommended)
conda env create -f https://raw.githubusercontent.com/zhaoqing-wang/CellJanus/main/environment.yml
conda activate celljanus
celljanus check   # All tools should show ✔ Found

# Option 2: Clone and install
git clone https://github.com/zhaoqing-wang/CellJanus.git
cd CellJanus
conda env create -f environment.yml
conda activate celljanus
celljanus check
```

> **Requirements**: [Conda](https://docs.conda.io/en/latest/miniconda.html) or [Mamba](https://mamba.readthedocs.io/). Works on **Linux / macOS / WSL2**.

<details>
<summary><b>Alternative methods</b></summary>

**pip only** (requires external tools already on PATH):
```bash
pip install celljanus
```

**Development install:**
```bash
git clone https://github.com/zhaoqing-wang/CellJanus.git
cd CellJanus && pip install -e ".[dev]"
```

**Docker:**
```bash
docker build -t celljanus . && docker run --rm celljanus celljanus check
```

**WSL2 (Ubuntu) with APT:**
```bash
# Install external tools
sudo apt-get update
sudo apt-get install -y fastp bowtie2 samtools kraken2

# Install Bracken from source
git clone https://github.com/jenniferlu717/Bracken.git /opt/Bracken
cd /opt/Bracken && bash install_bracken.sh
sudo ln -sf /opt/Bracken/bracken /usr/local/bin/bracken

# Install CellJanus
pip install celljanus --break-system-packages
celljanus check
```
</details>

---

## Bulk RNA-seq Mode

Full pipeline: QC → Host alignment → Microbial classification → Visualization.

### Quick Test

Run immediately with built-in test data — **no downloads required**:

```bash
celljanus run \
    --read1 testdata/reads_R1.fastq.gz \
    --read2 testdata/reads_R2.fastq.gz \
    --host-index testdata/refs/host_genome/host \
    --kraken2-db testdata/refs/kraken2_testdb \
    --output-dir test_results/bulk
```

**Test Results** (~4 seconds):

| Metric | Value |
|--------|-------|
| Input reads | 1,000 |
| QC-passed | 950 (95%) |
| Species detected | 5 |

| Species | Reads | Fraction |
|---------|------:|----------|
| *Staphylococcus aureus* | 109 | 37.1% |
| *Escherichia coli* | 79 | 26.9% |
| *Klebsiella pneumoniae* | 56 | 19.0% |
| *Acetitomaculum* | 31 | 10.5% |
| *Longispora* | 19 | 6.5% |

**Example Output:**

| Pipeline Dashboard |
|:--:|
| ![Dashboard](docs/pipeline_dashboard.png) |

| Abundance Bar | Abundance Pie | Abundance Heatmap |
|:--:|:--:|:--:|
| ![Bar](docs/abundance_bar.png) | ![Pie](docs/abundance_pie.png) | ![Heatmap](docs/abundance_heatmap.png) |

### Real Data

```bash
# 1. Download references (~13 GB total)
celljanus download hg38 -o ./refs           # Human genome + Bowtie2 index (~5 GB)
celljanus download kraken2 -o ./refs        # Kraken2 standard DB (~8 GB)

# 2. Run full pipeline (paired-end)
celljanus run \
    -1 sample_R1.fastq.gz \                  # Forward reads (required)
    -2 sample_R2.fastq.gz \                  # Reverse reads (optional for PE)
    -x ./refs/bowtie2_index/GRCh38_noalt_as \ # Bowtie2 index prefix
    -d ./refs/standard_8 \                   # Kraken2 database path
    -o ./results \                           # Output directory
    --threads 8                              # Number of threads (default: auto)
```

<details>
<summary><b>Advanced bulk options</b></summary>

```bash
# Run with custom QC and classification parameters
celljanus run \
    -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
    -x ./refs/bowtie2_index/GRCh38_noalt_as \
    -d ./refs/standard_8 \
    -o ./results \
    --threads 8 \
    --min-quality 20 \                       # Phred quality threshold (default: 15)
    --min-length 50 \                        # Min read length after trimming (default: 36)
    --confidence 0.1 \                       # Kraken2 confidence threshold (default: 0.05)
    --bracken-level G \                      # Taxonomic level: S=species, G=genus (default: S)
    --plot-format pdf \                      # Plot format: png, pdf, svg (default: png)
    --max-memory 16                          # Max memory in GB (default: 80% of RAM)

# Skip specific steps for partial re-runs
celljanus run \
    -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
    -x ./refs/bowtie2_index/GRCh38_noalt_as \
    -d ./refs/standard_8 \
    -o ./results \
    --skip-qc \                              # Skip QC if already trimmed
    --skip-visualize                         # Skip plot generation

# Single-end mode (omit --read2)
celljanus run \
    -1 sample_SE.fastq.gz \
    -x ./refs/bowtie2_index/GRCh38_noalt_as \
    -d ./refs/standard_8 \
    -o ./results_se
```

</details>

---

## scRNA-seq Mode

Per-cell microbial abundance tracking with barcode extraction.

### Quick Test

Run immediately with built-in test data — **no downloads required**:

```bash
celljanus scrnaseq \
    --read1 testdata/scrnaseq/scrna_R1.fastq.gz \
    --read2 testdata/scrnaseq/scrna_R2.fastq.gz \
    --kraken2-db testdata/refs/kraken2_testdb \
    --output-dir test_results/scrnaseq \
    --barcode-mode 10x \
    --min-reads 1
```

**Test Results** (~2 seconds):

| Metric | Value |
|--------|-------|
| Input reads | 15,000 |
| Cells processed | 300 |
| Species detected | 7 |
| Microbial reads | 2,395 |
| Mean reads/cell | 8.0 |

| Species | Reads | Cells | Prevalence |
|---------|------:|------:|-----------:|
| *Escherichia coli* | 460 | 235 | 78.3% |
| *Prevotella* | 426 | 234 | 78.0% |
| *Staphylococcus aureus* | 418 | 218 | 72.7% |
| *Klebsiella pneumoniae* | 405 | 225 | 75.0% |
| *Acetitomaculum* | 311 | 195 | 65.0% |
| *Longispora* | 228 | 134 | 44.7% |
| *Mobiluncus* | 147 | 108 | 36.0% |

**Example Output:**

| scRNA-seq Dashboard |
|:--:|
| ![scRNA-seq Dashboard](docs/scrnaseq_dashboard.png) |

| Species Abundance Pie | Microbial Summary |
|:--:|:--:|
| ![Abundance Pie](docs/scrnaseq_abundance_pie.png) | ![Summary](docs/cell_microbe_summary.png) |

### Real Data

#### Step 1: Download Kraken2 Database

scRNA-seq mode classifies microbial reads directly (no host alignment needed), so only a Kraken2 database is required:

```bash
# Download pre-built Kraken2 database (~8 GB)
celljanus download kraken2 -o ./refs         # Default: standard_8 database

# Or choose a different database variant
celljanus download kraken2 -o ./refs -n pluspfp_8  # Includes protozoa, fungi, plant
```

#### Step 2: Run scRNA-seq Pipeline

**10x Genomics data (most common):**

```bash
# Basic 10x Genomics paired-end run
celljanus scrnaseq \
    --read1 sample_R1.fastq.gz \             # R1: contains cell barcodes + UMI
    --read2 sample_R2.fastq.gz \             # R2: contains cDNA sequences
    --kraken2-db ./refs/standard_8 \         # Kraken2 database path
    --output-dir scrna_results \             # Output directory
    --barcode-mode 10x \                     # Barcode format: 10x, parse, or auto
    --threads 8                              # Number of threads (default: auto)

# With barcode whitelist filtering (recommended for real data)
# Whitelist files are typically provided by 10x Genomics, e.g.:
#   - 3M-february-2018.txt.gz (v3 chemistry)
#   - 737K-august-2016.txt   (v2 chemistry)
celljanus scrnaseq \
    --read1 sample_R1.fastq.gz \
    --read2 sample_R2.fastq.gz \
    --kraken2-db ./refs/standard_8 \
    --output-dir scrna_results_filtered \
    --barcode-mode 10x \
    --whitelist 3M-february-2018.txt.gz \    # Cell barcode whitelist (optional)
    --min-reads 5 \                          # Min reads per cell to include (default: 1)
    --threads 8
```

**Parse Biosciences data:**

```bash
# Parse Biosciences scRNA-seq (barcode in read name, colon-separated)
celljanus scrnaseq \
    --read1 parse_R1.fastq.gz \
    --read2 parse_R2.fastq.gz \
    --kraken2-db ./refs/standard_8 \
    --output-dir scrna_parse_results \
    --barcode-mode parse \                   # Parse Biosciences barcode format
    --min-reads 3 \
    --threads 8
```

**Auto-detect barcode format:**

```bash
# Let CellJanus auto-detect the barcode format
celljanus scrnaseq \
    --read1 sample_R1.fastq.gz \
    --read2 sample_R2.fastq.gz \
    --kraken2-db ./refs/standard_8 \
    --output-dir scrna_auto_results \
    --barcode-mode auto                      # Auto-detect: tries 10x → Parse → other
```

#### Step 3: Downstream Integration (Seurat / Scanpy)

The output CSV files integrate directly with standard scRNA-seq frameworks:

```r
# R / Seurat: Load normalized microbial abundance as an assay
library(Seurat)

# Load your existing Seurat object
sc <- readRDS("seurat_object.rds")

# Read CellJanus normalized matrix (cells × species, CPM-normalized)
microbe_mat <- read.csv("scrna_results/tables/cell_species_normalized.csv",
                        row.names = 1, check.names = FALSE)
microbe_mat <- t(as.matrix(microbe_mat))            # Seurat expects features × cells

# Add as a new assay
sc[["Microbe"]] <- CreateAssayObject(counts = microbe_mat)
```

```python
# Python / Scanpy: Load microbial abundance as an AnnData layer
import scanpy as sc
import pandas as pd

# Load your existing AnnData object
adata = sc.read_h5ad("adata.h5ad")

# Read CellJanus normalized matrix
microbe_df = pd.read_csv("scrna_results/tables/cell_species_normalized.csv",
                         index_col=0)

# Align cells and add as obsm
common_cells = adata.obs_names.intersection(microbe_df.index)
adata_sub = adata[common_cells, :]
adata_sub.obsm["X_microbe"] = microbe_df.loc[common_cells].values
```

### Supported Platforms

| Platform | Mode | Barcode Location |
|----------|------|------------------|
| 10x Genomics | `10x` | Header: `CB:Z:BARCODE UB:Z:UMI` |
| Parse Biosciences | `parse` | Read name (colon-separated) |
| Custom | `auto` | Auto-detect |

### Output Files

| File | Description |
|------|-------------|
| `cell_species_counts.csv` | Raw counts (cells × species) |
| `cell_species_normalized.csv` | CPM-normalized for Seurat/Scanpy |
| `cell_species_long.csv` | Tidy format for ggplot2/seaborn |
| `species_summary.csv` | Per-species statistics |
| `cell_summary.csv` | Per-cell diversity metrics |

---

## CLI Reference

| Command | Description |
|---------|-------------|
| `celljanus run` | Full pipeline: QC → Align → Classify → Visualize |
| `celljanus scrnaseq` | scRNA-seq mode with per-cell tracking |
| `celljanus qc` | Quality control (fastp) |
| `celljanus align` | Host alignment (Bowtie2) |
| `celljanus extract` | Extract unmapped (non-host) reads from BAM |
| `celljanus classify` | Taxonomic classification (Kraken2 + Bracken) |
| `celljanus visualize` | Generate abundance plots |
| `celljanus download` | Download reference databases |
| `celljanus check` | Verify tool installation |

### Key Options

| Option | Default | Description |
|--------|---------|-------------|
| `-1, --read1` | *required* | R1 FASTQ file |
| `-2, --read2` | — | R2 FASTQ (paired-end) |
| `-x, --host-index` | *required* | Bowtie2 index prefix |
| `-d, --kraken2-db` | *required* | Kraken2 database path |
| `-o, --output-dir` | `celljanus_output` | Output directory |
| `-t, --threads` | auto | Worker threads |
| `--min-quality` | 15 | Phred quality (bulk QC) |
| `--confidence` | 0.05 | Kraken2 confidence threshold |
| `--bracken-level` | S | Taxonomic level (D/P/C/O/F/G/S) |
| `--barcode-mode` | 10x | scRNA-seq barcode format (10x/parse/auto) |
| `-w, --whitelist` | — | Cell barcode whitelist (scRNA-seq) |
| `--min-reads` | 1 | Min reads per cell (scRNA-seq) |
| `--skip-qc` | — | Skip QC step (bulk) |
| `--plot-format` | png | Plot format: png, pdf, svg |

Run `celljanus <command> --help` for full details.

---

## Python API

### Bulk Pipeline

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

result = run_pipeline(Path("sample_R1.fastq.gz"), read2=Path("sample_R2.fastq.gz"), cfg=cfg)
result.bracken_df           # Species abundance DataFrame
result.qc_report.summary()  # QC statistics
```

### scRNA-seq Pipeline

```python
from pathlib import Path
from celljanus.config import CellJanusConfig
from celljanus.scrnaseq import BarcodeConfig, run_scrnaseq_classification

# Configure pipeline
cfg = CellJanusConfig(
    output_dir=Path("./scrna_results"),
    kraken2_db=Path("./refs/standard_8"),
    threads=8,
)

# Configure barcode extraction
barcode_cfg = BarcodeConfig(
    mode="10x",                              # Barcode format: 10x, parse, auto
    min_reads_per_cell=5,                     # Filter low-quality cells
    # whitelist_path=Path("barcodes.txt"),    # Optional: barcode whitelist
)

# Run scRNA-seq classification
result = run_scrnaseq_classification(
    Path("sample_R1.fastq.gz"),              # R1 FASTQ (barcodes)
    Path("./refs/standard_8"),               # Kraken2 database
    Path("./scrna_results"),                 # Output directory
    read2=Path("sample_R2.fastq.gz"),        # R2 FASTQ (cDNA)
    barcode_cfg=barcode_cfg,
    cfg=cfg,
)

# Access per-cell abundance data
abundance = result["abundance"]
abundance.to_matrix()             # Raw counts (cells × species DataFrame)
abundance.to_normalized_matrix()  # CPM-normalized (ready for Seurat/Scanpy)
abundance.to_long_format()        # Tidy format (for ggplot2/seaborn)
abundance.to_species_summary()    # Per-species statistics
abundance.to_cell_summary()       # Per-cell diversity metrics
abundance.summary()               # Pipeline summary dict

# Access output file paths
print(result["matrix_path"])          # cell_species_counts.csv
print(result["normalized_path"])      # cell_species_normalized.csv
print(result["long_path"])            # cell_species_long.csv
print(result["species_summary_path"]) # species_summary.csv
print(result["cell_summary_path"])    # cell_summary.csv
```

---

## Output Structure

### Bulk RNA-seq

```
test_results/bulk/
├── 01_qc/                        # Quality control (fastp)
│   ├── *_qc.fastq.gz            # Trimmed reads
│   ├── *_fastp.json              # QC metrics (machine-readable)
│   └── *_fastp.html              # QC report (interactive)
├── 02_alignment/                 # Host alignment (Bowtie2)
│   ├── host_aligned.sorted.bam   # Full alignment
│   ├── host_mapped.sorted.bam    # Host-only reads
│   ├── unmapped_R{1,2}.fastq.gz  # Non-host → classification
│   └── host_align_stats.txt      # Alignment statistics
├── 04_classification/            # Microbial classification
│   ├── kraken2_report.txt        # Taxonomic report
│   ├── kraken2_output.txt        # Per-read assignments
│   └── bracken_S.txt             # Species-level abundance
├── 05_visualisation/plots/       # Figures (4 PNG + 4 PDF)
│   ├── abundance_bar.*           # Top species bar chart
│   ├── abundance_pie.*           # Donut chart
│   ├── abundance_heatmap.*       # Heatmap (log₁₀)
│   └── pipeline_dashboard.*      # Summary dashboard
├── 06_tables/                    # Machine-readable results
│   ├── species_abundance.csv     # Species × reads × fraction
│   ├── pipeline_summary.csv      # Per-step metrics
│   └── output_manifest.csv       # File inventory
└── celljanus.log                 # Pipeline log
```

### scRNA-seq

```
test_results/scrnaseq/
├── classification/               # Kraken2 output
│   ├── kraken2_report.txt        # Taxonomic report
│   └── kraken2_output.txt        # Per-read assignments
├── tables/                       # 6 CSV exports
│   ├── cell_species_counts.csv   # Raw counts (cells × species)
│   ├── cell_species_normalized.csv # CPM-normalized (Seurat/Scanpy)
│   ├── cell_species_long.csv     # Tidy format (ggplot2/seaborn)
│   ├── species_summary.csv       # Per-species statistics
│   ├── cell_summary.csv          # Per-cell diversity metrics
│   └── pipeline_summary.csv      # Pipeline metrics
├── visualisation/plots/          # Figures (3 PNG + 3 PDF)
│   ├── scrnaseq_dashboard.*      # Summary dashboard
│   ├── scrnaseq_abundance_pie.*  # Community composition donut
│   └── cell_microbe_summary.*    # Distribution panels (3-panel)
└── celljanus.log                 # Pipeline log
```

---

## Advanced

### WSL2 Optimization

For best performance on WSL2, store data on the Linux filesystem:

```bash
# Copy to Linux filesystem (10-50× faster I/O)
mkdir -p ~/celljanus_work
cp /mnt/c/Data/sample*.fastq.gz ~/celljanus_work/

# Run from Linux path
celljanus run -1 ~/celljanus_work/sample_R1.fastq.gz ...
```

CellJanus auto-detects WSL2 and warns about slow cross-filesystem paths.

---

## Citation

```
Wang Z (2026). CellJanus: Dual-Perspective Deconvolution of Host and Microbial Transcriptomes from FASTQ Data.
https://github.com/zhaoqing-wang/CellJanus
```

## License

[MIT](LICENSE)

## Contact

**Author:** Zhaoqing Wang ([ORCID](https://orcid.org/0000-0001-8348-7245)) | **Email:** <zhaoqingwang@mail.sdu.edu.cn> | **Issues:** [CellJanus Issues](https://github.com/zhaoqing-wang/CellJanus/issues)
