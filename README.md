<table>
  <tr>
    <td>
      <h1>CellJanus: Dual-Perspective Deconvolution of Host and Microbial Transcriptomes from FASTQ Data</h1>
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

## Pipeline

### Bulk RNA-seq Mode

```
FASTQ ─→ fastp (QC) ─→ Bowtie2 (host alignment) ─→ unmapped reads
                             │                            │
                             ▼                            ▼
                    host_aligned.bam           Kraken2 + Bracken  ─→  plots (PNG/PDF) + CSV tables
                   (gene expression)        (microbial abundance)
```

### scRNA-seq Mode

```
10x FASTQ ─→ Extract CB+UMI ─→ Kraken2 classification ─→ Per-cell abundance
    │              │                    │                        │
    ▼              ▼                    ▼                        ▼
 R1+R2 reads   Cell barcode     Species per read      Cell × Species matrix
               + UMI tags                              Heatmaps + dot plots
```

**Features**: Per-cell barcode tracking • UMI deduplication • 6 CSV exports • Large dataset support (10,000+ cells)

## Contents

1. [Installation](#installation)
2. [Quick Start](#quick-start)
3. [scRNA-seq Mode](#scrna-seq-mode)
4. [Real Data](#real-data)
5. [CLI Reference](#cli-reference)
6. [Python API](#python-api)
7. [Output Structure](#output-structure)
8. [WSL2 Performance](#wsl2-performance-tips)
9. [STAR Installation](#star-installation-wsl2linux)
10. [Citation](#citation)

---

## Installation

### Quick Setup (recommended)

One command creates a complete environment with CellJanus **and** all external tools (fastp, Bowtie2, samtools, Kraken2, Bracken):

> Requires [Conda](https://docs.conda.io/en/latest/miniconda.html) or [Mamba](https://mamba.readthedocs.io/). Works on **Linux / macOS / WSL2**.

```bash
# Download the environment file and create the environment
conda env create -f https://raw.githubusercontent.com/zhaoqing-wang/CellJanus/main/environment.yml

# Activate
conda activate celljanus

# Verify — all tools should show ✔ Found
celljanus check
```

Or if you already cloned the repository:

```bash
conda env create -f environment.yml
conda activate celljanus
celljanus check
```

<details>
<summary><b>Alternative installation methods</b></summary>

#### pip only (Python package without external tools)

CellJanus itself is a pure-Python orchestrator. If you already have fastp, Bowtie2, etc. on your PATH:

```bash
pip install celljanus
celljanus --version
```

#### From source (development)

```bash
git clone https://github.com/zhaoqing-wang/CellJanus.git
cd CellJanus
pip install -e ".[dev]"
```

#### Docker

```bash
docker build -t celljanus .
docker run --rm celljanus celljanus check
```

</details>

All tools should show **✔ Found**. STAR is optional (for future RNA-seq alignment support).

---

## Quick Start

The repository includes test data and pre-built reference databases — run the full pipeline immediately with **no downloads required**.

```bash
conda activate celljanus
cd celljanus

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

#### Example Output

| Pipeline Dashboard |
|:--:|
| ![Pipeline Dashboard](docs/pipeline_dashboard.png) |
| *Summarises QC, alignment and classification metrics in a single view.* |

| Abundance Bar Chart | Abundance Donut Chart | Abundance Heatmap |
|:--:|:--:|:--:|
| ![Bar](docs/abundance_bar.png) | ![Pie](docs/abundance_pie.png) | ![Heatmap](docs/abundance_heatmap.png) |
| Top species ranked by read count. | Relative proportion of each species. | Log₁₀-scaled heatmap of species abundance. |

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

## scRNA-seq Mode

CellJanus supports per-cell microbial abundance tracking for single-cell RNA-seq data from 10x Genomics, Parse Biosciences, and other platforms.

**Key capabilities**:
- Extract cell barcodes (CB) and UMIs from 10x Genomics-style FASTQ headers
- Track microbial abundances per cell barcode with UMI deduplication
- Generate comprehensive exports: raw counts, normalized values, summaries
- Optimized for large datasets (10,000+ cells) with smart sampling
- WSL2 I/O optimization with cross-filesystem detection

### Quick Start (scRNA-seq)

```bash
celljanus scrnaseq \
    --read1 sample_S1_L001_R1_001.fastq.gz \
    --read2 sample_S1_L001_R2_001.fastq.gz \
    --kraken2-db ./refs/standard_8 \
    --output-dir scrna_results \
    --barcode-mode 10x \
    --threads 8
```

### Supported Barcode Formats

| Platform | Mode | Barcode Location |
|----------|------|------------------|
| 10x Genomics 3'/5' | `10x` | Read header: `CB:Z:BARCODE UB:Z:UMI` |
| Parse Biosciences | `parse` | Read name (colon-separated) |
| Custom | `auto` | Auto-detect from multiple patterns |

### scRNA-seq Output

```
scrna_results/
├── classification/
│   ├── kraken2_report.txt
│   └── kraken2_output.txt
├── tables/                           # Comprehensive data exports
│   ├── cell_species_counts.csv       # Raw counts (cells × species)
│   ├── cell_species_normalized.csv   # CPM-normalized (for Seurat/Scanpy)
│   ├── cell_species_long.csv         # Long format (for ggplot/seaborn)
│   ├── species_summary.csv           # Per-species statistics
│   ├── cell_summary.csv              # Per-cell summary (diversity, etc.)
│   └── pipeline_summary.csv          # Pipeline metrics
└── plots/
    ├── cell_species_heatmap.*        # Per-cell abundance heatmap
    ├── cell_microbe_summary.*        # Distribution summaries (3 panels)
    └── cell_bacteria_dotplot.*       # Cell–bacteria association
```

#### CSV Output Details

| File | Description | Use Case |
|------|-------------|----------|
| `cell_species_counts.csv` | Raw read/UMI counts per cell | Downstream analysis, integration |
| `cell_species_normalized.csv` | CPM-like normalized values | Seurat AddModuleScore, Scanpy |
| `cell_species_long.csv` | Tidy format (cell, species, count, fraction) | ggplot2, seaborn visualisation |
| `species_summary.csv` | Species: total reads, n_cells, prevalence, mean/median | Species-level statistics |
| `cell_summary.csv` | Cell: total reads, n_species, Shannon diversity | Cell-level QC metrics |
| `pipeline_summary.csv` | Overall pipeline statistics | Quality reports |

### Python API (scRNA-seq)

```python
from pathlib import Path
from celljanus.scrnaseq import BarcodeConfig, run_scrnaseq_classification
from celljanus.config import CellJanusConfig

cfg = CellJanusConfig(
    output_dir=Path("./scrna_results"),
    kraken2_db=Path("./refs/standard_8"),
    threads=8,
)

barcode_cfg = BarcodeConfig(
    mode="10x",
    min_reads_per_cell=10,
)

result = run_scrnaseq_classification(
    Path("sample_R1.fastq.gz"),
    Path("./refs/standard_8"),
    Path("./scrna_results"),
    read2=Path("sample_R2.fastq.gz"),
    barcode_cfg=barcode_cfg,
    cfg=cfg,
)

# Access per-cell data
abundance = result["abundance"]

# Different matrix formats
raw_matrix = abundance.to_matrix()           # Raw counts (cells × species)
normalized = abundance.to_normalized_matrix() # CPM-like normalized
long_format = abundance.to_long_format()     # Tidy format for plotting

# Summary statistics
species_stats = abundance.to_species_summary()  # Per-species statistics
cell_stats = abundance.to_cell_summary()        # Per-cell metrics (diversity, etc.)
summary = result["summary"]                     # Overall pipeline statistics
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
| `celljanus scrnaseq` | **NEW** scRNA-seq mode with per-cell barcode tracking |
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

## WSL2 Performance Tips

When running CellJanus on WSL2 (Windows Subsystem for Linux 2), file I/O can be a significant bottleneck when accessing Windows filesystem paths (e.g., `/mnt/c/`, `/mnt/d/`).

### Recommended Workflow

1. **Store data on Linux filesystem**: Copy your FASTQ files to a native Linux path (e.g., `/home/user/data/`) for 10-50× faster I/O performance.

   ```bash
   # Create workspace on Linux filesystem
   mkdir -p ~/celljanus_work
   
   # Copy data from Windows (one-time)
   cp /mnt/c/Data/sample_R1.fastq.gz ~/celljanus_work/
   cp /mnt/c/Data/sample_R2.fastq.gz ~/celljanus_work/
   
   # Run pipeline on Linux filesystem
   celljanus run \
       -1 ~/celljanus_work/sample_R1.fastq.gz \
       -2 ~/celljanus_work/sample_R2.fastq.gz \
       -x ~/refs/bowtie2_index/hg38 \
       -d ~/refs/kraken2_db \
       -o ~/celljanus_work/results
   ```

2. **CellJanus auto-detection**: The `scrnaseq` command automatically detects WSL2 and warns about cross-filesystem paths:

   ```
   ⚠️  Performance Warning: 2 path(s) are on Windows filesystem (/mnt/...).
   For better I/O performance, copy data to a native Linux path.
   ```

3. **Memory-mapped databases**: Store Kraken2 databases on Linux filesystem to enable efficient memory mapping.

---

## STAR Installation (WSL2/Linux)

STAR (Spliced Transcripts Alignment to a Reference) is optional but useful for RNA-seq alignment. Here's how to install it on WSL2 or Linux:

### Method 1: Conda (Recommended)

```bash
# Install via Bioconda (easiest)
conda install -c bioconda star=2.7.11b

# Verify installation
STAR --version
```

### Method 2: Pre-built Binary

```bash
# Download pre-compiled binary
wget https://github.com/alexdobin/STAR/releases/download/2.7.11b/STAR_2.7.11b.zip
unzip STAR_2.7.11b.zip

# Add to PATH
sudo cp STAR_2.7.11b/Linux_x86_64_static/STAR /usr/local/bin/
STAR --version
```

### Method 3: Compile from Source

```bash
# Clone and build (requires g++ and make)
git clone https://github.com/alexdobin/STAR.git
cd STAR/source

# Compile (adjust threads as needed)
make -j4 STAR

# Install
sudo cp STAR /usr/local/bin/
```

### STAR Hardware Requirements

| Genome | RAM Required | Note |
|--------|--------------|------|
| Human (hg38) | 32+ GB | Recommended: 64 GB for large datasets |
| Mouse (mm10) | 32+ GB | Similar to human |
| Bacteria | 4-8 GB | Much smaller genomes |

> **Note**: CellJanus currently uses Bowtie2 for host alignment by default. STAR integration for splice-aware RNA-seq alignment is planned for future versions.

---

## Citation

```
Wang Z (2026). CellJanus: A Dual-Perspective Tool for Deconvolving Host
Single-Cell and Microbial Transcriptomes. Python package.
https://github.com/zhaoqing-wang/CellJanus
```

## License

[MIT](LICENSE)
