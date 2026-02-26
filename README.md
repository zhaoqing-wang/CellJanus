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

**Bulk RNA-seq:** `FASTQ → fastp (QC) → Bowtie2 (host) → unmapped reads → Kraken2+Bracken → plots + CSV`

**scRNA-seq:** `10x FASTQ → Extract CB+UMI → Kraken2 → Per-cell abundance → Cell×Species matrix`

---

## Installation

```bash
conda env create -f https://raw.githubusercontent.com/zhaoqing-wang/CellJanus/main/environment.yml
conda activate celljanus
celljanus check   # All tools should show ✔ Found
```

> **Requirements**: [Conda](https://docs.conda.io/en/latest/miniconda.html) or [Mamba](https://mamba.readthedocs.io/). Works on **Linux / macOS / WSL2**.

<details>
<summary><b>Alternative installation methods</b></summary>

```bash
# pip only (requires fastp, bowtie2, samtools, kraken2, bracken on PATH)
pip install celljanus

# Development install
git clone https://github.com/zhaoqing-wang/CellJanus.git
cd CellJanus && pip install -e ".[dev]"

# Docker
docker build -t celljanus . && docker run --rm celljanus celljanus check
```

</details>

---

## Bulk RNA-seq Mode

Full pipeline: QC → Host alignment → Microbial classification → Visualization.

### Quick Test (no downloads required)

```bash
celljanus run \
    --read1 testdata/reads_R1.fastq.gz \
    --read2 testdata/reads_R2.fastq.gz \
    --host-index testdata/refs/host_genome/host \
    --kraken2-db testdata/refs/kraken2_testdb \
    --output-dir test_results/bulk
```

| Pipeline Dashboard |
|:--:|
| ![Dashboard](docs/pipeline_dashboard.png) |

| Abundance Bar | Abundance Pie | Abundance Heatmap |
|:--:|:--:|:--:|
| ![Bar](docs/abundance_bar.png) | ![Pie](docs/abundance_pie.png) | ![Heatmap](docs/abundance_heatmap.png) |

<details>
<summary><b>Test results (~4 seconds)</b></summary>

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

</details>

### Real Data

```bash
# 1. Download references (~13 GB total)
celljanus download hg38 -o ./refs           # Human genome + Bowtie2 index (~5 GB)
celljanus download kraken2 -o ./refs        # Kraken2 standard DB (~8 GB)

# 2. Run pipeline
celljanus run \
    -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
    -x ./refs/bowtie2_index/GRCh38_noalt_as \
    -d ./refs/standard_8 \
    -o ./results --threads 8
```

<details>
<summary><b>Advanced options (custom QC, skip steps, single-end)</b></summary>

```bash
# Custom QC and classification parameters
celljanus run \
    -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
    -x ./refs/bowtie2_index/GRCh38_noalt_as \
    -d ./refs/standard_8 -o ./results --threads 8 \
    --min-quality 20 --min-length 50 \
    --confidence 0.1 --bracken-level G \
    --plot-format pdf --max-memory 16

# Skip steps for partial re-runs
celljanus run ... --skip-qc --skip-visualize

# Single-end mode
celljanus run -1 sample_SE.fastq.gz -x ./refs/... -d ./refs/... -o ./results
```

</details>

---

## scRNA-seq Mode

Per-cell microbial abundance tracking. Only a Kraken2 database is needed (no host alignment).

### Quick Test (no downloads required)

```bash
celljanus scrnaseq \
    --read1 testdata/scrnaseq/scrna_R1.fastq.gz \
    --read2 testdata/scrnaseq/scrna_R2.fastq.gz \
    --kraken2-db testdata/refs/kraken2_testdb \
    --output-dir test_results/scrnaseq \
    --barcode-mode 10x --min-reads 1
```

| scRNA-seq Dashboard |
|:--:|
| ![scRNA-seq Dashboard](docs/scrnaseq_dashboard.png) |

| Species Abundance Pie | Microbial Summary |
|:--:|:--:|
| ![Abundance Pie](docs/scrnaseq_abundance_pie.png) | ![Summary](docs/cell_microbe_summary.png) |

<details>
<summary><b>Test results (~2 seconds)</b></summary>

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

</details>

### Real Data

```bash
# 1. Download Kraken2 database (~8 GB)
celljanus download kraken2 -o ./refs

# 2. Run scRNA-seq pipeline (10x Genomics)
celljanus scrnaseq \
    --read1 sample_R1.fastq.gz \             # R1: cell barcodes + UMI
    --read2 sample_R2.fastq.gz \             # R2: cDNA sequences
    --kraken2-db ./refs/standard_8 \
    --output-dir scrna_results \
    --barcode-mode 10x --threads 8

# With barcode whitelist filtering (recommended for real data)
celljanus scrnaseq \
    --read1 sample_R1.fastq.gz --read2 sample_R2.fastq.gz \
    --kraken2-db ./refs/standard_8 \
    --output-dir scrna_results \
    --barcode-mode 10x \
    --whitelist 3M-february-2018.txt.gz \    # 10x v3 barcode whitelist
    --min-reads 5 --threads 8
```

<details>
<summary><b>Parse Biosciences / auto-detect mode</b></summary>

```bash
# Parse Biosciences (barcode in read name)
celljanus scrnaseq \
    --read1 parse_R1.fastq.gz --read2 parse_R2.fastq.gz \
    --kraken2-db ./refs/standard_8 \
    --output-dir scrna_parse_results \
    --barcode-mode parse --min-reads 3 --threads 8

# Auto-detect barcode format
celljanus scrnaseq \
    --read1 sample_R1.fastq.gz --read2 sample_R2.fastq.gz \
    --kraken2-db ./refs/standard_8 \
    --output-dir scrna_auto_results \
    --barcode-mode auto
```

| Platform | Mode | Barcode Location |
|----------|------|------------------|
| 10x Genomics | `10x` | Header: `CB:Z:BARCODE UB:Z:UMI` |
| Parse Biosciences | `parse` | Read name (colon-separated) |
| Custom | `auto` | Auto-detect |

</details>

### Downstream Integration (Seurat / Scanpy)

```r
# R / Seurat
library(Seurat)
sc <- readRDS("seurat_object.rds")
microbe_mat <- read.csv("scrna_results/tables/cell_species_normalized.csv",
                        row.names = 1, check.names = FALSE)
sc[["Microbe"]] <- CreateAssayObject(counts = t(as.matrix(microbe_mat)))
```

```python
# Python / Scanpy
import scanpy as sc, pandas as pd
adata = sc.read_h5ad("adata.h5ad")
microbe_df = pd.read_csv("scrna_results/tables/cell_species_normalized.csv", index_col=0)
common = adata.obs_names.intersection(microbe_df.index)
adata[common].obsm["X_microbe"] = microbe_df.loc[common].values
```

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
| `celljanus run` | Full bulk pipeline: QC → Align → Classify → Visualize |
| `celljanus scrnaseq` | scRNA-seq per-cell microbial classification |
| `celljanus qc` | Quality control only (fastp) |
| `celljanus align` | Host alignment only (Bowtie2) |
| `celljanus extract` | Extract unmapped reads from BAM |
| `celljanus classify` | Taxonomic classification (Kraken2 + Bracken) |
| `celljanus visualize` | Generate abundance plots |
| `celljanus download` | Download reference databases (hg38, kraken2) |
| `celljanus check` | Verify external tool installation |

Run `celljanus <command> --help` for full option details.

<details>
<summary><b>Key options reference</b></summary>

| Option | Default | Description |
|--------|---------|-------------|
| `-1, --read1` | *required* | R1 FASTQ file |
| `-2, --read2` | — | R2 FASTQ (paired-end) |
| `-x, --host-index` | *required* | Bowtie2 index prefix (bulk only) |
| `-d, --kraken2-db` | *required* | Kraken2 database path |
| `-o, --output-dir` | `celljanus_output` | Output directory |
| `-t, --threads` | auto | Worker threads |
| `--min-quality` | 15 | Phred quality threshold (bulk QC) |
| `--confidence` | 0.05 | Kraken2 confidence threshold |
| `--bracken-level` | S | Taxonomic level (D/P/C/O/F/G/S) |
| `--barcode-mode` | 10x | Barcode format: 10x / parse / auto (scRNA-seq) |
| `-w, --whitelist` | — | Cell barcode whitelist (scRNA-seq) |
| `--min-reads` | 1 | Min reads per cell (scRNA-seq) |
| `--skip-qc` | — | Skip QC step (bulk) |
| `--plot-format` | png | Output format: png / pdf / svg |

</details>

---

## Python API

<details>
<summary><b>Bulk pipeline</b></summary>

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

</details>

<details>
<summary><b>scRNA-seq pipeline</b></summary>

```python
from pathlib import Path
from celljanus.config import CellJanusConfig
from celljanus.scrnaseq import BarcodeConfig, run_scrnaseq_classification

cfg = CellJanusConfig(
    output_dir=Path("./scrna_results"),
    kraken2_db=Path("./refs/standard_8"),
    threads=8,
)
barcode_cfg = BarcodeConfig(
    mode="10x",               # 10x / parse / auto
    min_reads_per_cell=5,     # Filter low-quality cells
)

result = run_scrnaseq_classification(
    Path("sample_R1.fastq.gz"),
    Path("./refs/standard_8"),
    Path("./scrna_results"),
    read2=Path("sample_R2.fastq.gz"),
    barcode_cfg=barcode_cfg,
    cfg=cfg,
)

# Access per-cell abundance data
abundance = result["abundance"]
abundance.to_matrix()             # Raw counts (cells × species)
abundance.to_normalized_matrix()  # CPM-normalized (Seurat/Scanpy)
abundance.to_long_format()        # Tidy format (ggplot2/seaborn)
abundance.to_species_summary()    # Per-species statistics
abundance.to_cell_summary()       # Per-cell diversity metrics
```

</details>

---

## Output Structure

<details>
<summary><b>Bulk RNA-seq output</b></summary>

```
results/
├── 01_qc/                        # fastp QC: *_qc.fastq.gz, *.json, *.html
├── 02_alignment/                 # Bowtie2: BAM, unmapped FASTQs, stats
├── 04_classification/            # Kraken2 report + Bracken species table
├── 05_visualisation/plots/       # 4 PNG + 4 PDF (bar, pie, heatmap, dashboard)
├── 06_tables/                    # species_abundance.csv, pipeline_summary.csv
└── celljanus.log
```

</details>

<details>
<summary><b>scRNA-seq output</b></summary>

```
scrna_results/
├── classification/               # Kraken2 report + per-read output
├── tables/                       # 6 CSVs: counts, normalized, long, species/cell summary
├── visualisation/plots/          # 3 PNG + 3 PDF (dashboard, pie, 3-panel summary)
└── celljanus.log
```

</details>

---

## Advanced

<details>
<summary><b>WSL2 optimization</b></summary>

For best performance on WSL2, store data on the Linux filesystem (10–50× faster I/O):

```bash
mkdir -p ~/celljanus_work
cp /mnt/c/Data/sample*.fastq.gz ~/celljanus_work/
celljanus run -1 ~/celljanus_work/sample_R1.fastq.gz ...
```

CellJanus auto-detects WSL2 and warns about slow cross-filesystem paths.

</details>

---

## Citation

```
Wang Z (2026). CellJanus: Dual-Perspective Deconvolution of Host and Microbial
Transcriptomes from FASTQ Data. https://github.com/zhaoqing-wang/CellJanus
```

## License

[MIT](LICENSE)

## Contact

**Author:** Zhaoqing Wang ([ORCID](https://orcid.org/0000-0001-8348-7245)) | **Email:** <zhaoqingwang@mail.sdu.edu.cn> | **Issues:** [CellJanus Issues](https://github.com/zhaoqing-wang/CellJanus/issues)
