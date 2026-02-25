# Changelog

All notable changes to CellJanus will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.5] — 2026-02-25

### Added

- **scRNA-seq support**: New `celljanus scrnaseq` command for per-cell microbial abundance tracking.
  - Parses cell barcodes (CB) and UMIs (UB) from 10x Genomics, Parse Biosciences, and custom formats.
  - Generates cell × species abundance matrix for integration with Seurat/Scanpy.
  - **6 comprehensive CSV outputs**:
    - `cell_species_counts.csv` — Raw read/UMI counts (cells × species matrix).
    - `cell_species_normalized.csv` — CPM-like normalized values for Seurat/Scanpy.
    - `cell_species_long.csv` — Tidy format for ggplot2/seaborn.
    - `species_summary.csv` — Per-species stats (total reads, n_cells, prevalence, mean/median).
    - `cell_summary.csv` — Per-cell metrics (total reads, n_species, Shannon diversity).
    - `pipeline_summary.csv` — Overall pipeline metrics.
- **scRNA-seq visualizations** (optimized for 10,000+ cells):
  - `cell_species_heatmap`: Per-cell microbial abundance heatmap with smart sampling (stratified/random/top).
  - `cell_microbe_summary`: 3-panel summary with log-scale support for large dynamic range.
  - `cell_bacteria_dotplot`: Dot plot with efficient filtering and dynamic figure sizing.
- **Large dataset optimizations**:
  - Smart cell sampling for visualizations (stratified sampling by read depth).
  - Automatic log-scale axes when data range exceeds 100×.
  - Memory-efficient figure sizing limits.
- **WSL2 I/O optimization**:
  - Auto-detection of WSL2 environment.
  - Warnings for cross-filesystem paths (`/mnt/c/`, `/mnt/d/`) with performance recommendations.
  - `wsl2_io_warning()` and `recommend_native_path()` helpers in `scrnaseq` module.
- **CellMicrobialAbundance class methods**:
  - `to_matrix()` — Raw counts matrix.
  - `to_normalized_matrix()` — CPM-like normalized (scale=10000).
  - `to_long_format()` — Tidy format with count and fraction columns.
  - `to_species_summary()` — Species-level statistics.
  - `to_cell_summary()` — Cell-level metrics with Shannon diversity.
- **New module**: `celljanus/scrnaseq.py` with barcode extraction, `BarcodeConfig`, `CellMicrobialAbundance` classes.
- **STAR installation guide** in README: conda, pre-built binary, and source compilation methods for WSL2/Linux.
- **scRNA-seq test data generator**: `generate_scrnaseq_fastq()` function in `tests/generate_test_data.py`.
- **Test suite expanded**: 26 tests covering scRNA-seq barcode extraction, cell abundance tracking, CLI commands.

### Changed

- README restructured with new sections: scRNA-seq Mode, WSL2 Performance Tips, STAR Installation.
- README updated with comprehensive CSV output documentation.
- Pipeline description updated to show both bulk RNA-seq and scRNA-seq workflows.
- CLI Reference updated with `celljanus scrnaseq` command.
- Python API documentation expanded with full scRNA-seq export examples.
- Version bumped to 0.1.5.

### Fixed

- `is_cross_filesystem_path()` now handles cross-platform execution correctly (returns False on native Windows, works properly in WSL2).
- Test suite platform compatibility fixed for WSL2 detection tests.
- Unicode subscript characters (`log₁₀`) replaced with LaTeX mathtext (`$\log_{10}$`) for font compatibility.

### Tested

- **Full WSL2 validation** (Ubuntu 24.04 on WSL2):
  - All 26 pytest tests passed.
  - Full pipeline test completed successfully with fastp, bowtie2, kraken2, bracken.
  - scRNA-seq module WSL2 detection verified.
  - Cross-filesystem path detection working correctly.
  - Tool versions tested: STAR 2.7.11b, fastp 0.23.4, bowtie2 2.5.2, kraken2 2.1.3, samtools 1.19.2.

## [0.1.4] — 2026-02-24

### Added

- **Published to PyPI**: `pip install celljanus` now works directly.
- **GitHub Actions CI/CD**: automatic PyPI publishing on new GitHub Release (`.github/workflows/publish.yml`).
- **PyPI badge** and **Downloads badge** in README.

### Changed

- README Installation section restructured: `pip install` as primary method; conda, source, and Docker as alternatives.
- Citation updated to version 0.1.4.

## [0.1.3] — 2026-02-22

### Added

- **Dual-format plot output**: every plot is now saved as both PNG (raster, 300 dpi) and PDF (vector) for publication use.
- **Result tables** in `06_tables/`:
  - `pipeline_summary.csv` — QC, alignment, and classification metrics (one row per metric).
  - `species_abundance.csv` — clean species table with taxonomy ID and fraction percentage.
  - `output_manifest.csv` — complete file inventory with byte sizes.
- Pipeline summary table in CLI now shows file sizes and plot/table counts.

### Changed

- Publication-quality plot styling: global `rcParams` for consistent fonts, spines, and layout.
- Bar chart: value labels with read counts and percentages displayed on bars.
- Donut chart: moved species labels to a side legend; total read count shown in centre.
- Heatmap: wider figure, rotated y-labels, improved colour bar.
- Dashboard: gridspec layout with colour-coded panels (blue/green/red) and structured key-value formatting.
- Version bumped to 0.1.3.

## [0.1.2] — 2026-02-23

### Fixed

- Bowtie2 alignment: replaced `--no-unal` with `--un-conc-gz` (PE) / `--un-gz` (SE) for direct unmapped FASTQ output — standard host-subtraction approach that correctly captures all non-host reads.
- `align_to_host()` now returns `(bam, unmapped_r1, unmapped_r2)` tuple; pipeline no longer needs a separate samtools extraction step for unmapped reads.
- Bracken: graceful handling of empty Kraken2 reports; falls back to Kraken2 report for summary stats.
- Dashboard plot: fixed `ValueError` in `plot_qc_dashboard()` when computing retained-read percentage.
- Default `align_extra_args` changed from `"--very-sensitive --no-unal"` to `""` to avoid silently discarding unmapped reads.

### Added

- Pre-built Bowtie2 index and Kraken2+Bracken database in `testdata/refs/` — full pipeline runs with zero setup.
- `testdata/build_test_refs.py`: script to regenerate test reference databases from NCBI taxonomy.
- Full end-to-end pipeline test validated: QC (90%), Alignment (66%), Classification (33%, 3 species), Visualisation (4 plots).

### Changed

- README fully rewritten: conda virtual environment as primary installation, full pipeline "Try It Now" section with real test results, accurate output file listing.
- Pipeline steps 2 + 3 merged — Bowtie2 handles unmapped extraction directly via `--un-conc-gz`.
- Version bumped to 0.1.2.

## [0.1.1] — 2026-02-22

### Added

- Shipped real-format FASTQ test data in `testdata/` (1,000 PE + SE reads with Illumina headers).
  - 600 human reads (HBB, ACTB), 300 microbial (S. aureus, E. coli, K. pneumoniae 16S), 100 low-quality.
  - Enables immediate `celljanus qc` testing without any database downloads.
- New test data generator: `tests/generate_testdata.py` with realistic Illumina-style read headers.
- Docker image now runs full test suite + QC smoke test during build.
- Docker image includes bundled `testdata/` for instant demos.

### Changed

- README restructured: Docker-first installation, "Try It Now" section with bundled test data.
- README simplified: removed verbose platform matrices, added step-by-step Docker verification.
- Version bumped to 0.1.1.

## [0.1.0] — 2026-02-22

### Added

- Initial release of CellJanus.
- Full pipeline: QC (fastp) → Host alignment (Bowtie2) → Extract unmapped (samtools) → Classify (Kraken2 + Bracken) → Visualize (matplotlib / seaborn).
- Click-based CLI with subcommands: `run`, `qc`, `align`, `extract`, `classify`, `visualize`, `download`, `check`.
- Reference data downloader for hg38 genome and Kraken2 databases.
- Publication-quality plots: abundance bar charts, donut plots, heatmaps, QC dashboards.
- Memory-efficient streaming I/O for large FASTQ files.
- Docker support with Miniforge + bioconda for reproducible deployments.
- Comprehensive test suite (18 tests).
- Platform support: Linux, macOS, Windows (WSL2).
