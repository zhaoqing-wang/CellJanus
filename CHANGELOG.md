# Changelog

All notable changes to CellJanus will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

> **Note**: When updating the package version, synchronize the version info in [pyproject.toml](pyproject.toml) (line 7), [Dockerfile](Dockerfile) (lines 6, 10, 13, 28), [docker-compose.yml](docker-compose.yml) (line 13), and [celljanus/__init__.py](celljanus/__init__.py) (line 13).

---

## [0.2.2] — 2026-03-11

### Fixed

- **scRNA-seq output consistency**: All six output CSV files (`cell_species_counts.csv`, `cell_species_normalized.csv`, `cell_species_long.csv`, `species_summary.csv`, `cell_summary.csv`, `pipeline_summary.csv`) and the CLI summary table now reflect the **same filtered cell set** determined by `--min-reads`. Previously, `pipeline_summary.csv` and `species_summary.csv` could report pre-filter counts while the count matrix was correctly filtered, causing metric mismatches.
- **Bulk Kraken2 summary totals**: `Classified reads` and `Unclassified` in bulk mode are now derived from the Kraken2 root and unclassified rows instead of assuming the first report row is the total.

### Added

- **Early cell filtering (`filter_cells`)**: After Kraken2 classification, cells below the `--min-reads` threshold are removed from the internal data structure *before* any CSV export. This reduces peak memory and I/O for large datasets (e.g., 2 M raw barcodes filtered to 50 K cells avoids building a 2 M × 3 K matrix in RAM).
- **`input_reads` metric**: `pipeline_summary.csv` now records the total number of reads in the input FASTQ.
- **Raw pre-filter traceability**: `pipeline_summary.csv` retains `total_cells_raw`, `species_detected_raw`, `total_microbial_reads_raw`, and `cells_filtered_out` alongside the filtered headline metrics.
- **Recommended `--min-reads` values**: README §3.5 now includes a guidance table for choosing `--min-reads` (1 for testdata, 5–10 with whitelist, 50+ without whitelist).

### Changed

- **Unified INFO logging**: All CLI and pipeline output (banners, summary tables, status messages, and external tool output from fastp/kraken2/bracken/samtools) now goes through the Python `logging` module with `[HH:MM:SS] INFO` prefix. Previously some output used `rich.console.Console.print()` (no prefix) and external tools wrote directly to stderr, causing inconsistent formatting. Added `log_renderable()` helper for Rich objects and changed `run_cmd()` to capture subprocess output line-by-line and route it through the logger.
- **CLI summary labels clarified**: "Cells passing --min-reads filter", "Species detected (filtered)", "Total microbial reads (filtered)", and "Mean reads / cell (filtered)" make it unambiguous that numbers refer to the post-filter cell set.
- **Dashboard plot label**: "Total Cells" renamed to "Cells (passing --min-reads)" in the scRNA-seq dashboard key-metrics panel.
- **README output file table** (§3.5): `cell_species_counts.csv` now described as "Count matrix (cells × species). Only cells with ≥ `--min-reads` microbial reads are included." instead of the misleading "Raw counts". Full `--min-reads` filtering explanation and metric definitions added.
- **README test results updated**: scRNA-seq testdata and real-reference result tables updated to match v0.3.0 output with host-taxa filtering metrics.
- **CHANGELOG consolidated**: Combined granular patch entries (0.1.0–0.2.1) into clearer milestone summaries to improve readability.

### Tested (WSL2 Ubuntu 24.04)

- **`celljanus check`**: All 10 tools detected, all output uniformly through `[HH:MM:SS] INFO`.
- **Bulk pipeline (testdata + kraken2_testdb)**: 1,000 PE reads → 950 QC-passed → 31.58% host alignment → 5 species detected (~3 s).
- **scRNA-seq (testdata + kraken2_testdb)**: 15,000 reads → 300 cells × 7 species, 2,395 microbial reads, 8.0 mean reads/cell (~2 s).
- **scRNA-seq (testdata + standard_8)**: 15,000 reads → 299 cells × 32 species, 1,543 retained reads after host filtering (~3 s).
- All pipeline_summary, species_summary, and matrix dimensions verified consistent.

---

## [0.2.1] — 2026-03-10

### Changed

- **Sequence-based barcode extraction**: When R1 headers lack `CB:Z:`/`UB:Z:` tags (standard for raw 10x Chromium FASTQ), barcodes are now extracted from the R1 sequence (first 16 bp = cell barcode, next 12 bp = UMI for v3 chemistry). Reads with `N` in the barcode are filtered out.

---

## [0.2.0] — 2026-03-10

### Changed

- **Memory-efficient scRNA-seq pipeline**: Rewrote barcode extraction to stream R1 FASTQ to a temp file instead of in-memory dict. Eliminates ~30 GB peak memory on large datasets (~170 M reads). New peak: ~1–2 GB.
- **Streaming Kraken2 join**: Position-based lockstep iteration replaces O(n) read-ID lookup, using O(cells × species) memory.
- **R2-only classification**: When R2 is provided, Kraken2 classifies R2 alone (cDNA insert). R1 is 28 bp for 10x v3 (shorter than k=35, zero k-mers), so R2-only gives identical results and halves I/O.
- **Default host/non-informative taxon removal**: `scrnaseq` now removes `Homo sapiens`, `cellular organisms`, `root`, and `other sequences` before per-cell aggregation. `--keep-host-taxa` restores unfiltered output.
- **Step duration logging**: Explicit timing logs for all pipeline steps and table export.
- **Progress logging**: Step 1 logs progress every 5 M reads.

---

## [0.1.9] — 2026-03-09

### Changed

- README test result tables updated to match FR paired-end and genomic DNA template fixes.
- Fixed scRNA-seq real-reference command paths in README.
- testdata/README.md updated for genomic DNA templates and dual header formats.

---

## [0.1.8] — 2026-03-08

### Fixed

- **Bowtie2 index download URL**: S3 prefix `bt2/` → `bt/` (404 fix).
- **Paired-end read orientation**: Test data now produces proper FR reads (was FF, causing all pairs to be discordant).
- **Host templates**: Replaced exon-junction cDNA templates with contiguous genomic DNA (TP53/MSH2) for correct Bowtie2 alignment.
- **Test datasets regenerated**: Validated RefSeq genome fragments for 7 species.

---

## [0.1.7] — 2026-03-08

### Changed

- Automated PyPI release on GitHub Release.
- Docker, pyproject.toml, and README presentation improvements.

---

## [0.1.6] — 2026-02-25

### Added

- **scRNA-seq Dashboard**: 4-panel summary (abundance bar, read depth histogram, species prevalence, key metrics).
- **scRNA-seq Abundance Pie**: Donut chart of community composition.
- **300-cell scRNA-seq test data**: 15,000 reads across 5 abundance profiles.

### Changed

- Visualizations streamlined to 3 plots (dashboard, pie, 3-panel summary); heatmap/dotplot removed (don't scale for 5,000+ cells).
- README unified: Quick Test → Results → Real Data structure for both modes.

---

## [0.1.5] — 2026-02-25

### Added

- **scRNA-seq support** (`celljanus scrnaseq`): Per-cell microbial abundance tracking with 10x Genomics, Parse Biosciences, and custom barcode extraction.
- **6 CSV outputs**: counts, normalized, long-format, species summary, cell summary, pipeline summary.
- **CellMicrobialAbundance class**: `to_matrix()`, `to_normalized_matrix()`, `to_long_format()`, `to_species_summary()`, `to_cell_summary()`.
- **WSL2 I/O optimization**: Auto-detection, cross-filesystem warnings, native-path recommendations.
- **26 tests** covering scRNA-seq barcode extraction, cell abundance tracking, CLI commands.

---

## [0.1.4] — 2026-02-24

### Added

- Published to PyPI (`pip install celljanus`).
- GitHub Actions CI/CD for automatic PyPI publishing.

---

## [0.1.3] — 2026-02-22

### Added

- Dual-format plot output (PNG 300 dpi + PDF vector).
- Result tables in `06_tables/`: `pipeline_summary.csv`, `species_abundance.csv`, `output_manifest.csv`.

### Changed

- Publication-quality plot styling with improved bar chart, donut, heatmap, and dashboard layouts.

---

## [0.1.2] — 2026-02-23

### Fixed

- Bowtie2 alignment: `--un-conc-gz`/`--un-gz` for proper unmapped FASTQ output.
- Bracken: graceful handling of empty reports.
- Dashboard `ValueError` fix.

### Added

- Pre-built test reference databases in `testdata/refs/`.
- Full end-to-end pipeline test (QC 90%, Alignment 66%, 3 species).

---

## [0.1.1] — 2026-02-22

### Added

- Real-format FASTQ test data (1,000 PE + SE reads).
- Docker image with bundled `testdata/` and build-time test suite.

---

## [0.1.0] — 2026-02-22

### Added

- Initial release: QC (fastp) → Host alignment (Bowtie2) → Extract unmapped (samtools) → Classify (Kraken2 + Bracken) → Visualize (matplotlib/seaborn).
- CLI: `run`, `qc`, `align`, `extract`, `classify`, `visualize`, `download`, `check`.
- Docker support, 18 tests, Linux/macOS/WSL2.
