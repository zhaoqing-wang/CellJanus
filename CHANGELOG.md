# Changelog

All notable changes to CellJanus will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

> **Note**: When updating the package version, synchronize the version info in [pyproject.toml](pyproject.toml) (line 7), [Dockerfile](Dockerfile) (lines 6, 10, 13, 28), [docker-compose.yml](docker-compose.yml) (line 13), and [celljanus/__init__.py](celljanus/__init__.py) (line 13).

---

## [0.2.4] — 2026-03-12

### Summary

- **Unified CLI naming**: Renamed `celljanus run` to `celljanus bulk` so both analysis modes have explicit, parallel commands (`celljanus bulk` / `celljanus scrnaseq`).
- **Mode-labeled sub-commands**: Bulk-only steps (`align`, `extract`) now show `[Bulk]` in `--help` output; shared and scRNA-seq commands are similarly clarified.
- **README & Dockerfile updated**: All examples and CLI reference table reflect the new `celljanus bulk` command.

---

## [0.2.3] — 2026-03-12

### Summary

- **scRNA-seq run completion output improved**: now prints a final green completion panel with elapsed time and output directory, aligned with bulk mode style.
- **System identification remains explicit in scRNA-seq mode**: startup output includes OS/CPU/RAM/Threads/MaxMem plus input/output panel for easier run auditing.
- **`cell_microbe_summary` layout optimized**: adjusted figure size/top spacing to prevent main title overlap with subplot titles in PNG/PDF exports.
- **`scrnaseq_dashboard` key-metrics readability improved**: long first metric label was shortened and the `--min-reads` note moved to a subtitle to avoid text/value overlap.
---

## [0.2.2] — 2026-03-12

### Summary

- **scRNA-seq metrics are now consistent**: all exports and CLI summary use the same post-`--min-reads` filtered cell set.
- **Bulk classification reporting is more robust**: classified/unclassified totals now come from correct Kraken2 rows.
- **Bracken failures no longer crash bulk runs**: when species-level estimation is unavailable, CellJanus falls back to Kraken2 genus/family summaries.
- **Host signal removal is optimized and default-enabled (scRNA-seq)**: host/non-informative taxa (e.g., `Homo sapiens`, `root`, `cellular organisms`) are removed by default before per-cell aggregation, improving microbial specificity.
- **Logging/output readability improved**: unified `[HH:MM:SS] INFO` logging, clearer scRNA labels (e.g., "Cells With Microbe"), and updated dashboard wording.
- **New summary metrics added**: `input_reads` and raw-vs-filtered traceability fields in `pipeline_summary.csv`.
- **Documentation updated**: README test-result tables and `--min-reads` recommendations were simplified and aligned with actual outputs.

### Validation (WSL2 Ubuntu 24.04)

- `celljanus check`: all required tools detected.
- Bulk testdata validated on both `kraken2_testdb` and `standard_8` (including graceful Bracken fallback path).
- scRNA-seq testdata validated on both `kraken2_testdb` and `standard_8` with expected filtered metrics.

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
