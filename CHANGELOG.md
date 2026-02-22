# Changelog

All notable changes to CellJanus will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
