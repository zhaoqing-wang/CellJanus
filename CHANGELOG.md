# Changelog

All notable changes to CellJanus will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
