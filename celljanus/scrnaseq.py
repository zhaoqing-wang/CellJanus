"""
Single-cell RNA-seq barcode processing module.

Extracts cell barcodes (CB) and UMIs (UB) from 10x Genomics-style FASTQ
reads and tracks per-cell microbial abundance for scRNA-seq integration.

Supports:
  • 10x Genomics 3' and 5' chemistry (CB:Z: / UB:Z: tags)
  • Parse-style barcodes (r1:linker:CB_UMI)
  • Custom barcode patterns via regex

Performance optimisations for WSL2:
  • Streaming processing to minimise memory footprint
  • Batched I/O with configurable chunk sizes
  • Support for native Linux paths to avoid cross-filesystem overhead
"""

from __future__ import annotations

import gzip
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator, Optional

import pandas as pd
import numpy as np

from celljanus.config import CellJanusConfig
from celljanus.utils import get_logger, file_size_human


# ---------------------------------------------------------------------------
# Barcode patterns for different scRNA-seq platforms
# ---------------------------------------------------------------------------

# 10x Genomics: CB:Z:AGCTAGCTAG UB:Z:ACTGACTG in read header comment
TENX_CB_PATTERN = re.compile(r"CB:Z:([ACGT]+)")
TENX_UB_PATTERN = re.compile(r"UB:Z:([ACGT]+)")

# Alternative: barcode in read name (e.g., @readname_CB_UMI)
UNDERSCORE_PATTERN = re.compile(r"_([ACGTN]{16})_([ACGTN]{10,12})$")


@dataclass
class BarcodeConfig:
    """Configuration for barcode extraction from scRNA-seq FASTQ."""

    # Barcode detection mode
    mode: str = "10x"  # "10x", "parse", "custom", "auto"

    # Custom regex patterns (when mode="custom")
    cb_pattern: Optional[str] = None  # regex for cell barcode
    ub_pattern: Optional[str] = None  # regex for UMI

    # 10x barcode parameters
    cb_length: int = 16  # Cell barcode length (16 for v3, 14 for v2)
    umi_length: int = 12  # UMI length (12 for v3, 10 for v2)

    # Whitelist filtering
    whitelist_path: Optional[Path] = None  # Path to barcode whitelist

    # Filtering thresholds
    min_reads_per_cell: int = 10  # Minimum reads to consider a cell valid
    min_umi_per_species: int = 1  # Minimum UMIs per species per cell


@dataclass
class CellBarcodeRead:
    """A FASTQ read with extracted cell barcode and UMI."""

    read_id: str
    sequence: str
    quality: str
    cell_barcode: Optional[str] = None
    umi: Optional[str] = None

    @property
    def has_barcode(self) -> bool:
        return self.cell_barcode is not None


def extract_barcode_10x(header: str) -> tuple[Optional[str], Optional[str]]:
    """
    Extract cell barcode and UMI from 10x Genomics-style header.

    Format: @read_id CB:Z:BARCODE UB:Z:UMI
    Also handles: @read_id:CB:UMI format
    """
    cb_match = TENX_CB_PATTERN.search(header)
    ub_match = TENX_UB_PATTERN.search(header)

    cb = cb_match.group(1) if cb_match else None
    umi = ub_match.group(1) if ub_match else None

    # Fallback: try underscore format (@readname_CB_UMI)
    if cb is None:
        underscore_match = UNDERSCORE_PATTERN.search(header.split()[0])
        if underscore_match:
            cb = underscore_match.group(1)
            umi = underscore_match.group(2)

    return cb, umi


def extract_barcode_auto(header: str) -> tuple[Optional[str], Optional[str]]:
    """
    Auto-detect barcode format and extract.

    Tries multiple patterns in order of likelihood.
    """
    # Try 10x format first (most common)
    cb, umi = extract_barcode_10x(header)
    if cb:
        return cb, umi

    # Try Parse Biosciences format: @read:sublibrary:CB_UMI:...
    parts = header.split(":")
    if len(parts) >= 4:
        # Look for 16bp barcode followed by 10-12bp UMI
        for part in parts:
            if len(part) >= 26 and all(c in "ACGTN" for c in part[:28]):
                return part[:16], part[16:28]

    return None, None


def parse_fastq_with_barcodes(
    fastq_path: Path,
    *,
    barcode_cfg: Optional[BarcodeConfig] = None,
    chunk_size: int = 100_000,
) -> Iterator[list[CellBarcodeRead]]:
    """
    Stream FASTQ file and extract cell barcodes in chunks.

    Yields chunks of CellBarcodeRead objects for memory-efficient
    processing of large scRNA-seq datasets.

    Parameters
    ----------
    fastq_path : Path
        Input FASTQ file (gzipped supported).
    barcode_cfg : BarcodeConfig
        Barcode extraction configuration.
    chunk_size : int
        Number of reads per chunk.

    Yields
    ------
    List of CellBarcodeRead objects.
    """
    if barcode_cfg is None:
        barcode_cfg = BarcodeConfig()

    # Select extraction function based on mode
    if barcode_cfg.mode == "10x":
        extract_fn = extract_barcode_10x
    elif barcode_cfg.mode == "auto":
        extract_fn = extract_barcode_auto
    elif barcode_cfg.mode == "custom" and barcode_cfg.cb_pattern:
        cb_re = re.compile(barcode_cfg.cb_pattern)
        ub_re = re.compile(barcode_cfg.ub_pattern) if barcode_cfg.ub_pattern else None

        def extract_fn(header):
            cb_m = cb_re.search(header)
            ub_m = ub_re.search(header) if ub_re else None
            return (cb_m.group(1) if cb_m else None, ub_m.group(1) if ub_m else None)
    else:
        extract_fn = extract_barcode_auto

    # Open file (handle gzip)
    opener = gzip.open if str(fastq_path).endswith(".gz") else open

    chunk: list[CellBarcodeRead] = []
    with opener(fastq_path, "rt") as fh:
        while True:
            # Read 4 lines (one FASTQ record)
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().strip()
            plus = fh.readline()  # + line
            qual = fh.readline().strip()

            # Extract barcode
            cb, umi = extract_fn(header)
            read_id = header.split()[0][1:]  # Remove @ prefix

            chunk.append(
                CellBarcodeRead(
                    read_id=read_id,
                    sequence=seq,
                    quality=qual,
                    cell_barcode=cb,
                    umi=umi,
                )
            )

            if len(chunk) >= chunk_size:
                yield chunk
                chunk = []

        if chunk:
            yield chunk


@dataclass
class CellMicrobialAbundance:
    """Per-cell microbial abundance data structure."""

    # cell_barcode -> species -> count
    cell_species_counts: dict[str, dict[str, int]] = field(
        default_factory=lambda: defaultdict(lambda: defaultdict(int))
    )

    # cell_barcode -> species -> set of UMIs (for deduplication)
    cell_species_umis: dict[str, dict[str, set]] = field(
        default_factory=lambda: defaultdict(lambda: defaultdict(set))
    )

    # Total reads per cell
    cell_total_reads: dict[str, int] = field(default_factory=lambda: defaultdict(int))

    # Cell barcodes without microbial reads
    cells_without_microbes: set[str] = field(default_factory=set)

    def add_classification(
        self,
        cell_barcode: str,
        species: str,
        umi: Optional[str] = None,
    ) -> None:
        """Record a microbial classification for a cell."""
        self.cell_total_reads[cell_barcode] += 1

        if umi:
            # UMI-based deduplication
            if umi not in self.cell_species_umis[cell_barcode][species]:
                self.cell_species_umis[cell_barcode][species].add(umi)
                self.cell_species_counts[cell_barcode][species] += 1
        else:
            # No UMI - just count reads
            self.cell_species_counts[cell_barcode][species] += 1

    def to_matrix(self, min_reads: int = 1) -> pd.DataFrame:
        """
        Convert to cells × species abundance matrix.

        Parameters
        ----------
        min_reads : int
            Minimum total reads for a cell to be included.

        Returns
        -------
        DataFrame with cells as rows and species as columns.
        """
        # Filter cells by minimum reads
        valid_cells = [cb for cb, n in self.cell_total_reads.items() if n >= min_reads]

        # Collect all species
        all_species = set()
        for cb in valid_cells:
            all_species.update(self.cell_species_counts[cb].keys())

        # Build matrix
        data = []
        for cb in valid_cells:
            row = {"cell_barcode": cb}
            for sp in all_species:
                row[sp] = self.cell_species_counts[cb].get(sp, 0)
            data.append(row)

        if not data:
            return pd.DataFrame(columns=["cell_barcode"])

        df = pd.DataFrame(data)
        df = df.set_index("cell_barcode")

        # Sort columns by total abundance
        col_sums = df.sum().sort_values(ascending=False)
        df = df[col_sums.index]

        return df

    def to_long_format(self, min_reads: int = 1) -> pd.DataFrame:
        """
        Convert to long-format DataFrame for plotting.

        Columns: cell_barcode, species, count, fraction
        """
        rows = []
        for cb, species_counts in self.cell_species_counts.items():
            if self.cell_total_reads[cb] < min_reads:
                continue
            total = sum(species_counts.values())
            for sp, count in species_counts.items():
                rows.append(
                    {
                        "cell_barcode": cb,
                        "species": sp,
                        "count": count,
                        "fraction": count / total if total > 0 else 0,
                    }
                )

        return pd.DataFrame(rows)

    def summary(self) -> dict:
        """Return summary statistics."""
        n_cells = len(self.cell_total_reads)
        n_cells_with_microbes = len(
            [cb for cb, counts in self.cell_species_counts.items() if counts]
        )

        all_species = set()
        for counts in self.cell_species_counts.values():
            all_species.update(counts.keys())

        total_microbial_reads = sum(
            sum(counts.values()) for counts in self.cell_species_counts.values()
        )

        return {
            "total_cells": n_cells,
            "cells_with_microbes": n_cells_with_microbes,
            "cells_without_microbes": n_cells - n_cells_with_microbes,
            "species_detected": len(all_species),
            "total_microbial_reads": total_microbial_reads,
            "mean_reads_per_cell": total_microbial_reads / n_cells if n_cells > 0 else 0,
        }


def parse_kraken2_output_with_barcodes(
    kraken2_output: Path,
    barcode_reads: dict[str, CellBarcodeRead],
) -> CellMicrobialAbundance:
    """
    Parse Kraken2 output and map classifications to cell barcodes.

    Parameters
    ----------
    kraken2_output : Path
        Kraken2 --output file (tab-separated).
    barcode_reads : dict
        Mapping of read_id -> CellBarcodeRead.

    Returns
    -------
    CellMicrobialAbundance with per-cell species counts.
    """
    log = get_logger()
    abundance = CellMicrobialAbundance()

    n_classified = 0
    n_with_barcode = 0

    with open(kraken2_output, "r") as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue

            classified = parts[0] == "C"
            read_id = parts[1]
            taxon = parts[2]  # species name or taxid

            if not classified:
                continue

            n_classified += 1

            # Look up barcode for this read
            if read_id in barcode_reads:
                br = barcode_reads[read_id]
                if br.cell_barcode:
                    abundance.add_classification(
                        br.cell_barcode,
                        taxon,
                        umi=br.umi,
                    )
                    n_with_barcode += 1

    log.info(
        f"Parsed Kraken2 output: {n_classified:,} classified reads, "
        f"{n_with_barcode:,} with cell barcodes"
    )

    return abundance


def run_scrnaseq_classification(
    read1: Path,
    kraken2_db: Path,
    output_dir: Path,
    *,
    read2: Optional[Path] = None,
    barcode_cfg: Optional[BarcodeConfig] = None,
    cfg: Optional[CellJanusConfig] = None,
) -> dict:
    """
    Run full scRNA-seq microbial classification pipeline.

    Pipeline:
    1. Parse FASTQ and extract cell barcodes
    2. Run Kraken2 classification
    3. Map classifications to cell barcodes
    4. Generate per-cell abundance matrix

    Parameters
    ----------
    read1 : Path
        R1 FASTQ with cell barcodes (typically from 10x).
    kraken2_db : Path
        Kraken2 database path.
    output_dir : Path
        Output directory.
    read2 : Path, optional
        R2 FASTQ for paired-end data.
    barcode_cfg : BarcodeConfig
        Barcode extraction settings.
    cfg : CellJanusConfig
        Pipeline configuration.

    Returns
    -------
    dict with keys:
        - abundance: CellMicrobialAbundance object
        - matrix_path: Path to cells × species CSV
        - long_path: Path to long-format CSV
        - summary: dict of statistics
    """
    from celljanus.classify import run_kraken2

    log = get_logger()
    if cfg is None:
        cfg = CellJanusConfig()
    if barcode_cfg is None:
        barcode_cfg = BarcodeConfig()

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    log.info(f"[scRNA-seq] Processing {read1.name} ({file_size_human(read1)})")

    # Step 1: Parse barcodes from R1 (or R2 if R1 has no barcodes)
    log.info("[scRNA-seq] Step 1/3: Extracting cell barcodes and UMIs")
    barcode_reads: dict[str, CellBarcodeRead] = {}
    n_with_barcode = 0
    total_reads = 0

    for chunk in parse_fastq_with_barcodes(
        read1, barcode_cfg=barcode_cfg, chunk_size=cfg.chunk_size
    ):
        for br in chunk:
            total_reads += 1
            barcode_reads[br.read_id] = br
            if br.has_barcode:
                n_with_barcode += 1

    log.info(
        f"[scRNA-seq] Extracted barcodes: {n_with_barcode:,} / {total_reads:,} reads ({100 * n_with_barcode / total_reads:.1f}%)"
    )

    # Step 2: Run Kraken2 classification
    log.info("[scRNA-seq] Step 2/3: Running Kraken2 classification")
    classify_dir = output_dir / "classification"
    report_path, output_path = run_kraken2(read1, kraken2_db, classify_dir, read2=read2, cfg=cfg)

    # Step 3: Map classifications to barcodes
    log.info("[scRNA-seq] Step 3/3: Mapping classifications to cell barcodes")
    abundance = parse_kraken2_output_with_barcodes(output_path, barcode_reads)

    # Export results
    matrix_path = output_dir / "cell_species_matrix.csv"
    long_path = output_dir / "cell_species_long.csv"

    matrix_df = abundance.to_matrix(min_reads=barcode_cfg.min_reads_per_cell)
    matrix_df.to_csv(matrix_path)
    log.info(f"[scRNA-seq] Saved cell × species matrix → {matrix_path} ({len(matrix_df)} cells)")

    long_df = abundance.to_long_format(min_reads=barcode_cfg.min_reads_per_cell)
    long_df.to_csv(long_path, index=False)
    log.info(f"[scRNA-seq] Saved long-format table → {long_path}")

    summary = abundance.summary()
    log.info(
        f"[scRNA-seq] Summary: {summary['cells_with_microbes']} cells with microbes, "
        f"{summary['species_detected']} species detected"
    )

    return {
        "abundance": abundance,
        "matrix_path": matrix_path,
        "long_path": long_path,
        "barcode_reads": barcode_reads,
        "summary": summary,
    }


# ---------------------------------------------------------------------------
# WSL2 I/O optimisation
# ---------------------------------------------------------------------------


def detect_wsl2() -> bool:
    """Detect if running on WSL2."""
    import platform
    import os

    if platform.system() != "Linux":
        return False

    # Check for WSL indicators
    try:
        with open("/proc/version", "r") as f:
            version = f.read().lower()
            return "microsoft" in version or "wsl" in version
    except FileNotFoundError:
        return False


def is_cross_filesystem_path(path: Path) -> bool:
    """
    Check if path crosses Windows/Linux filesystem boundary in WSL2.

    Paths under /mnt/c, /mnt/d, etc. are Windows filesystems with
    slower I/O performance.

    Note: This function is designed to run inside WSL2. On native Windows,
    it will always return False as the path resolution differs.
    """
    import platform

    # On native Windows, this check doesn't apply
    if platform.system() == "Windows":
        return False

    # Use string representation, try both resolved and original
    path_str = str(path)
    try:
        resolved_str = str(path.resolve())
    except (OSError, ValueError):
        resolved_str = path_str

    return path_str.startswith("/mnt/") or resolved_str.startswith("/mnt/")


def recommend_native_path(path: Path, workspace: str = "/home") -> Path:
    """
    Suggest a native Linux path for better I/O in WSL2.

    Parameters
    ----------
    path : Path
        Original path (possibly on Windows filesystem).
    workspace : str
        Preferred Linux workspace directory.

    Returns
    -------
    Suggested native Linux path.
    """
    if not is_cross_filesystem_path(path):
        return path

    # Extract the filename/relative part
    name = path.name

    # Suggest path in Linux filesystem
    return Path(workspace) / "celljanus_work" / name


def wsl2_io_warning(paths: list[Path]) -> Optional[str]:
    """
    Generate I/O performance warning if using Windows filesystem in WSL2.

    Returns warning message or None if no issue.
    """
    if not detect_wsl2():
        return None

    cross_fs = [p for p in paths if is_cross_filesystem_path(p)]
    if not cross_fs:
        return None

    return (
        f"⚠️  Performance Warning: {len(cross_fs)} path(s) are on Windows filesystem "
        f"(/mnt/...). For better I/O performance, copy data to a native Linux path "
        f"(e.g., /home/user/data/) before running CellJanus."
    )
