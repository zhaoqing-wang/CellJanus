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
import time
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
TAXID_PATTERN = re.compile(r"\(taxid\s+(\d+)\)")

# Alternative: barcode in read name (e.g., @readname_CB_UMI)
UNDERSCORE_PATTERN = re.compile(r"_([ACGTN]{16})_([ACGTN]{10,12})$")

DEFAULT_EXCLUDED_TAXIDS = frozenset({1, 9606, 28384, 131567})
DEFAULT_EXCLUDED_TAXON_NAMES = frozenset(
    {
        "homo sapiens",
        "cellular organisms",
        "root",
        "other sequences",
    }
)


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

    # Taxon filtering
    remove_host_taxa: bool = True  # Remove host/root/non-informative taxa by default

    # Filtering thresholds
    min_reads_per_cell: int = 1  # Minimum reads to consider a cell valid
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

    def valid_cells(self, min_reads: int = 1) -> list[str]:
        """Return cell barcodes that pass the minimum-read filter."""
        return [cb for cb, n in self.cell_total_reads.items() if n >= min_reads]

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
        valid_cells = self.valid_cells(min_reads=min_reads)

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

    def to_normalized_matrix(self, min_reads: int = 1, scale: float = 10000.0) -> pd.DataFrame:
        """
        Convert to normalized cells × species matrix (CPM-like).

        Each cell's counts are normalized to sum to `scale` (default 10,000).
        Suitable for downstream analysis with Seurat/Scanpy.

        Parameters
        ----------
        min_reads : int
            Minimum total reads for a cell to be included.
        scale : float
            Target sum for each cell (default: 10000 for CPM-like).

        Returns
        -------
        DataFrame with normalized abundances.
        """
        raw = self.to_matrix(min_reads=min_reads)
        if len(raw) == 0:
            return raw

        row_sums = raw.sum(axis=1)
        # Avoid division by zero
        row_sums = row_sums.replace(0, 1)
        normalized = raw.div(row_sums, axis=0) * scale

        return normalized

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

    def to_species_summary(self, min_reads: int = 1) -> pd.DataFrame:
        """
        Generate species-level summary statistics.

        Returns DataFrame with columns:
        - species: Species name
        - total_reads: Total reads across all cells
        - total_umis: Total UMIs (if available)
        - n_cells: Number of cells with this species
        - mean_per_cell: Mean count per positive cell
        - prevalence: Fraction of cells with this species
        """
        species_stats = defaultdict(
            lambda: {
                "total_reads": 0,
                "total_umis": 0,
                "n_cells": 0,
                "counts_per_cell": [],
            }
        )

        valid_cells = self.valid_cells(min_reads=min_reads)
        valid_cell_set = set(valid_cells)
        total_cells = len(valid_cells)

        for cb, species_counts in self.cell_species_counts.items():
            if cb not in valid_cell_set:
                continue
            for sp, count in species_counts.items():
                species_stats[sp]["total_reads"] += count
                species_stats[sp]["n_cells"] += 1
                species_stats[sp]["counts_per_cell"].append(count)
                # Count UMIs if available
                if cb in self.cell_species_umis and sp in self.cell_species_umis[cb]:
                    species_stats[sp]["total_umis"] += len(self.cell_species_umis[cb][sp])

        rows = []
        for sp, stats in species_stats.items():
            counts = stats["counts_per_cell"]
            rows.append(
                {
                    "species": sp,
                    "total_reads": stats["total_reads"],
                    "total_umis": stats["total_umis"],
                    "n_cells": stats["n_cells"],
                    "mean_per_cell": np.mean(counts) if counts else 0,
                    "median_per_cell": np.median(counts) if counts else 0,
                    "max_per_cell": max(counts) if counts else 0,
                    "prevalence_pct": 100 * stats["n_cells"] / total_cells
                    if total_cells > 0
                    else 0,
                }
            )

        df = pd.DataFrame(rows)
        if len(df) > 0:
            df = df.sort_values("total_reads", ascending=False).reset_index(drop=True)
        return df

    def to_cell_summary(self, min_reads: int = 1) -> pd.DataFrame:
        """
        Generate cell-level summary statistics.

        Returns DataFrame with columns:
        - cell_barcode: Cell barcode
        - total_microbial_reads: Total microbial reads
        - n_species: Number of species detected
        - dominant_species: Most abundant species
        - dominant_fraction: Fraction of dominant species
        - shannon_diversity: Shannon diversity index
        """
        rows = []

        for cb, species_counts in self.cell_species_counts.items():
            if self.cell_total_reads[cb] < min_reads:
                continue

            total = sum(species_counts.values())
            n_species = len(species_counts)

            if total == 0:
                continue

            # Find dominant species
            dominant_sp = max(species_counts.items(), key=lambda x: x[1])
            dominant_frac = dominant_sp[1] / total

            # Calculate Shannon diversity
            shannon = 0.0
            for count in species_counts.values():
                if count > 0:
                    p = count / total
                    shannon -= p * np.log(p)

            rows.append(
                {
                    "cell_barcode": cb,
                    "total_microbial_reads": total,
                    "n_species": n_species,
                    "dominant_species": dominant_sp[0],
                    "dominant_fraction": round(dominant_frac, 4),
                    "shannon_diversity": round(shannon, 4),
                }
            )

        df = pd.DataFrame(rows)
        if len(df) > 0:
            df = df.sort_values("total_microbial_reads", ascending=False).reset_index(drop=True)
        return df

    def summary(self, min_reads: int = 1) -> dict:
        """Return summary statistics aligned with exported tables."""
        raw_cells = len(self.cell_total_reads)
        valid_cells = self.valid_cells(min_reads=min_reads)
        valid_cell_set = set(valid_cells)

        raw_species = set()
        filtered_species = set()
        for cb, counts in self.cell_species_counts.items():
            raw_species.update(counts.keys())
            if cb in valid_cell_set:
                filtered_species.update(counts.keys())

        raw_total_reads = sum(sum(counts.values()) for counts in self.cell_species_counts.values())
        filtered_total_reads = sum(
            sum(self.cell_species_counts[cb].values()) for cb in valid_cells
        )

        filtered_cells = len(valid_cells)

        return {
            "total_cells": filtered_cells,
            "cells_with_microbes": filtered_cells,
            "cells_filtered_out": raw_cells - filtered_cells,
            "species_detected": len(filtered_species),
            "total_microbial_reads": filtered_total_reads,
            "mean_reads_per_cell": filtered_total_reads / filtered_cells if filtered_cells > 0 else 0,
            "total_cells_raw": raw_cells,
            "species_detected_raw": len(raw_species),
            "total_microbial_reads_raw": raw_total_reads,
            "min_reads_per_cell": min_reads,
        }


def _should_exclude_taxon(taxon: str, *, barcode_cfg: BarcodeConfig) -> bool:
    """Return True if a classification should be excluded before cell aggregation."""
    if not barcode_cfg.remove_host_taxa:
        return False

    taxon_lc = taxon.strip().lower()
    if taxon_lc in DEFAULT_EXCLUDED_TAXON_NAMES:
        return True

    taxid_match = TAXID_PATTERN.search(taxon)
    if taxid_match and int(taxid_match.group(1)) in DEFAULT_EXCLUDED_TAXIDS:
        return True

    return False


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


# ---------------------------------------------------------------------------
# Memory-efficient helpers for large-scale scRNA-seq (v0.2.0)
# ---------------------------------------------------------------------------


def _extract_barcodes_to_file(
    fastq_path: Path,
    output_path: Path,
    *,
    barcode_cfg: Optional[BarcodeConfig] = None,
) -> tuple[int, int]:
    """
    Stream R1 FASTQ and write cell-barcode + UMI to a lightweight temp file.

    One line per read, in input order:
        cell_barcode<TAB>umi
    Empty strings when barcode/UMI is absent.

    Memory: O(1) — only GZIP decompression buffer is held.

    Returns (total_reads, n_with_barcode).
    """
    log = get_logger()
    if barcode_cfg is None:
        barcode_cfg = BarcodeConfig()

    # Select extraction function
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

    cb_len = barcode_cfg.cb_length   # default 16 for 10x v3
    umi_len = barcode_cfg.umi_length  # default 12 for 10x v3
    min_seq_len = cb_len + umi_len    # 28 for 10x v3

    opener = gzip.open if str(fastq_path).endswith(".gz") else open
    output_path.parent.mkdir(parents=True, exist_ok=True)

    total_reads = 0
    n_with_barcode = 0
    n_from_seq = 0

    with opener(fastq_path, "rt") as fh, open(output_path, "w") as out:
        while True:
            header = fh.readline()
            if not header:
                break
            sequence = fh.readline().rstrip("\n")  # keep for sequence-based extraction
            fh.readline()   # + line   — discard
            fh.readline()   # quality  — discard

            # Try header-based extraction first (CB:Z:/UB:Z: tags)
            cb, umi = extract_fn(header)

            # Fallback: extract barcode & UMI from R1 sequence
            # Real 10x Chromium R1: first cb_len bp = cell barcode,
            # next umi_len bp = UMI (e.g. 16+12 = 28 bp for v3)
            if cb is None and len(sequence) >= min_seq_len:
                cb = sequence[:cb_len]
                umi = sequence[cb_len:cb_len + umi_len]
                # Reject if barcode contains N
                if "N" in cb:
                    cb = None
                    umi = None
                else:
                    n_from_seq += 1

            total_reads += 1

            if cb:
                n_with_barcode += 1
                out.write(f"{cb}\t{umi or ''}\n")
            else:
                out.write("\t\n")

            # Progress every 5 million reads
            if total_reads % 5_000_000 == 0:
                log.info(
                    f"[scRNA-seq]   ... {total_reads:,} reads processed "
                    f"({n_with_barcode:,} with barcodes, "
                    f"{100 * n_with_barcode / total_reads:.1f}%)"
                )

    if n_from_seq > 0:
        log.info(
            f"[scRNA-seq] Barcodes extracted from sequence: {n_from_seq:,} "
            f"(header tags: {n_with_barcode - n_from_seq:,})"
        )

    return total_reads, n_with_barcode


def _parse_kraken2_with_barcode_file(
    kraken2_output: Path,
    barcode_file: Path,
    *,
    barcode_cfg: BarcodeConfig,
) -> tuple[CellMicrobialAbundance, dict[str, int]]:
    """
    Stream Kraken2 output and barcode temp file in lockstep.

    Both files have one line per input read in the same order
    (Kraken2 preserves input order), so position-based joining
    uses O(cells × species) memory instead of O(total_reads).
    """
    log = get_logger()
    abundance = CellMicrobialAbundance()

    n_classified = 0
    n_excluded = 0
    n_with_barcode = 0
    n_total = 0

    with open(kraken2_output, "r") as k2_fh, open(barcode_file, "r") as bc_fh:
        for k2_line, bc_line in zip(k2_fh, bc_fh):
            n_total += 1
            parts = k2_line.strip().split("\t")
            if len(parts) < 5 or parts[0] != "C":
                continue

            n_classified += 1
            taxon = parts[2]

            if _should_exclude_taxon(taxon, barcode_cfg=barcode_cfg):
                n_excluded += 1
                continue

            bc_parts = bc_line.strip().split("\t")
            cb = bc_parts[0] if bc_parts[0] else None
            umi = bc_parts[1] if len(bc_parts) > 1 and bc_parts[1] else None

            if cb:
                abundance.add_classification(cb, taxon, umi=umi)
                n_with_barcode += 1

    log.info(
        f"Parsed Kraken2 output: {n_classified:,} classified reads, "
        f"{n_excluded:,} excluded by taxon filter, "
        f"{n_with_barcode:,} retained with cell barcodes"
    )

    return abundance, {
        "classified_reads_raw": n_classified,
        "excluded_taxon_reads": n_excluded,
        "classified_reads_retained": n_classified - n_excluded,
        "classified_reads_with_barcodes": n_with_barcode,
    }


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

    # Step 1: Extract barcodes from R1 → temp file  (streaming, ~0 RAM)
    log.info("[scRNA-seq] Step 1/3: Extracting cell barcodes and UMIs")
    step_t0 = time.perf_counter()
    temp_dir = output_dir / "tmp"
    temp_dir.mkdir(parents=True, exist_ok=True)
    barcode_temp = temp_dir / "barcode_map.tsv"

    total_reads, n_with_barcode = _extract_barcodes_to_file(
        read1, barcode_temp, barcode_cfg=barcode_cfg
    )

    if total_reads == 0:
        raise ValueError(f"No reads found in {read1}")

    log.info(
        f"[scRNA-seq] Extracted barcodes: {n_with_barcode:,} / {total_reads:,} reads "
        f"({100 * n_with_barcode / total_reads:.1f}%)"
    )
    log.info(f"[scRNA-seq] Step 1/3 completed in {time.perf_counter() - step_t0:.1f}s")

    # Step 2: Classify biological reads with Kraken2
    # For scRNA-seq R2 carries the cDNA insert; R1 is barcode+UMI (28 bp
    # for 10x v3 — shorter than Kraken2's default k=35, yielding 0 k-mers).
    # Classifying R2 alone gives identical results and halves Kraken2 I/O.
    log.info("[scRNA-seq] Step 2/3: Running Kraken2 classification")
    step_t0 = time.perf_counter()
    classify_dir = output_dir / "classification"
    classify_input = read2 if read2 else read1
    report_path, output_path = run_kraken2(
        classify_input, kraken2_db, classify_dir, cfg=cfg
    )
    log.info(f"[scRNA-seq] Step 2/3 completed in {time.perf_counter() - step_t0:.1f}s")

    # Step 3: Position-based join  (streaming, ~0 RAM)
    # Both barcode_map.tsv and Kraken2 output have one line per input read
    # in the same order, so we iterate them in lockstep.
    log.info("[scRNA-seq] Step 3/3: Mapping classifications to cell barcodes")
    step_t0 = time.perf_counter()
    abundance, parse_stats = _parse_kraken2_with_barcode_file(
        output_path,
        barcode_temp,
        barcode_cfg=barcode_cfg,
    )
    log.info(f"[scRNA-seq] Step 3/3 completed in {time.perf_counter() - step_t0:.1f}s")

    # Clean up temp barcode file
    try:
        barcode_temp.unlink()
    except OSError:
        pass

    # Export results - comprehensive CSV outputs
    tables_dir = output_dir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    export_t0 = time.perf_counter()

    # 1. Raw count matrix (cells × species)
    matrix_path = tables_dir / "cell_species_counts.csv"
    matrix_df = abundance.to_matrix(min_reads=barcode_cfg.min_reads_per_cell)
    matrix_df.to_csv(matrix_path)
    log.info(
        f"[scRNA-seq] Saved raw count matrix → {matrix_path} ({len(matrix_df)} cells × {len(matrix_df.columns)} species)"
    )

    # 2. Normalized abundance matrix (CPM-like, for Seurat/Scanpy)
    normalized_path = tables_dir / "cell_species_normalized.csv"
    normalized_df = abundance.to_normalized_matrix(min_reads=barcode_cfg.min_reads_per_cell)
    normalized_df.to_csv(normalized_path)
    log.info(f"[scRNA-seq] Saved normalized matrix → {normalized_path}")

    # 3. Long-format table (for ggplot2/seaborn)
    long_path = tables_dir / "cell_species_long.csv"
    long_df = abundance.to_long_format(min_reads=barcode_cfg.min_reads_per_cell)
    long_df.to_csv(long_path, index=False)
    log.info(f"[scRNA-seq] Saved long-format table → {long_path}")

    # 4. Species summary statistics
    species_summary_path = tables_dir / "species_summary.csv"
    species_summary_df = abundance.to_species_summary(min_reads=barcode_cfg.min_reads_per_cell)
    species_summary_df.to_csv(species_summary_path, index=False)
    log.info(
        f"[scRNA-seq] Saved species summary → {species_summary_path} ({len(species_summary_df)} species)"
    )

    # 5. Cell summary statistics
    cell_summary_path = tables_dir / "cell_summary.csv"
    cell_summary_df = abundance.to_cell_summary(min_reads=barcode_cfg.min_reads_per_cell)
    cell_summary_df.to_csv(cell_summary_path, index=False)
    log.info(f"[scRNA-seq] Saved cell summary → {cell_summary_path}")

    # 6. Pipeline summary
    summary = abundance.summary(min_reads=barcode_cfg.min_reads_per_cell)
    summary.update(parse_stats)
    summary_path = tables_dir / "pipeline_summary.csv"
    summary_rows = [{"metric": k, "value": v} for k, v in summary.items()]
    pd.DataFrame(summary_rows).to_csv(summary_path, index=False)
    log.info(f"[scRNA-seq] Table export completed in {time.perf_counter() - export_t0:.1f}s")

    log.info(
        f"[scRNA-seq] Summary: {summary['total_cells']:,} cells passed min-reads filter, "
        f"{summary['species_detected']} species detected"
    )

    return {
        "abundance": abundance,
        "matrix_path": matrix_path,
        "normalized_path": normalized_path,
        "long_path": long_path,
        "species_summary_path": species_summary_path,
        "cell_summary_path": cell_summary_path,
        "matrix_df": matrix_df,
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
