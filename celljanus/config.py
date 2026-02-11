"""
Configuration management for CellJanus.

Centralises default parameters, external tool paths, resource limits,
and pipeline-wide settings so that every module shares a single source
of truth.
"""

from __future__ import annotations

import multiprocessing
import os
import platform
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import psutil


# ---------------------------------------------------------------------------
# External tool discovery
# ---------------------------------------------------------------------------


def _find_tool(name: str) -> Optional[str]:
    """Return the absolute path to *name* if it is on PATH, else None."""
    return shutil.which(name)


def find_tools() -> dict[str, Optional[str]]:
    """Scan PATH for every external binary CellJanus may call."""
    names = [
        "fastp",
        "bowtie2",
        "bowtie2-build",
        "samtools",
        "kraken2",
        "kraken2-build",
        "bracken",
        "bracken-build",
        "star",
        "STAR",
    ]
    found: dict[str, Optional[str]] = {}
    for n in names:
        found[n] = _find_tool(n)
    # Normalise STAR (case-insensitive on Windows)
    if found.get("star") is None and found.get("STAR"):
        found["star"] = found["STAR"]
    return found


def require_tool(name: str) -> str:
    """Return the path to *name* or raise with a helpful message."""
    path = _find_tool(name)
    if path is None:
        raise EnvironmentError(
            f"Required external tool '{name}' was not found on PATH.\n"
            f"Please install it and make sure it is accessible.\n"
            f"  - fastp:    https://github.com/OpenGene/fastp\n"
            f"  - bowtie2:  https://github.com/BenLangmead/bowtie2\n"
            f"  - samtools: https://github.com/samtools/samtools\n"
            f"  - kraken2:  https://github.com/DerrickWood/kraken2\n"
            f"  - bracken:  https://github.com/jenniferlu717/Bracken\n"
        )
    return path


# ---------------------------------------------------------------------------
# System resource helpers
# ---------------------------------------------------------------------------


def available_memory_gb() -> float:
    """Return available physical memory in GiB."""
    return psutil.virtual_memory().available / (1024**3)


def total_memory_gb() -> float:
    """Return total physical memory in GiB."""
    return psutil.virtual_memory().total / (1024**3)


def default_threads() -> int:
    """Sensible default thread count (leave 1–2 cores free)."""
    n = multiprocessing.cpu_count()
    return max(1, n - 2)


# ---------------------------------------------------------------------------
# Pipeline configuration dataclass
# ---------------------------------------------------------------------------


@dataclass
class CellJanusConfig:
    """Master configuration object passed through the pipeline."""

    # --- I/O ---
    output_dir: Path = field(default_factory=lambda: Path("celljanus_output"))
    temp_dir: Optional[Path] = None  # defaults to output_dir / "tmp"
    log_file: Optional[Path] = None  # defaults to output_dir / "celljanus.log"

    # --- Computing resources ---
    threads: int = field(default_factory=default_threads)
    max_memory_gb: float = field(default_factory=lambda: min(total_memory_gb() * 0.8, 28.0))
    chunk_size: int = 500_000  # FASTQ records per chunk for streaming

    # --- QC (fastp) ---
    qc_min_quality: int = 15  # Phred; gentle to keep microbial reads
    qc_min_length: int = 36  # bp; very permissive
    qc_complexity_threshold: float = 0.0  # 0 = disabled; avoids discarding novel seqs
    qc_adapter_trimming: bool = True
    qc_extra_args: str = ""  # arbitrary fastp flags

    # --- Alignment (Bowtie2 to host genome) ---
    aligner: str = "bowtie2"  # "bowtie2" or "star"
    host_index: Optional[Path] = None  # path to pre-built index prefix
    align_extra_args: str = "--very-sensitive --no-unal"
    align_seed_length: int = 22

    # --- Microbial classification (Kraken2 + Bracken) ---
    kraken2_db: Optional[Path] = None
    kraken2_confidence: float = 0.05
    bracken_read_length: int = 150
    bracken_level: str = "S"  # S=species, G=genus, etc.

    # --- Visualisation ---
    top_n_species: int = 20
    plot_format: str = "png"  # png | pdf | svg
    dpi: int = 300

    def __post_init__(self) -> None:
        self.output_dir = Path(self.output_dir)
        if self.temp_dir is None:
            self.temp_dir = self.output_dir / "tmp"
        else:
            self.temp_dir = Path(self.temp_dir)
        if self.log_file is None:
            self.log_file = self.output_dir / "celljanus.log"
        else:
            self.log_file = Path(self.log_file)

    def ensure_dirs(self) -> None:
        """Create output and temp directories if they do not exist."""
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.temp_dir.mkdir(parents=True, exist_ok=True)

    @property
    def system_summary(self) -> str:
        mem = total_memory_gb()
        cpu = multiprocessing.cpu_count()
        return (
            f"OS={platform.system()} {platform.release()}  "
            f"CPUs={cpu}  RAM={mem:.1f} GiB  "
            f"Threads={self.threads}  MaxMem={self.max_memory_gb:.1f} GiB"
        )
