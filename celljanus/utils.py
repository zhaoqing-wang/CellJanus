"""
Shared utility helpers for CellJanus.

Covers subprocess execution with real-time logging, FASTQ streaming,
file-size helpers, and elapsed-time formatting.
"""

from __future__ import annotations

import gzip
import logging
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import IO, Iterator, Optional, Sequence

from rich.console import Console
from rich.logging import RichHandler

console = Console(stderr=True)

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

_logger: Optional[logging.Logger] = None


def get_logger(log_file: Optional[Path] = None) -> logging.Logger:
    """Return (and lazily configure) the package-wide logger."""
    global _logger
    if _logger is not None:
        return _logger

    _logger = logging.getLogger("celljanus")
    _logger.setLevel(logging.DEBUG)

    # Rich console handler (INFO+)
    rh = RichHandler(console=console, show_path=False, markup=True)
    rh.setLevel(logging.INFO)
    _logger.addHandler(rh)

    # File handler (DEBUG+) — optional
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(log_file, encoding="utf-8")
        fh.setLevel(logging.DEBUG)
        fmt = logging.Formatter("%(asctime)s  %(levelname)-8s  %(message)s")
        fh.setFormatter(fmt)
        _logger.addHandler(fh)

    return _logger


# ---------------------------------------------------------------------------
# Subprocess runner
# ---------------------------------------------------------------------------


def run_cmd(
    cmd: Sequence[str],
    *,
    desc: str = "",
    capture: bool = False,
    check: bool = True,
    env: Optional[dict] = None,
    cwd: Optional[Path] = None,
    timeout: Optional[int] = None,
) -> subprocess.CompletedProcess:
    """
    Execute an external command with logging and error handling.

    Parameters
    ----------
    cmd : list of str
        The command and its arguments.
    desc : str
        Human-readable description printed before execution.
    capture : bool
        If True, capture stdout/stderr and return them.
    check : bool
        If True, raise on non-zero exit code.
    env : dict, optional
        Extra environment variables (merged with os.environ).
    cwd : Path, optional
        Working directory.
    timeout : int, optional
        Maximum seconds to wait.

    Returns
    -------
    subprocess.CompletedProcess
    """
    log = get_logger()

    cmd_str = " ".join(str(c) for c in cmd)
    if desc:
        log.info(f"[bold cyan]{desc}[/bold cyan]")
    log.debug(f"CMD: {cmd_str}")

    merged_env = {**os.environ, **(env or {})}
    start = time.perf_counter()

    try:
        result = subprocess.run(
            [str(c) for c in cmd],
            stdout=subprocess.PIPE if capture else None,
            stderr=subprocess.PIPE if capture else subprocess.STDOUT,
            text=True,
            check=check,
            env=merged_env,
            cwd=str(cwd) if cwd else None,
            timeout=timeout,
        )
    except FileNotFoundError:
        log.error(f"Command not found: {cmd[0]}")
        raise
    except subprocess.CalledProcessError as exc:
        log.error(f"Command failed (exit {exc.returncode}): {cmd_str}")
        if exc.stderr:
            log.error(exc.stderr[:2000])
        raise
    except subprocess.TimeoutExpired:
        log.error(f"Command timed out after {timeout}s: {cmd_str}")
        raise

    elapsed = time.perf_counter() - start
    log.debug(f"Finished in {fmt_elapsed(elapsed)}")
    return result


# ---------------------------------------------------------------------------
# FASTQ streaming helpers (memory-efficient)
# ---------------------------------------------------------------------------


def count_fastq_reads(path: Path) -> int:
    """Count reads in a FASTQ file (4 lines per record). Handles .gz."""
    opener = gzip.open if str(path).endswith(".gz") else open
    n_lines = 0
    with opener(path, "rt") as fh:
        for _ in fh:
            n_lines += 1
    return n_lines // 4


def iter_fastq_chunks(
    path: Path,
    chunk_size: int = 500_000,
) -> Iterator[list[str]]:
    """
    Yield chunks of *chunk_size* FASTQ records (each record = 4 lines)
    as a list of raw lines.  Memory usage ≈ chunk_size × avg_line_len.
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    buf: list[str] = []
    record_count = 0
    with opener(path, "rt") as fh:
        for line in fh:
            buf.append(line)
            if len(buf) % 4 == 0:
                record_count += 1
                if record_count >= chunk_size:
                    yield buf
                    buf = []
                    record_count = 0
    if buf:
        yield buf


# ---------------------------------------------------------------------------
# File helpers
# ---------------------------------------------------------------------------


def file_size_human(path: Path) -> str:
    """Return human-readable file size string."""
    size = path.stat().st_size
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if size < 1024:
            return f"{size:.1f} {unit}"
        size /= 1024
    return f"{size:.1f} PB"


def fmt_elapsed(seconds: float) -> str:
    """Format seconds into H:MM:SS or M:SS."""
    h, rem = divmod(int(seconds), 3600)
    m, s = divmod(rem, 60)
    if h:
        return f"{h}:{m:02d}:{s:02d}"
    return f"{m}:{s:02d}"


def estimate_time(total_reads: int, reads_per_sec: float) -> str:
    """Estimate remaining time given throughput."""
    if reads_per_sec <= 0:
        return "unknown"
    return fmt_elapsed(total_reads / reads_per_sec)


def ensure_parent(path: Path) -> Path:
    """Create parent directories and return *path*."""
    path.parent.mkdir(parents=True, exist_ok=True)
    return path
