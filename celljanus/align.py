"""
Host genome alignment module.

Aligns QC'd FASTQ reads against a host reference genome (typically hg38)
using Bowtie2.  The module produces a sorted BAM file for downstream
unmapped-read extraction.

Key design choices:
  • ``--very-sensitive`` preset for thorough host removal.
  • ``--no-unal`` flag to skip writing unaligned reads to the BAM
    (they are extracted separately via samtools in the extract module).
  • Multi-threaded via ``-p`` for Bowtie2 and ``-@ `` for samtools sort.
  • Streams Bowtie2 SAM output directly into ``samtools sort`` via a pipe
    to avoid writing a large intermediate SAM file to disk.
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Optional

from celljanus.config import CellJanusConfig, require_tool
from celljanus.utils import file_size_human, get_logger, run_cmd, fmt_elapsed

import time


def align_to_host(
    read1: Path,
    host_index_prefix: Path,
    output_dir: Path,
    *,
    read2: Optional[Path] = None,
    cfg: Optional[CellJanusConfig] = None,
) -> Path:
    """
    Align reads to the host genome and produce a coordinate-sorted BAM.

    Parameters
    ----------
    read1 : Path
        QC'd R1 FASTQ (.fastq.gz).
    host_index_prefix : Path
        Bowtie2 index prefix (e.g. ``/db/hg38/GRCh38_noalt_as``).
    output_dir : Path
        Directory for output BAM and alignment stats.
    read2 : Path, optional
        QC'd R2 FASTQ for paired-end data.
    cfg : CellJanusConfig, optional
        Configuration object.

    Returns
    -------
    Path to the sorted BAM file (``<output_dir>/host_aligned.sorted.bam``).
    """
    log = get_logger()
    if cfg is None:
        cfg = CellJanusConfig()

    bowtie2 = require_tool("bowtie2")
    samtools = require_tool("samtools")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    sorted_bam = output_dir / "host_aligned.sorted.bam"
    stats_file = output_dir / "host_align_stats.txt"

    if sorted_bam.exists():
        log.info(f"Sorted BAM already exists, skipping alignment: {sorted_bam.name}")
        return sorted_bam

    threads = cfg.threads

    # ----- Build Bowtie2 command -----
    bt2_cmd = [
        bowtie2,
        "-x",
        str(host_index_prefix),
        "-p",
        str(threads),
        "--very-sensitive",
        "--seed",
        "42",
        "--met-stderr",
    ]
    if read2 is not None:
        bt2_cmd.extend(["-1", str(read1), "-2", str(read2)])
    else:
        bt2_cmd.extend(["-U", str(read1)])

    # Extra user flags
    if cfg.align_extra_args:
        for arg in cfg.align_extra_args.split():
            if arg not in bt2_cmd:  # avoid duplicating --very-sensitive etc.
                bt2_cmd.append(arg)

    # ----- Pipe Bowtie2 → samtools sort (avoids intermediate SAM) -----
    sort_cmd = [
        samtools,
        "sort",
        "-@",
        str(threads),
        "-m",
        f"{max(1, int(cfg.max_memory_gb / threads))}G",
        "-o",
        str(sorted_bam),
        "-",  # read from stdin
    ]

    log.info(
        f"Aligning {read1.name} ({file_size_human(read1)}) "
        f"to host index {host_index_prefix.name}  [threads={threads}]"
    )
    start = time.perf_counter()

    # Use shell pipe: bowtie2 | samtools sort
    bt2_proc = subprocess.Popen(
        [str(c) for c in bt2_cmd],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    sort_proc = subprocess.Popen(
        [str(c) for c in sort_cmd],
        stdin=bt2_proc.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    # Allow bt2_proc to receive SIGPIPE if sort_proc exits
    bt2_proc.stdout.close()

    sort_stdout, sort_stderr = sort_proc.communicate()
    bt2_stderr_bytes = bt2_proc.stderr.read()
    bt2_proc.wait()

    elapsed = time.perf_counter() - start
    bt2_stderr = bt2_stderr_bytes.decode("utf-8", errors="replace")

    # Save alignment stats
    stats_file.write_text(bt2_stderr, encoding="utf-8")
    log.info(f"Host alignment completed in {fmt_elapsed(elapsed)}")

    # Parse overall alignment rate from Bowtie2 stderr
    for line in bt2_stderr.splitlines():
        if "overall alignment rate" in line:
            log.info(f"  {line.strip()}")
            break

    if sort_proc.returncode != 0:
        err_msg = sort_stderr.decode("utf-8", errors="replace") if sort_stderr else ""
        raise RuntimeError(f"samtools sort failed:\n{err_msg}")

    # ----- Index the BAM -----
    run_cmd(
        [samtools, "index", "-@", str(threads), str(sorted_bam)],
        desc="Indexing BAM",
    )

    log.info(f"Sorted BAM: {sorted_bam}  ({file_size_human(sorted_bam)})")
    return sorted_bam


def get_alignment_stats(stats_file: Path) -> dict:
    """
    Parse Bowtie2 alignment statistics from the saved stderr output.

    Returns a dict with keys like 'total_reads', 'aligned_0_times',
    'aligned_1_time', 'aligned_gt1_times', 'overall_alignment_rate'.
    """
    text = stats_file.read_text(encoding="utf-8")
    stats: dict = {}
    for line in text.splitlines():
        line = line.strip()
        if "reads; of these:" in line:
            stats["total_reads"] = int(line.split()[0])
        elif "aligned concordantly 0 times" in line or "aligned 0 times" in line:
            stats["aligned_0_times"] = int(line.split()[0])
        elif "aligned concordantly exactly 1 time" in line or "aligned exactly 1 time" in line:
            stats["aligned_1_time"] = int(line.split()[0])
        elif "aligned concordantly >1 times" in line or "aligned >1 times" in line:
            stats["aligned_gt1_times"] = int(line.split()[0])
        elif "overall alignment rate" in line:
            stats["overall_alignment_rate"] = line.split()[0]
    return stats
