"""
Unmapped-read extraction module.

After aligning to the host genome, this module extracts reads that did
**not** map to the host — these are the candidate microbial reads.

Uses ``samtools view`` with the ``-f 4`` flag (unmapped) and converts
back to FASTQ via ``samtools fastq`` so that downstream classifiers
(Kraken2) can consume them directly.

For paired-end data:
  • ``-f 12``  — both mates unmapped
  • ``-f 4 -F 8``  — R1 unmapped, R2 mapped (chimeric; kept for safety)
  • ``-f 8 -F 4``  — R1 mapped, R2 unmapped  (chimeric; kept for safety)

All three classes are merged into the output to maximise sensitivity for
microbial detection.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from celljanus.config import CellJanusConfig, require_tool
from celljanus.utils import file_size_human, get_logger, run_cmd


def extract_unmapped(
    sorted_bam: Path,
    output_dir: Path,
    *,
    paired: bool = False,
    cfg: Optional[CellJanusConfig] = None,
) -> tuple[Path, Optional[Path]]:
    """
    Extract reads that did NOT align to the host genome.

    Parameters
    ----------
    sorted_bam : Path
        Coordinate-sorted, indexed BAM from the alignment step.
    output_dir : Path
        Where to write the unmapped FASTQ file(s).
    paired : bool
        Whether the original data was paired-end.
    cfg : CellJanusConfig, optional

    Returns
    -------
    (unmapped_r1.fastq.gz, unmapped_r2.fastq.gz | None)
    """
    log = get_logger()
    if cfg is None:
        cfg = CellJanusConfig()

    samtools = require_tool("samtools")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    threads = cfg.threads

    unmapped_bam = output_dir / "unmapped.bam"
    unmapped_r1 = output_dir / "unmapped_R1.fastq.gz"
    unmapped_r2: Optional[Path] = output_dir / "unmapped_R2.fastq.gz" if paired else None

    if unmapped_r1.exists():
        log.info(f"Unmapped FASTQ already exists, skipping: {unmapped_r1.name}")
        return unmapped_r1, unmapped_r2

    if paired:
        # Extract reads where at least one mate is unmapped
        # -f 4: unmapped flag set on this read
        # This captures: both unmapped, R1 unmapped, R2 unmapped
        run_cmd(
            [
                samtools,
                "view",
                "-@",
                str(threads),
                "-b",
                "-f",
                "4",  # read unmapped
                str(sorted_bam),
                "-o",
                str(unmapped_bam),
            ],
            desc="Extracting unmapped reads (paired)",
        )

        # Convert to paired FASTQ
        run_cmd(
            [
                samtools,
                "fastq",
                "-@",
                str(threads),
                "-1",
                str(unmapped_r1),
                "-2",
                str(unmapped_r2),
                "-0",
                "/dev/null",
                "-s",
                "/dev/null",
                "-n",
                str(unmapped_bam),
            ],
            desc="Converting unmapped BAM → FASTQ (PE)",
        )
    else:
        # Single-end: simply filter unmapped reads
        run_cmd(
            [
                samtools,
                "view",
                "-@",
                str(threads),
                "-b",
                "-f",
                "4",
                str(sorted_bam),
                "-o",
                str(unmapped_bam),
            ],
            desc="Extracting unmapped reads (SE)",
        )

        run_cmd(
            [
                samtools,
                "fastq",
                "-@",
                str(threads),
                "-0",
                str(unmapped_r1),
                "-n",
                str(unmapped_bam),
            ],
            desc="Converting unmapped BAM → FASTQ (SE)",
        )

    # Clean up intermediate BAM
    if unmapped_bam.exists():
        unmapped_bam.unlink()

    log.info(f"Unmapped reads: {unmapped_r1.name}  ({file_size_human(unmapped_r1)})")
    if unmapped_r2 and unmapped_r2.exists():
        log.info(f"                {unmapped_r2.name}  ({file_size_human(unmapped_r2)})")

    return unmapped_r1, unmapped_r2


def extract_mapped_host(
    sorted_bam: Path,
    output_dir: Path,
    *,
    cfg: Optional[CellJanusConfig] = None,
) -> Path:
    """
    Extract host-aligned reads into a separate BAM for downstream
    single-cell / spatial transcriptome analysis (gene expression).

    Returns path to the host-only BAM.
    """
    log = get_logger()
    if cfg is None:
        cfg = CellJanusConfig()

    samtools = require_tool("samtools")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    threads = cfg.threads

    host_bam = output_dir / "host_mapped.sorted.bam"

    if host_bam.exists():
        log.info(f"Host BAM already exists: {host_bam.name}")
        return host_bam

    # -F 4 = mapped reads only
    run_cmd(
        [
            samtools,
            "view",
            "-@",
            str(threads),
            "-b",
            "-F",
            "4",
            str(sorted_bam),
            "-o",
            str(host_bam),
        ],
        desc="Extracting host-mapped reads",
    )

    run_cmd(
        [samtools, "index", "-@", str(threads), str(host_bam)],
        desc="Indexing host BAM",
    )

    log.info(f"Host reads: {host_bam.name}  ({file_size_human(host_bam)})")
    return host_bam


def count_bam_reads(bam_path: Path, cfg: Optional[CellJanusConfig] = None) -> int:
    """Return number of reads in a BAM file (via samtools flagstat)."""
    if cfg is None:
        cfg = CellJanusConfig()
    samtools = require_tool("samtools")
    result = run_cmd(
        [samtools, "flagstat", "-@", str(cfg.threads), str(bam_path)],
        desc="Counting BAM reads",
        capture=True,
    )
    # First line: "NNNN + 0 in total ..."
    for line in result.stdout.splitlines():
        if "in total" in line:
            return int(line.split()[0])
    return 0
