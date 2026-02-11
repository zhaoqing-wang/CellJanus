"""
Quality control module — wraps **fastp** for ultra-fast FASTQ QC.

Design choices for host+microbe dual analysis:
  • Very gentle quality trimming (Q15) to avoid discarding low-abundance
    microbial reads that may have slightly lower quality.
  • Adapter trimming ON (Illumina universal adapters auto-detected).
  • Low-complexity filtering DISABLED by default so that short but real
    microbial sequences are retained.
  • Poly-G / poly-X tail trimming enabled (common Illumina artefact).
  • JSON + HTML reports are always generated for audit trail.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from celljanus.config import CellJanusConfig, require_tool
from celljanus.utils import file_size_human, get_logger, run_cmd


@dataclass
class QCReport:
    """Parsed summary of a fastp QC run."""

    total_reads_before: int = 0
    total_reads_after: int = 0
    total_bases_before: int = 0
    total_bases_after: int = 0
    q20_rate_before: float = 0.0
    q20_rate_after: float = 0.0
    q30_rate_before: float = 0.0
    q30_rate_after: float = 0.0
    adapter_trimmed_reads: int = 0
    duplication_rate: float = 0.0
    gc_content: float = 0.0

    @property
    def reads_retained_pct(self) -> float:
        if self.total_reads_before == 0:
            return 0.0
        return 100.0 * self.total_reads_after / self.total_reads_before

    def summary(self) -> str:
        return (
            f"Reads: {self.total_reads_before:,} → {self.total_reads_after:,} "
            f"({self.reads_retained_pct:.1f}% retained)\n"
            f"Bases: {self.total_bases_before:,} → {self.total_bases_after:,}\n"
            f"Q20:  {self.q20_rate_before:.2%} → {self.q20_rate_after:.2%}\n"
            f"Q30:  {self.q30_rate_before:.2%} → {self.q30_rate_after:.2%}\n"
            f"Adapters trimmed: {self.adapter_trimmed_reads:,}  "
            f"Dup rate: {self.duplication_rate:.2%}  GC: {self.gc_content:.2%}"
        )


def _parse_fastp_json(json_path: Path) -> QCReport:
    """Extract key metrics from fastp's JSON report."""
    with open(json_path) as fh:
        data = json.load(fh)

    s = data.get("summary", {})
    bf = s.get("before_filtering", {})
    af = s.get("after_filtering", {})
    ar = data.get("adapter_cutting", {})
    dup = data.get("duplication", {})

    return QCReport(
        total_reads_before=bf.get("total_reads", 0),
        total_reads_after=af.get("total_reads", 0),
        total_bases_before=bf.get("total_bases", 0),
        total_bases_after=af.get("total_bases", 0),
        q20_rate_before=bf.get("q20_rate", 0.0),
        q20_rate_after=af.get("q20_rate", 0.0),
        q30_rate_before=bf.get("q30_rate", 0.0),
        q30_rate_after=af.get("q30_rate", 0.0),
        adapter_trimmed_reads=ar.get("adapter_trimmed_reads", 0),
        duplication_rate=dup.get("rate", 0.0),
        gc_content=af.get("gc_content", 0.0),
    )


def run_qc(
    read1: Path,
    output_dir: Path,
    *,
    read2: Optional[Path] = None,
    cfg: Optional[CellJanusConfig] = None,
) -> tuple[Path, Optional[Path], QCReport]:
    """
    Run fastp quality control on single-end or paired-end FASTQ files.

    Parameters
    ----------
    read1 : Path
        Path to the R1 FASTQ file (.fastq or .fastq.gz).
    output_dir : Path
        Directory where trimmed reads and reports are written.
    read2 : Path, optional
        Path to the R2 FASTQ file for paired-end data.
    cfg : CellJanusConfig, optional
        Pipeline configuration; defaults are used if not provided.

    Returns
    -------
    (out_r1, out_r2_or_None, QCReport)
    """
    log = get_logger()
    if cfg is None:
        cfg = CellJanusConfig()

    fastp = require_tool("fastp")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    stem = read1.stem.replace(".fastq", "").replace(".fq", "")
    out_r1 = output_dir / f"{stem}_qc.fastq.gz"
    out_r2: Optional[Path] = None
    json_report = output_dir / f"{stem}_fastp.json"
    html_report = output_dir / f"{stem}_fastp.html"

    cmd: list[str] = [
        fastp,
        "--in1",
        str(read1),
        "--out1",
        str(out_r1),
        # Quality
        "--qualified_quality_phred",
        str(cfg.qc_min_quality),
        "--length_required",
        str(cfg.qc_min_length),
        # Adapters
        *(["--detect_adapter_for_pe"] if read2 else []),
        # Poly-G trimming (NextSeq / NovaSeq artefact)
        "--trim_poly_g",
        # Poly-X trimming
        "--trim_poly_x",
        # Low-complexity filter — DISABLED by default to keep microbe reads
        *(
            [
                "--low_complexity_filter",
                "--complexity_threshold",
                str(int(cfg.qc_complexity_threshold)),
            ]
            if cfg.qc_complexity_threshold > 0
            else []
        ),
        # Threading
        "--thread",
        str(min(cfg.threads, 16)),  # fastp max 16
        # Reports
        "--json",
        str(json_report),
        "--html",
        str(html_report),
    ]

    if read2 is not None:
        stem2 = read2.stem.replace(".fastq", "").replace(".fq", "")
        out_r2 = output_dir / f"{stem2}_qc.fastq.gz"
        cmd.extend(["--in2", str(read2), "--out2", str(out_r2)])

    # Adapter trimming toggle
    if not cfg.qc_adapter_trimming:
        cmd.append("--disable_adapter_trimming")

    # Extra user-supplied flags
    if cfg.qc_extra_args:
        cmd.extend(cfg.qc_extra_args.split())

    log.info(
        f"Running fastp QC on {read1.name}"
        + (f" + {read2.name}" if read2 else "")
        + f"  ({file_size_human(read1)})"
    )
    run_cmd(cmd, desc="fastp quality control")

    report = _parse_fastp_json(json_report)
    log.info(f"QC complete — {report.reads_retained_pct:.1f}% reads retained")
    log.info(report.summary())

    return out_r1, out_r2, report
