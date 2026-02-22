"""
Full pipeline orchestrator.

Chains every module together:
  QC → Align → Extract → Classify → Visualise

Provides runtime estimation, progress tracking through Rich, and
optional step-skipping for partial reruns.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from celljanus.align import align_to_host, get_alignment_stats
from celljanus.classify import classify_and_quantify
from celljanus.config import CellJanusConfig
from celljanus.extract import extract_mapped_host
from celljanus.qc import QCReport, run_qc
from celljanus.utils import file_size_human, fmt_elapsed, get_logger
from celljanus.visualize import generate_all_plots, plot_qc_dashboard

console = Console(stderr=True)


@dataclass
class PipelineResult:
    """Container for all pipeline outputs."""

    qc_r1: Optional[Path] = None
    qc_r2: Optional[Path] = None
    qc_report: Optional[QCReport] = None
    host_bam: Optional[Path] = None
    host_mapped_bam: Optional[Path] = None
    unmapped_r1: Optional[Path] = None
    unmapped_r2: Optional[Path] = None
    kraken2_report: Optional[Path] = None
    bracken_path: Optional[Path] = None
    bracken_df: Optional[pd.DataFrame] = None
    classify_summary: str = ""
    plots: list[Path] = field(default_factory=list)
    elapsed_seconds: float = 0.0

    def summary_table(self) -> str:
        tbl = Table(title="CellJanus Pipeline Results", show_lines=True)
        tbl.add_column("Step", style="bold cyan")
        tbl.add_column("Output", style="green")
        tbl.add_column("Details")

        if self.qc_r1:
            tbl.add_row("QC", str(self.qc_r1), self.qc_report.summary() if self.qc_report else "")
        if self.host_bam:
            tbl.add_row("Alignment", str(self.host_bam), file_size_human(self.host_bam))
        if self.unmapped_r1:
            tbl.add_row("Unmapped", str(self.unmapped_r1), file_size_human(self.unmapped_r1))
        if self.bracken_path:
            n = len(self.bracken_df) if self.bracken_df is not None else 0
            tbl.add_row("Classification", str(self.bracken_path), f"{n} species detected")
        if self.plots:
            tbl.add_row(
                "Plots", str(len(self.plots)) + " files", ", ".join(p.name for p in self.plots)
            )
        tbl.add_row("Total time", fmt_elapsed(self.elapsed_seconds), "")

        buf = console.export_text() if hasattr(console, "export_text") else ""
        return f"Pipeline completed in {fmt_elapsed(self.elapsed_seconds)}"


def run_pipeline(
    read1: Path,
    *,
    read2: Optional[Path] = None,
    cfg: Optional[CellJanusConfig] = None,
    skip_qc: bool = False,
    skip_align: bool = False,
    skip_classify: bool = False,
    skip_visualize: bool = False,
) -> PipelineResult:
    """
    Execute the complete CellJanus pipeline.

    Parameters
    ----------
    read1 : Path
        Path to R1 FASTQ (or single-end FASTQ).
    read2 : Path, optional
        Path to R2 FASTQ for paired-end data.
    cfg : CellJanusConfig
        Pipeline configuration.
    skip_* : bool
        Skip individual pipeline steps (for partial reruns).

    Returns
    -------
    PipelineResult with paths to every output.
    """
    log = get_logger()
    if cfg is None:
        cfg = CellJanusConfig()

    cfg.ensure_dirs()
    log.info(f"System: {cfg.system_summary}")
    result = PipelineResult()
    t0 = time.perf_counter()
    paired = read2 is not None

    # ---- Banner ----
    console.print(
        Panel.fit(
            "[bold magenta]CellJanus[/bold magenta]: "
            "Dual-Perspective Host–Microbe Deconvolution\n"
            f"Input: {read1.name}" + (f" + {read2.name}" if read2 else "") + "\n"
            f"Output: {cfg.output_dir}",
            border_style="blue",
        )
    )

    # ================================================================
    # Step 1: Quality Control
    # ================================================================
    if not skip_qc:
        log.info("[bold]Step 1/5: Quality Control (fastp)[/bold]")
        qc_dir = cfg.output_dir / "01_qc"
        qc_r1, qc_r2, qc_report = run_qc(read1, qc_dir, read2=read2, cfg=cfg)
        result.qc_r1, result.qc_r2, result.qc_report = qc_r1, qc_r2, qc_report
        read1_for_align = qc_r1
        read2_for_align = qc_r2
    else:
        log.info("Skipping QC (--skip-qc)")
        read1_for_align = read1
        read2_for_align = read2

    # ================================================================
    # Step 2: Host Genome Alignment + Unmapped Extraction
    # ================================================================
    if not skip_align:
        log.info("[bold]Step 2/5: Host Genome Alignment (Bowtie2)[/bold]")
        if cfg.host_index is None:
            raise ValueError(
                "Host genome index not set. Use --host-index or 'celljanus download hg38' first."
            )
        align_dir = cfg.output_dir / "02_alignment"
        sorted_bam, unmapped_r1, unmapped_r2 = align_to_host(
            read1_for_align,
            cfg.host_index,
            align_dir,
            read2=read2_for_align,
            cfg=cfg,
        )
        result.host_bam = sorted_bam

        # Step 3: Unmapped reads are captured directly by Bowtie2
        # (--un-conc-gz for PE, --un-gz for SE)
        log.info("[bold]Step 3/5: Unmapped Read Extraction[/bold]")
        result.unmapped_r1, result.unmapped_r2 = unmapped_r1, unmapped_r2

        # Count unmapped reads for logging
        if unmapped_r1 and unmapped_r1.exists():
            log.info(f"Unmapped reads (microbial candidates): {file_size_human(unmapped_r1)}")

        # Also produce a host-only BAM (for single-cell gene expression)
        host_mapped = extract_mapped_host(sorted_bam, align_dir, cfg=cfg)
        result.host_mapped_bam = host_mapped
    else:
        log.info("Skipping alignment (--skip-align)")
        # If user skips alignment, they must provide unmapped reads
        result.unmapped_r1 = read1
        result.unmapped_r2 = read2

    # ================================================================
    # Step 4: Microbial Classification
    # ================================================================
    if not skip_classify:
        log.info("[bold]Step 4/5: Microbial Classification (Kraken2 + Bracken)[/bold]")
        if cfg.kraken2_db is None:
            raise ValueError(
                "Kraken2 database not set. Use --kraken2-db or "
                "'celljanus download kraken2-db' first."
            )
        classify_dir = cfg.output_dir / "04_classification"

        # Check if we have unmapped reads to classify
        if result.unmapped_r1 is None or not result.unmapped_r1.exists():
            log.warning("No unmapped reads found — skipping classification.")
        else:
            cls_result = classify_and_quantify(
                result.unmapped_r1,
                cfg.kraken2_db,
                classify_dir,
                read2=result.unmapped_r2,
                cfg=cfg,
            )
            result.kraken2_report = cls_result["kraken2_report_path"]
            result.bracken_path = cls_result["bracken_path"]
            result.bracken_df = cls_result["bracken_df"]
            result.classify_summary = cls_result["summary"]
    else:
        log.info("Skipping classification (--skip-classify)")

    # ================================================================
    # Step 5: Visualisation
    # ================================================================
    if not skip_visualize and result.bracken_df is not None:
        log.info("[bold]Step 5/5: Generating Visualisations[/bold]")
        vis_dir = cfg.output_dir / "05_visualisation"
        plots = generate_all_plots(result.bracken_df, vis_dir, cfg=cfg)
        result.plots = plots

        # QC dashboard
        if result.qc_report is not None:
            align_stats = {}
            stats_file = cfg.output_dir / "02_alignment" / "host_align_stats.txt"
            if stats_file.exists():
                align_stats = get_alignment_stats(stats_file)
            dashboard = plot_qc_dashboard(
                qc_report=result.qc_report.__dict__,
                align_stats=align_stats,
                classify_summary=result.classify_summary,
                output_path=vis_dir / "plots" / f"pipeline_dashboard.{cfg.plot_format}",
                cfg=cfg,
            )
            result.plots.append(dashboard)
    else:
        log.info("Skipping visualisation")

    # ================================================================
    # Done
    # ================================================================
    result.elapsed_seconds = time.perf_counter() - t0
    console.print(
        Panel.fit(
            f"[bold green]Pipeline completed in {fmt_elapsed(result.elapsed_seconds)}[/bold green]\n"
            f"Output directory: {cfg.output_dir}",
            border_style="green",
        )
    )
    return result
