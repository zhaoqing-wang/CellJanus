"""
Command-line interface for CellJanus.

Usage examples
--------------
# Download reference databases
celljanus download hg38   --output-dir ./refs
celljanus download kraken2 --output-dir ./refs --db-name standard_8

# Run full bulk pipeline
celljanus bulk \\
    --read1 sample_R1.fastq.gz \\
    --read2 sample_R2.fastq.gz \\
    --host-index ./refs/bowtie2_index/GRCh38_noalt_as \\
    --kraken2-db ./refs/standard_8 \\
    --output-dir ./results \\
    --threads 8

# Run individual steps
celljanus qc   --read1 sample_R1.fastq.gz --output-dir ./results/01_qc
celljanus align --read1 qc_R1.fastq.gz   --host-index ./refs/bt2/hg38 ...
celljanus classify --read1 unmapped.fastq.gz --kraken2-db ./refs/k2 ...
celljanus visualize --bracken-file results.txt --output-dir ./plots

# Check tool availability
celljanus check
"""

from __future__ import annotations

from pathlib import Path

import click
from rich.panel import Panel
from rich.table import Table

from celljanus import __version__
from celljanus.config import CellJanusConfig, find_tools
from celljanus.utils import get_logger, log_renderable

# Banner
BANNER = r"""
   ____     _ _     _
  / ___|___| | |   | | __ _ _ __  _   _ ___
 | |   / _ \ | |_  | |/ _` | '_ \| | | / __|
 | |__|  __/ | | |_| | (_| | | | | |_| \__ \
  \____\___|_|_|\___/ \__,_|_| |_|\__,_|___/
  Dual-Perspective Host–Microbe Deconvolution
"""


# ======================================================================
# Top-level group
# ======================================================================


@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(__version__, prog_name="CellJanus")
def main():
    """CellJanus: Dual-Perspective Host–Microbe Deconvolution Pipeline."""
    pass


# ======================================================================
# celljanus check — verify external tools
# ======================================================================


@main.command()
def check():
    """Check that all required external tools are installed."""
    log = get_logger()
    log.info(BANNER)
    tools = find_tools()
    tbl = Table(title="External Tool Availability", show_lines=True)
    tbl.add_column("Tool", style="bold")
    tbl.add_column("Status")
    tbl.add_column("Path")

    optional = {"bowtie2-build"}
    all_required_ok = True
    for name, path in tools.items():
        is_opt = name in optional
        label = f"{name} (optional)" if is_opt else name
        if path:
            tbl.add_row(label, "[green]✔ Found[/green]", path)
        else:
            if is_opt:
                tbl.add_row(label, "[dim]— Not found[/dim]", "—")
            else:
                tbl.add_row(label, "[red]✘ Missing[/red]", "—")
                all_required_ok = False

    log_renderable(tbl)
    if all_required_ok:
        log.info("[bold green]All tools available![/bold green]")
    else:
        log.info(
            "[yellow]Some required tools are missing. Install them and add to PATH.[/yellow]\n"
            "Required: fastp, bowtie2, samtools, kraken2, bracken\n"
            "Optional: bowtie2-build (for celljanus download hg38)"
        )


# ======================================================================
# celljanus download — download reference data
# ======================================================================


@main.group()
def download():
    """Download reference genomes and databases."""
    pass


@download.command("hg38")
@click.option(
    "--output-dir", "-o", required=True, type=click.Path(), help="Directory for downloaded files."
)
@click.option(
    "--skip-index/--with-index", default=False, help="Skip downloading pre-built Bowtie2 index."
)
def download_hg38(output_dir: str, skip_index: bool):
    """Download the human genome hg38 FASTA and Bowtie2 index."""
    from celljanus.download import download_hg38_genome

    log = get_logger(Path(output_dir) / "download.log")
    result = download_hg38_genome(Path(output_dir), download_index=not skip_index)
    log.info(f"[green]hg38 downloaded to {output_dir}[/green]")
    for k, v in result.items():
        log.info(f"  {k}: {v}")


@download.command("kraken2")
@click.option(
    "--output-dir", "-o", required=True, type=click.Path(), help="Directory for downloaded DB."
)
@click.option(
    "--db-name",
    "-n",
    default="standard_8",
    type=click.Choice(["minikraken2_8gb", "standard_8", "pluspfp_8"]),
    help="Kraken2 database variant.",
)
def download_kraken2(output_dir: str, db_name: str):
    """Download a pre-built Kraken2 database."""
    from celljanus.download import download_kraken2_db

    log = get_logger(Path(output_dir) / "download.log")
    db_dir = download_kraken2_db(Path(output_dir), db_name=db_name)
    log.info(f"[green]Kraken2 DB '{db_name}' downloaded to {db_dir}[/green]")


@download.command("refseq")
@click.option("--taxon", "-t", required=True, help="Taxon name, e.g. 'bacteria'.")
@click.option("--output-dir", "-o", required=True, type=click.Path(), help="Output directory.")
def download_refseq(taxon: str, output_dir: str):
    """Download RefSeq genomes for a given taxon."""
    from celljanus.download import download_refseq_genomes

    log = get_logger(Path(output_dir) / "download.log")
    out = download_refseq_genomes(taxon, Path(output_dir))
    log.info(f"[green]RefSeq genomes for '{taxon}' \u2192 {out}[/green]")


# ======================================================================
# celljanus qc — quality control
# ======================================================================


@main.command("qc")
@click.option("--read1", "-1", required=True, type=click.Path(exists=True), help="R1 FASTQ file.")
@click.option(
    "--read2", "-2", default=None, type=click.Path(exists=True), help="R2 FASTQ file (paired-end)."
)
@click.option("--output-dir", "-o", required=True, type=click.Path(), help="Output directory.")
@click.option(
    "--min-quality", default=15, help="Min Phred quality for trimming (default: 15, gentle)."
)
@click.option("--min-length", default=36, help="Min read length after trimming (default: 36).")
@click.option("--threads", "-t", default=None, type=int, help="Number of threads.")
def qc_cmd(read1, read2, output_dir, min_quality, min_length, threads):
    """Run quality control on FASTQ files using fastp."""
    from celljanus.qc import run_qc

    cfg = CellJanusConfig(
        output_dir=output_dir,
        qc_min_quality=min_quality,
        qc_min_length=min_length,
    )
    if threads:
        cfg.threads = threads
    get_logger(cfg.log_file)

    out_r1, out_r2, report = run_qc(
        Path(read1),
        Path(output_dir),
        read2=Path(read2) if read2 else None,
        cfg=cfg,
    )
    log_renderable(Panel(report.summary(), title="QC Report", border_style="blue"))


# ======================================================================
# celljanus align — host genome alignment
# ======================================================================


@main.command("align")
@click.option("--read1", "-1", required=True, type=click.Path(exists=True), help="QC'd R1 FASTQ.")
@click.option("--read2", "-2", default=None, type=click.Path(exists=True), help="QC'd R2 FASTQ.")
@click.option(
    "--host-index",
    "-x",
    required=True,
    type=click.Path(),
    help="Bowtie2 index prefix for host genome.",
)
@click.option("--output-dir", "-o", required=True, type=click.Path(), help="Output directory.")
@click.option("--threads", "-t", default=None, type=int, help="Number of threads.")
def align_cmd(read1, read2, host_index, output_dir, threads):
    """[Bulk] Align reads to host genome (hg38) using Bowtie2."""
    from celljanus.align import align_to_host

    cfg = CellJanusConfig(output_dir=output_dir, host_index=Path(host_index))
    if threads:
        cfg.threads = threads
    get_logger(cfg.log_file)

    bam, unmapped_r1, unmapped_r2 = align_to_host(
        Path(read1),
        Path(host_index),
        Path(output_dir),
        read2=Path(read2) if read2 else None,
        cfg=cfg,
    )
    log = get_logger()
    log.info(f"[green]Alignment complete: {bam}[/green]")
    if unmapped_r1 and unmapped_r1.exists():
        log.info(f"[green]Unmapped reads: {unmapped_r1}[/green]")


# ======================================================================
# celljanus extract — extract unmapped reads
# ======================================================================


@main.command("extract")
@click.option(
    "--bam",
    "-b",
    required=True,
    type=click.Path(exists=True),
    help="Sorted BAM from alignment step.",
)
@click.option("--output-dir", "-o", required=True, type=click.Path(), help="Output directory.")
@click.option("--paired/--single", default=False, help="Paired-end or single-end.")
@click.option("--threads", "-t", default=None, type=int, help="Number of threads.")
def extract_cmd(bam, output_dir, paired, threads):
    """[Bulk] Extract reads that did NOT align to the host genome."""
    from celljanus.extract import extract_unmapped

    cfg = CellJanusConfig(output_dir=output_dir)
    if threads:
        cfg.threads = threads
    get_logger(cfg.log_file)

    r1, r2 = extract_unmapped(Path(bam), Path(output_dir), paired=paired, cfg=cfg)
    get_logger().info(f"[green]Unmapped reads: {r1}[/green]")


# ======================================================================
# celljanus classify — microbial classification
# ======================================================================


@main.command("classify")
@click.option(
    "--read1", "-1", required=True, type=click.Path(exists=True), help="Unmapped R1 FASTQ."
)
@click.option(
    "--read2", "-2", default=None, type=click.Path(exists=True), help="Unmapped R2 FASTQ."
)
@click.option(
    "--kraken2-db",
    "-d",
    required=True,
    type=click.Path(exists=True),
    help="Kraken2 database directory.",
)
@click.option("--output-dir", "-o", required=True, type=click.Path(), help="Output directory.")
@click.option("--confidence", default=0.05, help="Kraken2 confidence threshold (default: 0.05).")
@click.option(
    "--bracken-level",
    default="S",
    type=click.Choice(["D", "P", "C", "O", "F", "G", "S"]),
    help="Bracken taxonomic level (default: S=species).",
)
@click.option("--threads", "-t", default=None, type=int, help="Number of threads.")
def classify_cmd(read1, read2, kraken2_db, output_dir, confidence, bracken_level, threads):
    """Classify microbial reads using Kraken2 + Bracken."""
    from celljanus.classify import classify_and_quantify

    cfg = CellJanusConfig(
        output_dir=output_dir,
        kraken2_db=Path(kraken2_db),
        kraken2_confidence=confidence,
        bracken_level=bracken_level,
    )
    if threads:
        cfg.threads = threads
    get_logger(cfg.log_file)

    result = classify_and_quantify(
        Path(read1),
        Path(kraken2_db),
        Path(output_dir),
        read2=Path(read2) if read2 else None,
        cfg=cfg,
    )
    log_renderable(Panel(result["summary"], title="Classification Results", border_style="green"))


# ======================================================================
# celljanus visualize — generate plots
# ======================================================================


@main.command("visualize")
@click.option(
    "--bracken-file",
    "-b",
    required=True,
    type=click.Path(exists=True),
    help="Bracken output file (TSV).",
)
@click.option(
    "--output-dir", "-o", required=True, type=click.Path(), help="Output directory for plots."
)
@click.option("--top-n", default=20, help="Number of top species to plot.")
@click.option(
    "--format",
    "fmt",
    default="png",
    type=click.Choice(["png", "pdf", "svg"]),
    help="Plot output format.",
)
def visualize_cmd(bracken_file, output_dir, top_n, fmt):
    """Generate abundance plots from Bracken results."""
    from celljanus.classify import parse_bracken_output
    from celljanus.visualize import generate_all_plots

    cfg = CellJanusConfig(output_dir=output_dir, top_n_species=top_n, plot_format=fmt)
    get_logger(cfg.log_file)

    df = parse_bracken_output(Path(bracken_file))
    plots = generate_all_plots(df, Path(output_dir), cfg=cfg)
    get_logger().info(f"[green]Generated {len(plots)} plots in {output_dir}[/green]")


# ======================================================================
# celljanus scrnaseq — scRNA-seq microbial classification
# ======================================================================


@main.command("scrnaseq")
@click.option(
    "--read1",
    "-1",
    required=True,
    type=click.Path(exists=True),
    help="R1 FASTQ with cell barcodes.",
)
@click.option(
    "--read2",
    "-2",
    default=None,
    type=click.Path(exists=True),
    help="R2 FASTQ with cDNA sequences.",
)
@click.option(
    "--kraken2-db",
    "-d",
    required=True,
    type=click.Path(exists=True),
    help="Kraken2 database directory.",
)
@click.option(
    "--output-dir",
    "-o",
    default="celljanus_scrna_output",
    type=click.Path(),
    help="Output directory.",
)
@click.option(
    "--barcode-mode",
    default="10x",
    type=click.Choice(["10x", "parse", "auto"]),
    help="Barcode extraction mode.",
)
@click.option(
    "--whitelist", "-w", default=None, type=click.Path(exists=True), help="Cell barcode whitelist."
)
@click.option(
    "--remove-host-taxa/--keep-host-taxa",
    default=True,
    help="Remove host/root/non-informative taxa from scRNA-seq output before cell aggregation.",
)
@click.option("--min-reads", default=1, help="Minimum reads per cell to include (default: 1).")
@click.option("--threads", "-t", default=None, type=int, help="Number of threads.")
def scrnaseq_cmd(
    read1,
    read2,
    kraken2_db,
    output_dir,
    barcode_mode,
    whitelist,
    remove_host_taxa,
    min_reads,
    threads,
):
    """Run scRNA-seq per-cell microbial classification (10x / Parse)."""
    import time

    from celljanus.scrnaseq import (
        BarcodeConfig,
        run_scrnaseq_classification,
        detect_wsl2,
        wsl2_io_warning,
    )
    from celljanus.utils import fmt_elapsed
    from celljanus.visualize import generate_scrnaseq_plots
    from rich.panel import Panel

    t0 = time.perf_counter()

    cfg = CellJanusConfig(
        output_dir=Path(output_dir),
        kraken2_db=Path(kraken2_db),
    )
    if threads:
        cfg.threads = threads
    log = get_logger(cfg.log_file)

    log.info(BANNER)
    log.info(f"System: {cfg.system_summary}")
    log_renderable(
        Panel.fit(
            "[bold magenta]CellJanus[/bold magenta]: "
            "scRNA-seq Per-Cell Microbial Classification\n"
            f"Input: {Path(read1).name}" + (f" + {Path(read2).name}" if read2 else "") + "\n"
            f"Output: {output_dir}",
            border_style="blue",
        )
    )

    # Check WSL2 I/O performance
    paths_to_check = [Path(read1)]
    if read2:
        paths_to_check.append(Path(read2))

    if detect_wsl2():
        warning = wsl2_io_warning(paths_to_check)
        if warning:
            log.info(f"[yellow]{warning}[/yellow]")

    # Configure barcode extraction
    barcode_cfg = BarcodeConfig(
        mode=barcode_mode,
        whitelist_path=Path(whitelist) if whitelist else None,
        remove_host_taxa=remove_host_taxa,
        min_reads_per_cell=min_reads,
    )

    # Run scRNA-seq classification
    result = run_scrnaseq_classification(
        Path(read1),
        Path(kraken2_db),
        Path(output_dir),
        read2=Path(read2) if read2 else None,
        barcode_cfg=barcode_cfg,
        cfg=cfg,
    )

    matrix_df = result["matrix_df"]
    if len(matrix_df) > 0:
        vis_dir = Path(output_dir) / "02_visualisation"
        plots = generate_scrnaseq_plots(
            matrix_df,
            vis_dir,
            pipeline_summary=result["summary"],
            cfg=cfg,
        )
        log.info(f"[green]Generated {len(plots)} scRNA-seq plots[/green]")
    else:
        log.info(
            "[yellow]Warning: No cells passed filtering. Try lowering --min-reads.[/yellow]"
        )

    # Summary table
    summary = result["summary"]
    tbl = Table(title="scRNA-seq Microbial Classification Results", show_lines=True)
    tbl.add_column("Metric", style="bold cyan")
    tbl.add_column("Value", style="green")
    tbl.add_row("Input reads", f"{summary.get('input_reads', 0):,}")
    tbl.add_row("Raw barcodes with retained classifications", f"{summary['total_cells_raw']:,}")
    tbl.add_row("Cells With Microbe (passing --min-reads)", f"{summary['total_cells']:,}")
    tbl.add_row("Cells filtered out", f"{summary['cells_filtered_out']:,}")
    tbl.add_row("Min reads per cell (--min-reads)", f"{summary['min_reads_per_cell']}")
    tbl.add_row("Species detected (filtered)", f"{summary['species_detected']}")
    tbl.add_row("Total microbial reads (filtered)", f"{summary['total_microbial_reads']:,}")
    tbl.add_row("Mean reads / cell (filtered)", f"{summary['mean_reads_per_cell']:.1f}")
    log_renderable(tbl)

    elapsed = time.perf_counter() - t0
    log_renderable(
        Panel.fit(
            f"[bold green]Pipeline completed in {fmt_elapsed(elapsed)}[/bold green]\n"
            f"Output directory: {output_dir}",
            border_style="green",
        )
    )


# ======================================================================
# celljanus bulk — full bulk RNA-seq pipeline
# ======================================================================


@main.command("bulk")
@click.option("--read1", "-1", required=True, type=click.Path(exists=True), help="R1 FASTQ file.")
@click.option(
    "--read2", "-2", default=None, type=click.Path(exists=True), help="R2 FASTQ file (paired-end)."
)
@click.option(
    "--host-index", "-x", required=True, type=click.Path(), help="Bowtie2 index prefix (e.g. hg38)."
)
@click.option(
    "--kraken2-db",
    "-d",
    required=True,
    type=click.Path(exists=True),
    help="Kraken2 database directory.",
)
@click.option(
    "--output-dir", "-o", default="celljanus_output", type=click.Path(), help="Output directory."
)
@click.option("--threads", "-t", default=None, type=int, help="Number of threads (default: auto).")
@click.option(
    "--max-memory", "-m", default=None, type=float, help="Max memory in GB (default: 80% of RAM)."
)
@click.option("--min-quality", default=15, help="Min Phred quality for QC (default: 15).")
@click.option("--min-length", default=36, help="Min read length for QC (default: 36).")
@click.option("--confidence", default=0.05, help="Kraken2 confidence threshold.")
@click.option(
    "--bracken-level",
    default="S",
    type=click.Choice(["D", "P", "C", "O", "F", "G", "S"]),
    help="Bracken taxonomic level.",
)
@click.option("--plot-format", default="png", type=click.Choice(["png", "pdf", "svg"]))
@click.option("--skip-qc", is_flag=True, help="Skip QC step.")
@click.option("--skip-classify", is_flag=True, help="Skip classification step.")
@click.option("--skip-visualize", is_flag=True, help="Skip visualisation step.")
def run_cmd(
    read1,
    read2,
    host_index,
    kraken2_db,
    output_dir,
    threads,
    max_memory,
    min_quality,
    min_length,
    confidence,
    bracken_level,
    plot_format,
    skip_qc,
    skip_classify,
    skip_visualize,
):
    """Run the complete CellJanus bulk RNA-seq pipeline (QC → Align → Classify → Visualize)."""
    cfg = CellJanusConfig(
        output_dir=Path(output_dir),
        host_index=Path(host_index),
        kraken2_db=Path(kraken2_db),
        qc_min_quality=min_quality,
        qc_min_length=min_length,
        kraken2_confidence=confidence,
        bracken_level=bracken_level,
        plot_format=plot_format,
    )
    if threads:
        cfg.threads = threads
    if max_memory:
        cfg.max_memory_gb = max_memory

    log = get_logger(cfg.log_file)
    log.info(BANNER)

    from celljanus.pipeline import run_pipeline

    result = run_pipeline(
        Path(read1),
        read2=Path(read2) if read2 else None,
        cfg=cfg,
        skip_qc=skip_qc,
        skip_classify=skip_classify,
        skip_visualize=skip_visualize,
    )

    # Final summary
    tbl = Table(title="Pipeline Results", show_lines=True)
    tbl.add_column("Output", style="bold cyan")
    tbl.add_column("Path", style="green")
    tbl.add_column("Size", style="dim")
    for label, p in [
        ("QC'd reads", result.qc_r1),
        ("Host BAM", result.host_bam),
        ("Unmapped FASTQ", result.unmapped_r1),
        ("Kraken2 report", result.kraken2_report),
        ("Bracken output", result.bracken_path),
    ]:
        if p and Path(p).exists():
            from celljanus.utils import file_size_human

            tbl.add_row(label, str(p), file_size_human(Path(p)))
        elif p:
            tbl.add_row(label, str(p), "")
    if result.plots:
        png_count = sum(1 for p in result.plots if str(p).endswith(".png"))
        pdf_count = sum(1 for p in result.plots if str(p).endswith(".pdf"))
        tbl.add_row(
            "Plots",
            f"{output_dir}/05_visualisation/",  # plots stored directly (no plots/ subfolder)
            f"{png_count} PNG + {pdf_count} PDF",
        )

    # Tables
    tables_dir = Path(output_dir) / "06_tables"
    if tables_dir.exists():
        table_files = sorted(tables_dir.glob("*.csv"))
        tbl.add_row(
            "Result tables",
            str(tables_dir),
            f"{len(table_files)} CSV files",
        )

    log_renderable(tbl)


if __name__ == "__main__":
    main()
