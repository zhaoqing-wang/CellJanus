"""
Command-line interface for CellJanus.

Usage examples
--------------
# Download reference databases
celljanus download hg38   --output-dir ./refs
celljanus download kraken2 --output-dir ./refs --db-name standard_8

# Run full pipeline
celljanus run \\
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
from typing import Optional

import click
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from celljanus import __version__
from celljanus.config import CellJanusConfig, find_tools
from celljanus.utils import get_logger

console = Console(stderr=True)

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
    console.print(BANNER, style="bold magenta")
    tools = find_tools()
    tbl = Table(title="External Tool Availability", show_lines=True)
    tbl.add_column("Tool", style="bold")
    tbl.add_column("Status")
    tbl.add_column("Path")

    all_ok = True
    for name, path in tools.items():
        if path:
            tbl.add_row(name, "[green]✔ Found[/green]", path)
        else:
            tbl.add_row(name, "[red]✘ Missing[/red]", "—")
            all_ok = False

    console.print(tbl)
    if all_ok:
        console.print("[bold green]All tools available![/bold green]")
    else:
        console.print(
            "[yellow]Some tools are missing. Install them and add to PATH.[/yellow]\n"
            "Required: fastp, bowtie2, bowtie2-build, samtools, kraken2, bracken"
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

    get_logger(Path(output_dir) / "download.log")
    result = download_hg38_genome(Path(output_dir), download_index=not skip_index)
    console.print(f"[green]hg38 downloaded to {output_dir}[/green]")
    for k, v in result.items():
        console.print(f"  {k}: {v}")


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

    get_logger(Path(output_dir) / "download.log")
    db_dir = download_kraken2_db(Path(output_dir), db_name=db_name)
    console.print(f"[green]Kraken2 DB '{db_name}' downloaded to {db_dir}[/green]")


@download.command("refseq")
@click.option("--taxon", "-t", required=True, help="Taxon name, e.g. 'bacteria'.")
@click.option("--output-dir", "-o", required=True, type=click.Path(), help="Output directory.")
def download_refseq(taxon: str, output_dir: str):
    """Download RefSeq genomes for a given taxon."""
    from celljanus.download import download_refseq_genomes

    get_logger(Path(output_dir) / "download.log")
    out = download_refseq_genomes(taxon, Path(output_dir))
    console.print(f"[green]RefSeq genomes for '{taxon}' → {out}[/green]")


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
    console.print(Panel(report.summary(), title="QC Report", border_style="blue"))


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
    """Align reads to host genome (hg38) using Bowtie2."""
    from celljanus.align import align_to_host

    cfg = CellJanusConfig(output_dir=output_dir, host_index=Path(host_index))
    if threads:
        cfg.threads = threads
    get_logger(cfg.log_file)

    bam = align_to_host(
        Path(read1),
        Path(host_index),
        Path(output_dir),
        read2=Path(read2) if read2 else None,
        cfg=cfg,
    )
    console.print(f"[green]Alignment complete: {bam}[/green]")


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
    """Extract reads that did NOT align to the host genome."""
    from celljanus.extract import extract_unmapped

    cfg = CellJanusConfig(output_dir=output_dir)
    if threads:
        cfg.threads = threads
    get_logger(cfg.log_file)

    r1, r2 = extract_unmapped(Path(bam), Path(output_dir), paired=paired, cfg=cfg)
    console.print(f"[green]Unmapped reads: {r1}[/green]")


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
    console.print(Panel(result["summary"], title="Classification Results", border_style="green"))


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
    console.print(f"[green]Generated {len(plots)} plots in {output_dir}[/green]")


# ======================================================================
# celljanus run — full pipeline
# ======================================================================


@main.command("run")
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
    """Run the complete CellJanus pipeline."""
    console.print(BANNER, style="bold magenta")

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

    get_logger(cfg.log_file)

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
    for label, p in [
        ("QC'd reads", result.qc_r1),
        ("Host BAM", result.host_bam),
        ("Unmapped FASTQ", result.unmapped_r1),
        ("Kraken2 report", result.kraken2_report),
        ("Bracken output", result.bracken_path),
    ]:
        if p:
            tbl.add_row(label, str(p))
    if result.plots:
        tbl.add_row("Plots", f"{len(result.plots)} files in {output_dir}/05_visualisation/")
    console.print(tbl)


if __name__ == "__main__":
    main()
