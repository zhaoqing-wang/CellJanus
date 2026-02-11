"""
Microbial classification and quantification module.

Wraps **Kraken2** for taxonomic classification and **Bracken** for
Bayesian re-estimation of species-level abundance.  This mirrors
MetaScope's ``metascope_id()`` Bayesian reassignment step but uses
purpose-built C++ tools that are orders of magnitude faster.

Output:
  • Kraken2 report  (``kraken2_report.txt``)
  • Kraken2 output  (``kraken2_output.txt``)
  • Bracken report   (``bracken_species.txt``) — species-level counts
  • Parsed pandas DataFrame with columns:
    [name, taxonomy_id, taxonomy_lvl, kraken_assigned, bracken_estimated, fraction]
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd

from celljanus.config import CellJanusConfig, require_tool
from celljanus.utils import file_size_human, get_logger, run_cmd


def run_kraken2(
    read1: Path,
    kraken2_db: Path,
    output_dir: Path,
    *,
    read2: Optional[Path] = None,
    cfg: Optional[CellJanusConfig] = None,
) -> tuple[Path, Path]:
    """
    Classify reads using Kraken2.

    Returns (report_path, output_path).
    """
    log = get_logger()
    if cfg is None:
        cfg = CellJanusConfig()

    kraken2 = require_tool("kraken2")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    report_path = output_dir / "kraken2_report.txt"
    output_path = output_dir / "kraken2_output.txt"

    if report_path.exists() and output_path.exists():
        log.info("Kraken2 output already exists, skipping classification.")
        return report_path, output_path

    cmd = [
        kraken2,
        "--db",
        str(kraken2_db),
        "--threads",
        str(cfg.threads),
        "--confidence",
        str(cfg.kraken2_confidence),
        "--report",
        str(report_path),
        "--output",
        str(output_path),
        "--use-names",
        "--gzip-compressed",  # input is .fastq.gz
        "--memory-mapping",  # memory-map the DB for lower RAM usage
    ]

    if read2 is not None:
        cmd.extend(["--paired", str(read1), str(read2)])
    else:
        cmd.append(str(read1))

    log.info(
        f"Running Kraken2 on {read1.name} ({file_size_human(read1)}) with DB {kraken2_db.name}"
    )
    run_cmd(cmd, desc="Kraken2 taxonomic classification")
    log.info(f"Kraken2 report: {report_path}")
    return report_path, output_path


def run_bracken(
    kraken2_report: Path,
    kraken2_db: Path,
    output_dir: Path,
    *,
    cfg: Optional[CellJanusConfig] = None,
) -> Path:
    """
    Run Bracken for Bayesian abundance re-estimation at the species level.

    Returns path to the Bracken output file.
    """
    log = get_logger()
    if cfg is None:
        cfg = CellJanusConfig()

    bracken = require_tool("bracken")
    output_dir = Path(output_dir)
    bracken_out = output_dir / f"bracken_{cfg.bracken_level}.txt"

    if bracken_out.exists():
        log.info("Bracken output already exists, skipping.")
        return bracken_out

    cmd = [
        bracken,
        "-d",
        str(kraken2_db),
        "-i",
        str(kraken2_report),
        "-o",
        str(bracken_out),
        "-r",
        str(cfg.bracken_read_length),
        "-l",
        cfg.bracken_level,
        "-t",
        "10",  # min reads threshold
    ]

    run_cmd(cmd, desc=f"Bracken abundance estimation (level={cfg.bracken_level})")
    log.info(f"Bracken output: {bracken_out}")
    return bracken_out


def parse_kraken2_report(report_path: Path) -> pd.DataFrame:
    """
    Parse a Kraken2-style report into a tidy DataFrame.

    Columns: pct, reads_clade, reads_taxon, rank, taxid, name
    """
    rows = []
    with open(report_path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            rows.append(
                {
                    "pct": float(parts[0].strip()),
                    "reads_clade": int(parts[1].strip()),
                    "reads_taxon": int(parts[2].strip()),
                    "rank": parts[3].strip(),
                    "taxid": int(parts[4].strip()),
                    "name": parts[5].strip(),
                }
            )
    df = pd.DataFrame(rows)
    return df


def parse_bracken_output(bracken_path: Path) -> pd.DataFrame:
    """
    Parse a Bracken output table into a tidy DataFrame.

    Columns: name, taxonomy_id, taxonomy_lvl, kraken_assigned,
             added_reads, new_est_reads, fraction_total_reads
    """
    df = pd.read_csv(bracken_path, sep="\t")
    # Standardise column names
    col_map = {
        "name": "name",
        "taxonomy_id": "taxonomy_id",
        "taxonomy_lvl": "taxonomy_lvl",
        "kraken_assigned_reads": "kraken_assigned",
        "added_reads": "added_reads",
        "new_est_reads": "bracken_estimated",
        "fraction_total_reads": "fraction",
    }
    df = df.rename(columns=col_map)
    df = df.sort_values("bracken_estimated", ascending=False).reset_index(drop=True)
    return df


def classify_and_quantify(
    read1: Path,
    kraken2_db: Path,
    output_dir: Path,
    *,
    read2: Optional[Path] = None,
    cfg: Optional[CellJanusConfig] = None,
) -> dict:
    """
    Run the full Kraken2 → Bracken pipeline and return parsed results.

    Returns
    -------
    dict with keys:
      - kraken2_report_path
      - kraken2_output_path
      - bracken_path
      - kraken2_df   (pd.DataFrame)
      - bracken_df   (pd.DataFrame)
      - summary      (str)
    """
    log = get_logger()
    if cfg is None:
        cfg = CellJanusConfig()

    # Step 1: Kraken2
    report_path, output_path = run_kraken2(read1, kraken2_db, output_dir, read2=read2, cfg=cfg)
    kraken2_df = parse_kraken2_report(report_path)

    # Step 2: Bracken
    bracken_path = run_bracken(report_path, kraken2_db, output_dir, cfg=cfg)
    bracken_df = parse_bracken_output(bracken_path)

    # Summary statistics
    total_classified = kraken2_df.loc[kraken2_df["rank"] == "U", "reads_clade"]
    unclassified = int(total_classified.iloc[0]) if len(total_classified) > 0 else 0
    total_reads = kraken2_df["reads_clade"].iloc[0] if len(kraken2_df) > 0 else 0

    n_species = len(bracken_df)
    top5 = bracken_df.head(5)[["name", "bracken_estimated", "fraction"]].to_string(index=False)

    summary = (
        f"Classified reads: {total_reads - unclassified:,} / {total_reads:,}\n"
        f"Unclassified: {unclassified:,}\n"
        f"Species detected: {n_species}\n"
        f"\nTop 5 species:\n{top5}"
    )
    log.info(f"Classification summary:\n{summary}")

    return {
        "kraken2_report_path": report_path,
        "kraken2_output_path": output_path,
        "bracken_path": bracken_path,
        "kraken2_df": kraken2_df,
        "bracken_df": bracken_df,
        "summary": summary,
    }
