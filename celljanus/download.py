"""
Reference data download module.

Handles downloading and indexing of:
  • Human genome hg38 (UCSC / NCBI)
  • Kraken2 / Bracken microbial databases
  • Custom RefSeq genome sets
"""

from __future__ import annotations

import hashlib
import shutil
import tarfile
from pathlib import Path
from typing import Optional

import requests
from rich.progress import (
    BarColumn,
    DownloadColumn,
    Progress,
    TextColumn,
    TimeRemainingColumn,
    TransferSpeedColumn,
)

from celljanus.config import CellJanusConfig, require_tool
from celljanus.utils import get_logger, run_cmd

# ---------------------------------------------------------------------------
# Public URLs for reference genomes and databases
# ---------------------------------------------------------------------------

HG38_URLS = {
    # Primary assembly FASTA (soft-masked) — ~900 MB compressed
    "fasta": ("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"),
    # Pre-built Bowtie2 index — ~4 GB
    "bowtie2_index": ("https://genome-idx.s3.amazonaws.com/bt2/GRCh38_noalt_as.zip"),
}

KRAKEN2_DB_URLS = {
    # MiniKraken2 v2 8 GB — fits in 32 GB RAM
    "minikraken2_8gb": ("https://genome-idx.s3.amazonaws.com/kraken/minikraken2_v2_8GB_201904.tgz"),
    # Standard-8 PlusPF (smaller, ~8 GB)
    "standard_8": ("https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240605.tar.gz"),
    # PlusPFP-8 (plus protozoa, fungi, plant — ~8 GB)
    "pluspfp_8": ("https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_08gb_20240605.tar.gz"),
}

# ---------------------------------------------------------------------------
# Streaming download with progress
# ---------------------------------------------------------------------------


def _download_file(
    url: str,
    dest: Path,
    *,
    desc: str = "Downloading",
    chunk_size: int = 1024 * 1024,  # 1 MB chunks
) -> Path:
    """Stream-download *url* to *dest* with a Rich progress bar."""
    log = get_logger()
    dest.parent.mkdir(parents=True, exist_ok=True)

    if dest.exists():
        log.info(f"File already exists, skipping download: {dest.name}")
        return dest

    log.info(f"Downloading {desc}: {url}")
    resp = requests.get(url, stream=True, timeout=60)
    resp.raise_for_status()
    total = int(resp.headers.get("content-length", 0))

    with Progress(
        TextColumn("[bold blue]{task.description}"),
        BarColumn(),
        DownloadColumn(),
        TransferSpeedColumn(),
        TimeRemainingColumn(),
    ) as progress:
        task = progress.add_task(desc, total=total or None)
        with open(dest, "wb") as fh:
            for chunk in resp.iter_content(chunk_size=chunk_size):
                fh.write(chunk)
                progress.advance(task, len(chunk))

    log.info(f"Saved to {dest}  ({dest.stat().st_size / 1e6:.1f} MB)")
    return dest


def _extract_archive(archive: Path, dest_dir: Path) -> Path:
    """Extract .tar.gz / .tgz / .zip archive into *dest_dir*."""
    log = get_logger()
    dest_dir.mkdir(parents=True, exist_ok=True)

    if archive.suffix == ".zip" or str(archive).endswith(".zip"):
        import zipfile

        log.info(f"Extracting ZIP: {archive.name} → {dest_dir}")
        with zipfile.ZipFile(archive) as zf:
            zf.extractall(dest_dir)
    elif archive.suffixes[-2:] == [".tar", ".gz"] or archive.suffix in (".tgz",):
        log.info(f"Extracting tar.gz: {archive.name} → {dest_dir}")
        with tarfile.open(archive, "r:gz") as tf:
            tf.extractall(dest_dir)
    else:
        log.warning(f"Unknown archive format: {archive}; skipping extraction")
    return dest_dir


# ---------------------------------------------------------------------------
# High-level download functions
# ---------------------------------------------------------------------------


def download_hg38_genome(
    out_dir: Path,
    *,
    download_index: bool = True,
) -> dict[str, Path]:
    """
    Download the human genome hg38 FASTA and (optionally) a pre-built
    Bowtie2 index.

    Returns a dict with keys 'fasta' and/or 'bowtie2_index_dir'.
    """
    log = get_logger()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    result: dict[str, Path] = {}

    # --- FASTA ---
    fasta_dest = out_dir / "hg38.fa.gz"
    _download_file(HG38_URLS["fasta"], fasta_dest, desc="hg38 genome FASTA")
    result["fasta"] = fasta_dest

    # --- Pre-built Bowtie2 index (recommended — saves hours of indexing) ---
    if download_index:
        idx_zip = out_dir / "GRCh38_noalt_as.zip"
        _download_file(HG38_URLS["bowtie2_index"], idx_zip, desc="hg38 Bowtie2 index")
        idx_dir = out_dir / "bowtie2_index"
        _extract_archive(idx_zip, idx_dir)
        result["bowtie2_index_dir"] = idx_dir
        log.info(f"Bowtie2 index extracted to {idx_dir}")

    return result


def download_kraken2_db(
    out_dir: Path,
    *,
    db_name: str = "standard_8",
) -> Path:
    """
    Download a Kraken2 database suitable for 32 GB RAM laptops.

    Parameters
    ----------
    db_name : str
        One of 'minikraken2_8gb', 'standard_8', 'pluspfp_8'.

    Returns
    -------
    Path to the extracted database directory.
    """
    log = get_logger()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if db_name not in KRAKEN2_DB_URLS:
        raise ValueError(f"Unknown Kraken2 DB '{db_name}'. Choose from: {list(KRAKEN2_DB_URLS)}")

    url = KRAKEN2_DB_URLS[db_name]
    archive_name = url.rsplit("/", 1)[-1]
    archive_path = out_dir / archive_name

    _download_file(url, archive_path, desc=f"Kraken2 DB ({db_name})")

    db_dir = out_dir / db_name
    if not db_dir.exists():
        _extract_archive(archive_path, db_dir)
    log.info(f"Kraken2 database ready at {db_dir}")
    return db_dir


def build_bowtie2_index(
    fasta: Path,
    index_dir: Path,
    index_name: str = "host",
    *,
    threads: int = 4,
) -> Path:
    """
    Build a Bowtie2 index from a FASTA file.

    Returns the index prefix path (for use with ``bowtie2 -x``).
    """
    log = get_logger()
    bt2_build = require_tool("bowtie2-build")
    index_dir = Path(index_dir)
    index_dir.mkdir(parents=True, exist_ok=True)
    prefix = index_dir / index_name

    # Skip if index already built
    expected = prefix.parent / f"{index_name}.1.bt2"
    expected_large = prefix.parent / f"{index_name}.1.bt2l"
    if expected.exists() or expected_large.exists():
        log.info(f"Bowtie2 index already exists: {prefix}")
        return prefix

    run_cmd(
        [bt2_build, "--threads", str(threads), str(fasta), str(prefix)],
        desc="Building Bowtie2 index (this may take a while)",
    )
    log.info(f"Bowtie2 index built at {prefix}")
    return prefix


def download_refseq_genomes(
    taxon: str,
    out_dir: Path,
    *,
    assembly_level: str = "complete_genome",
) -> Path:
    """
    Download RefSeq genomes for a given taxon via NCBI datasets CLI or
    direct FTP.  This is a simplified placeholder — for production use
    consider NCBI datasets v2 CLI.

    Returns the directory containing downloaded FASTA files.
    """
    log = get_logger()
    out_dir = Path(out_dir) / taxon.replace(" ", "_")
    out_dir.mkdir(parents=True, exist_ok=True)

    log.info(
        f"RefSeq download for '{taxon}' → {out_dir}\n"
        f"  For large-scale downloads, use:\n"
        f'    datasets download genome taxon "{taxon}" '
        f"--assembly-level {assembly_level} --include genome\n"
        f"  Then extract into {out_dir}"
    )
    return out_dir
