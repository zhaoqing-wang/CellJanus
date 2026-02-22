#!/usr/bin/env python3
"""
Build minimal reference files for end-to-end testing.

Creates:
  testdata/refs/
  ├── host_genome/
  │   ├── host.fa                       # Tiny FASTA (human gene fragments)
  │   └── host.*.bt2                    # Bowtie2 index files
  └── kraken2_testdb/
      ├── taxo.k2d / hash.k2d / opts.k2d   # Kraken2 DB files
      └── database150mers.kmer_distrib      # Bracken kmer distribution

Requires: bowtie2-build, kraken2-build, bracken-build in PATH.
Run once — the resulting files are committed to testdata/refs/.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REFS = ROOT / "refs"

# ── Human gene fragments (same as reads in generate_testdata.py) ─────
HUMAN_HBB = (
    "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATG"
    "AAGTTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTC"
    "CTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAG"
    "TGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCTCACTGCAGT"
)
HUMAN_ACTB = (
    "ATGGATGATGATATCGCCGCGCTCGTCGTCGACAACGGCTCCGGCATGTGCAAAGCCGGCTTCGCGGG"
    "CGACGATGCCCCCCGGGCCGTCTTCCCCTCCATCGTGGGGCGCCCCAGGCACCAGGGCGTGATGGTGGG"
    "CATGGGTCAGAAGGATTCCTATGTGGGCGACGAGGCCCAGAGCAAGAGAGGCATCCTCACCCTGAAGTA"
)

# ── Microbial 16S fragments ─────────────────────────────────────────
STAPH_16S = (
    "TGGAGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAAC"
    "AGATAAGGAGCTTGCTCCTTTGACGTTAGCGGCGGACGGGTGAGTAACACGTGGATAACCTACCTATA"
    "AGACTGGGATAACTTCGGGAAACCGGAGCTAATACCGGATAATATATTGAACCGCATGGTTCAATAGTG"
    "AAAGACGGTTTTGCTGTCACTTATAGATGGATCCGCGCCGTATTAGCTAGTTGGTAAGGTAACGGCTTA"
)
ECOLI_16S = (
    "AAATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACG"
    "GTAACAGGAAGCAGCTTGCTGCTTCGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCC"
    "TGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGA"
    "CCTTCGGGCCTCTTGCCATCAGATGTGCCCAGATGGGATTAGCTAGTAGGTGGGGTAACGGCTCACCTA"
)
KPNEU_16S = (
    "AGAGTTTGATCCTGGCTCAGATTGAACGCTGGCGGCAGGCTTAACACATGCAAGTCGAGCGGTAGCAC"
    "AGGGAGCTTGCTCCCTGGGTGACGAGCGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAG"
    "GGGGATAACTACTGGAAACGGTAGCTAATACCGCATAATGTCGCAAGACCAAAGAGGGGGACCTTAGGG"
)


def _write_fasta(path: Path, records: dict[str, str]):
    """Write a simple FASTA file."""
    with open(path, "w") as fh:
        for name, seq in records.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i : i + 80] + "\n")


def _run(cmd: list[str], desc: str):
    print(f"  -> {desc}: {' '.join(cmd[:4])}...")
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  FAILED (exit {r.returncode}):")
        print(f"    stdout: {r.stdout[:300]}")
        print(f"    stderr: {r.stderr[:500]}")
        sys.exit(1)
    return r


def build_host_index():
    """Build a tiny Bowtie2 index from host gene fragments."""
    host_dir = REFS / "host_genome"
    host_dir.mkdir(parents=True, exist_ok=True)

    fasta = host_dir / "host.fa"
    _write_fasta(
        fasta,
        {
            "chr11_HBB_fragment": HUMAN_HBB,
            "chr7_ACTB_fragment": HUMAN_ACTB,
        },
    )

    prefix = host_dir / "host"
    _run(
        ["bowtie2-build", "--quiet", str(fasta), str(prefix)],
        "Building Bowtie2 index",
    )
    print(f"  OK Host index: {prefix}")
    return prefix


def build_kraken2_db():
    """
    Build a minimal Kraken2 database from the 16S fragments.

    Uses the official kraken2-build workflow:
      1. Download NCBI taxonomy (~60 MB)
      2. Add library sequences with proper taxid headers
      3. Build the hash table
      4. Build Bracken kmer distribution
    """
    db_dir = REFS / "kraken2_testdb"
    if db_dir.exists():
        shutil.rmtree(db_dir)
    db_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Download NCBI taxonomy
    print("  (this step downloads ~60 MB taxonomy from NCBI)")
    _run(
        ["kraken2-build", "--download-taxonomy", "--db", str(db_dir)],
        "Downloading NCBI taxonomy",
    )

    # Step 2: Create library FASTA with taxid headers
    # Kraken2 expects: >seqid|kraken:taxid|TAXID  description
    lib_fasta = db_dir / "library_seqs.fa"
    with open(lib_fasta, "w") as fh:
        for taxid, name, seq in [
            (1280, "Staphylococcus_aureus_16S", STAPH_16S),
            (562, "Escherichia_coli_16S", ECOLI_16S),
            (573, "Klebsiella_pneumoniae_16S", KPNEU_16S),
        ]:
            fh.write(f">{name}|kraken:taxid|{taxid}  {name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i : i + 80] + "\n")

    # Add to library via official method
    _run(
        ["kraken2-build", "--add-to-library", str(lib_fasta), "--db", str(db_dir)],
        "Adding 16S library sequences",
    )

    # Step 3: Build the Kraken2 DB
    _run(
        ["kraken2-build", "--build", "--db", str(db_dir), "--threads", "2"],
        "Building Kraken2 database",
    )

    # Step 4: Build Bracken kmer distribution
    _run(
        ["bracken-build", "-d", str(db_dir), "-t", "2", "-k", "35", "-l", "150"],
        "Building Bracken kmer distribution",
    )

    # Step 5: Clean up large intermediate files
    for cleanup in ["library", "taxonomy"]:
        p = db_dir / cleanup
        if p.exists():
            shutil.rmtree(p)
    for pattern in ["*.kraken", "library_seqs.fa", "seqid2taxid.map"]:
        for f in db_dir.glob(pattern):
            f.unlink()

    print(f"  OK Kraken2 DB: {db_dir}")
    return db_dir


if __name__ == "__main__":
    print("Building test reference files for CellJanus...")
    print()

    print("[1/2] Host genome index (Bowtie2)")
    host_prefix = build_host_index()
    print()

    print("[2/2] Kraken2 + Bracken database")
    k2_db = build_kraken2_db()
    print()

    print("Done! Reference files are in testdata/refs/")
    print(f"  Host index prefix: {host_prefix}")
    print(f"  Kraken2 DB:        {k2_db}")
    print()
    print("Test the full pipeline with:")
    print("  celljanus run \\")
    print("    --read1 testdata/reads_R1.fastq.gz \\")
    print("    --read2 testdata/reads_R2.fastq.gz \\")
    print(f"    --host-index {host_prefix} \\")
    print(f"    --kraken2-db {k2_db} \\")
    print("    --output-dir test_results")
