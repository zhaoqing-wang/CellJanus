"""
Test data generator for CellJanus.

Creates minimal synthetic FASTQ files with reads from:
  • Human genome (chr1 fragments) — should be filtered as host
  • Bacterial sequences — should be classified as microbes
  • Low-quality reads — should be trimmed by QC

This uses BioPython to generate realistic-looking test data without
requiring any database downloads.
"""

from __future__ import annotations

import gzip
import random
import string
from pathlib import Path

random.seed(42)

# ---------------------------------------------------------------------------
# Sequence generators
# ---------------------------------------------------------------------------


def _random_seq(length: int) -> str:
    """Generate a random DNA sequence."""
    return "".join(random.choices("ACGT", k=length))


def _random_quality(length: int, min_q: int = 20, max_q: int = 40) -> str:
    """Generate a random Phred+33 quality string."""
    return "".join(chr(random.randint(min_q + 33, max_q + 33)) for _ in range(length))


def _low_quality_string(length: int) -> str:
    """Mostly low-quality scores."""
    return "".join(chr(random.randint(2 + 33, 12 + 33)) for _ in range(length))


# ---------------------------------------------------------------------------
# Known test sequences (fragments for classification testing)
# ---------------------------------------------------------------------------

# Partial 16S rRNA of Staphylococcus aureus (for Kraken2 hit)
STAPH_AUREUS_16S = (
    "TGGAGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAAC"
    "AGATAAGGAGCTTGCTCCTTTGACGTTAGCGGCGGACGGGTGAGTAACACGTGGATAACCTACCTATA"
)

# Partial human β-globin gene (expected host alignment)
HUMAN_BGLOBIN = (
    "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATG"
    "AAGTTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTC"
)

# E. coli rRNA fragment
ECOLI_16S = (
    "AAATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACG"
    "GTAACAGGAAGCAGCTTGCTGCTTCGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGC"
)


def generate_test_fastq(
    output_dir: Path,
    *,
    n_host_reads: int = 500,
    n_microbe_reads: int = 200,
    n_low_quality: int = 50,
    read_length: int = 150,
    paired: bool = True,
) -> dict[str, Path]:
    """
    Generate synthetic FASTQ test files.

    Returns
    -------
    dict with keys: 'read1', 'read2' (if paired).
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    r1_path = output_dir / "test_R1.fastq.gz"
    r2_path = output_dir / "test_R2.fastq.gz" if paired else None

    records_r1: list[str] = []
    records_r2: list[str] = []

    read_id = 0

    # --- Host reads (human-like) ---
    for i in range(n_host_reads):
        read_id += 1
        # Use human β-globin fragments with mutations
        seq = HUMAN_BGLOBIN[:read_length]
        # Introduce ~2% mutation rate
        seq_list = list(seq)
        for j in range(len(seq_list)):
            if random.random() < 0.02:
                seq_list[j] = random.choice("ACGT")
        seq = "".join(seq_list)
        if len(seq) < read_length:
            seq += _random_seq(read_length - len(seq))

        qual = _random_quality(len(seq), min_q=25, max_q=38)
        records_r1.append(f"@HOST_READ_{read_id}/1\n{seq}\n+\n{qual}\n")
        if paired:
            seq2 = _random_seq(read_length)
            qual2 = _random_quality(read_length, min_q=25, max_q=38)
            records_r2.append(f"@HOST_READ_{read_id}/2\n{seq2}\n+\n{qual2}\n")

    # --- Microbial reads ---
    microbe_templates = [STAPH_AUREUS_16S, ECOLI_16S]
    for i in range(n_microbe_reads):
        read_id += 1
        template = random.choice(microbe_templates)
        seq = template[:read_length]
        seq_list = list(seq)
        for j in range(len(seq_list)):
            if random.random() < 0.01:
                seq_list[j] = random.choice("ACGT")
        seq = "".join(seq_list)
        if len(seq) < read_length:
            seq += _random_seq(read_length - len(seq))

        qual = _random_quality(len(seq), min_q=20, max_q=35)
        species = "STAPH" if template == STAPH_AUREUS_16S else "ECOLI"
        records_r1.append(f"@MICROBE_{species}_{read_id}/1\n{seq}\n+\n{qual}\n")
        if paired:
            seq2 = _random_seq(read_length)
            qual2 = _random_quality(read_length, min_q=20, max_q=35)
            records_r2.append(f"@MICROBE_{species}_{read_id}/2\n{seq2}\n+\n{qual2}\n")

    # --- Low quality reads ---
    for i in range(n_low_quality):
        read_id += 1
        seq = _random_seq(read_length)
        qual = _low_quality_string(read_length)
        records_r1.append(f"@LOWQ_READ_{read_id}/1\n{seq}\n+\n{qual}\n")
        if paired:
            seq2 = _random_seq(read_length)
            qual2 = _low_quality_string(read_length)
            records_r2.append(f"@LOWQ_READ_{read_id}/2\n{seq2}\n+\n{qual2}\n")

    # Shuffle to simulate real interleaving
    combined = list(range(len(records_r1)))
    random.shuffle(combined)

    # Write compressed FASTQ
    with gzip.open(r1_path, "wt") as fh:
        for idx in combined:
            fh.write(records_r1[idx])

    result = {"read1": r1_path}

    if paired and r2_path:
        with gzip.open(r2_path, "wt") as fh:
            for idx in combined:
                fh.write(records_r2[idx])
        result["read2"] = r2_path

    return result


# ---------------------------------------------------------------------------
# scRNA-seq test data generation (10x Genomics style)
# ---------------------------------------------------------------------------

# Valid 10x barcodes (subset for testing)
VALID_BARCODES = [
    "AAACCTGAGCGATGAC",
    "AAACCTGCATCATCCC",
    "AAACCTGGTAAATGTG",
    "AAACCTGTCAACACCA",
    "AAACGGGAGTAGCAAT",
    "AAACGGGCACAGGCCT",
    "AAACGGGCATCGGACC",
    "AAACGGGTCAGCGATT",
    "AAACGGGTCCACGTTC",
    "AAACGGGTCCCAAGTA",
]


def _random_umi(length: int = 12) -> str:
    """Generate a random UMI."""
    return "".join(random.choices("ACGT", k=length))


def generate_scrnaseq_fastq(
    output_dir: Path,
    *,
    n_cells: int = 10,
    reads_per_cell: int = 50,
    n_microbe_reads_per_cell: int = 5,
    read_length: int = 91,
) -> dict:
    """
    Generate 10x Genomics-style scRNA-seq FASTQ test files.

    R1 contains the cell barcode (16bp) + UMI (12bp)
    R2 contains the cDNA sequence

    Returns dict with 'read1', 'read2', 'barcodes' keys.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    r1_path = output_dir / "scrna_R1.fastq.gz"
    r2_path = output_dir / "scrna_R2.fastq.gz"

    records_r1: list[str] = []
    records_r2: list[str] = []
    cell_barcodes_used = []

    read_id = 0

    microbe_templates = [STAPH_AUREUS_16S, ECOLI_16S]
    species_names = ["Staphylococcus_aureus", "Escherichia_coli"]

    for cell_idx in range(min(n_cells, len(VALID_BARCODES))):
        cell_barcode = VALID_BARCODES[cell_idx]
        cell_barcodes_used.append(cell_barcode)

        # Host reads for this cell
        n_host = reads_per_cell - n_microbe_reads_per_cell
        for _ in range(n_host):
            read_id += 1
            umi = _random_umi(12)

            # R1: barcode + UMI + filler
            r1_seq = cell_barcode + umi + _random_seq(read_length - 28)
            r1_qual = _random_quality(len(r1_seq), 25, 38)

            # R2: host-like sequence (human gene fragment)
            r2_seq = HUMAN_BGLOBIN[:read_length]
            if len(r2_seq) < read_length:
                r2_seq += _random_seq(read_length - len(r2_seq))
            r2_qual = _random_quality(len(r2_seq), 25, 38)

            # 10x-style header with CB and UB tags
            header = f"@A00123:456:ABCDEFGHI:1:1101:{1000 + read_id}:{2000 + read_id} CB:Z:{cell_barcode} UB:Z:{umi}"

            records_r1.append(f"{header}\n{r1_seq}\n+\n{r1_qual}\n")
            records_r2.append(f"{header}\n{r2_seq}\n+\n{r2_qual}\n")

        # Microbial reads for this cell
        for _ in range(n_microbe_reads_per_cell):
            read_id += 1
            umi = _random_umi(12)

            # Pick random microbe
            tmpl_idx = random.randint(0, len(microbe_templates) - 1)
            template = microbe_templates[tmpl_idx]

            # R1: barcode + UMI
            r1_seq = cell_barcode + umi + _random_seq(read_length - 28)
            r1_qual = _random_quality(len(r1_seq), 20, 35)

            # R2: microbial sequence
            r2_seq = template[:read_length]
            if len(r2_seq) < read_length:
                r2_seq += _random_seq(read_length - len(r2_seq))
            # Add slight mutations
            r2_list = list(r2_seq)
            for j in range(len(r2_list)):
                if random.random() < 0.01:
                    r2_list[j] = random.choice("ACGT")
            r2_seq = "".join(r2_list)
            r2_qual = _random_quality(len(r2_seq), 20, 35)

            header = f"@A00123:456:ABCDEFGHI:1:1101:{1000 + read_id}:{2000 + read_id} CB:Z:{cell_barcode} UB:Z:{umi}"

            records_r1.append(f"{header}\n{r1_seq}\n+\n{r1_qual}\n")
            records_r2.append(f"{header}\n{r2_seq}\n+\n{r2_qual}\n")

    # Shuffle records
    indices = list(range(len(records_r1)))
    random.shuffle(indices)

    with gzip.open(r1_path, "wt") as fh:
        for idx in indices:
            fh.write(records_r1[idx])

    with gzip.open(r2_path, "wt") as fh:
        for idx in indices:
            fh.write(records_r2[idx])

    # Write barcode whitelist
    whitelist_path = output_dir / "barcodes.txt"
    with open(whitelist_path, "w") as fh:
        for bc in cell_barcodes_used:
            fh.write(f"{bc}\n")

    return {
        "read1": r1_path,
        "read2": r2_path,
        "whitelist": whitelist_path,
        "n_cells": len(cell_barcodes_used),
        "total_reads": len(records_r1),
        "barcodes": cell_barcodes_used,
    }


# ---------------------------------------------------------------------------
# CLI entry point for generating test data
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import sys

    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("test_data")
    paths = generate_test_fastq(out)
    print(f"Generated test data:")
    for k, v in paths.items():
        print(f"  {k}: {v}  ({v.stat().st_size} bytes)")

    # Also generate scRNA-seq test data
    scrna_out = out / "scrnaseq"
    scrna_paths = generate_scrnaseq_fastq(scrna_out)
    print(f"\nGenerated scRNA-seq test data:")
    for k, v in scrna_paths.items():
        if isinstance(v, Path):
            print(f"  {k}: {v}  ({v.stat().st_size} bytes)")
        else:
            print(f"  {k}: {v}")
