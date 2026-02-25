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

# Klebsiella pneumoniae 16S rRNA fragment
KPNEU_16S = (
    "AGAGTTTGATCCTGGCTCAGATTGAACGCTGGCGGCAGGCTTAACACATGCAAGTCGAGCGGTAGCAC"
    "AGGGAGCTTGCTCCCTGGGTGACGAGCGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAG"
)

# Prevotella 16S rRNA fragment (Prevotellaceae representative)
PREVO_16S = (
    "AGAGTTTGATCCTGGCTCAGGATGAACGCTAGCTACAGGCTTAACACATGCAAGTCGAGGGGCAGCAT"
    "GAGTTTAGCTTGCTAAGGCTGATGGCGACCGGCGCACGGGTGAGTAACGCGTATGCAACCTACCTTCG"
)

# Acetitomaculum 16S rRNA fragment
ACETI_16S = (
    "CAGAGTTTGATCCTGGCTCAGGATGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAACGAAGCAC"
    "TTTAAGGATGCTTGCTTGCCTCGGTGACGAGTGGCGGACGGGTGAGTAACACGTGAGCAATCTGTCCC"
)

# Longispora 16S rRNA fragment
LONGI_16S = (
    "TGGAGTTTGATCCTGGCTCAGGGCGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAGCGATGAAG"
    "CCTTTCGGGGTGGATTAGCGGCGAACGGGTGAGTAACACGTGGGCAATCTGCCCTGCACTCTGGGACA"
)

# Mobiluncus 16S rRNA fragment
MOBIL_16S = (
    "TGGAGTTTGATCCTGGCTCAGAACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAACGGAAAGG"
    "CCCCTTCGGGGGTACTCTTGATGCCGACGAGTGGCGGACGGGTGAGTAACACGTAAGTAACCTGCCCC"
)


def generate_test_fastq(
    output_dir: Path,
    *,
    n_host_reads: int = 600,
    n_microbe_reads: int = 350,
    n_low_quality: int = 50,
    read_length: int = 150,
    paired: bool = True,
) -> dict[str, Path]:
    """
    Generate synthetic FASTQ test files with 7 microbial species.

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

    # --- Microbial reads (7 species with varied abundances) ---
    microbe_templates = [
        (STAPH_AUREUS_16S, "STAPH", 30),  # 30% abundance
        (ECOLI_16S, "ECOLI", 25),  # 25% abundance
        (KPNEU_16S, "KPNEU", 15),  # 15% abundance
        (PREVO_16S, "PREVO", 12),  # 12% abundance
        (ACETI_16S, "ACETI", 8),  # 8% abundance
        (LONGI_16S, "LONGI", 6),  # 6% abundance
        (MOBIL_16S, "MOBIL", 4),  # 4% abundance
    ]
    # Create weighted list for random selection
    weighted_templates = []
    for template, name, weight in microbe_templates:
        weighted_templates.extend([(template, name)] * weight)

    for i in range(n_microbe_reads):
        read_id += 1
        template, species = random.choice(weighted_templates)
        seq = template[:read_length]
        seq_list = list(seq)
        for j in range(len(seq_list)):
            if random.random() < 0.01:
                seq_list[j] = random.choice("ACGT")
        seq = "".join(seq_list)
        if len(seq) < read_length:
            seq += _random_seq(read_length - len(seq))

        qual = _random_quality(len(seq), min_q=20, max_q=35)
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


def _generate_barcodes(n: int = 300) -> list[str]:
    """Generate n unique 16bp barcodes in 10x Genomics style."""
    barcodes = set()
    # Use deterministic seed for reproducibility
    rng = random.Random(12345)
    while len(barcodes) < n:
        # Generate 16bp barcode with realistic nucleotide distribution
        bc = "".join(rng.choices("ACGT", k=16))
        barcodes.add(bc)
    return sorted(list(barcodes))


# Generate 300 valid barcodes for testing
VALID_BARCODES = _generate_barcodes(300)


def _random_umi(length: int = 12) -> str:
    """Generate a random UMI."""
    return "".join(random.choices("ACGT", k=length))


def generate_scrnaseq_fastq(
    output_dir: Path,
    *,
    n_cells: int = 300,
    reads_per_cell: int = 50,
    n_microbe_reads_per_cell: int = 8,
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

    # 7 species with different abundances per cell group
    # Create varied abundance patterns for different cell groups
    microbe_templates = [
        (STAPH_AUREUS_16S, "Staphylococcus_aureus"),
        (ECOLI_16S, "Escherichia_coli"),
        (KPNEU_16S, "Klebsiella_pneumoniae"),
        (PREVO_16S, "Prevotella"),
        (ACETI_16S, "Acetitomaculum"),
        (LONGI_16S, "Longispora"),
        (MOBIL_16S, "Mobiluncus"),
    ]

    # Define abundance profiles for different cell groups
    # Each profile is (weights for each species)
    abundance_profiles = [
        [30, 25, 15, 12, 8, 6, 4],  # Profile A: Staph/Ecoli dominant
        [10, 35, 20, 15, 10, 5, 5],  # Profile B: Ecoli dominant
        [15, 15, 30, 20, 10, 5, 5],  # Profile C: Klebsiella/Prevotella dominant
        [20, 10, 10, 25, 20, 10, 5],  # Profile D: Prevotella/Aceti dominant
        [5, 10, 15, 15, 20, 25, 10],  # Profile E: Longi dominant
    ]

    for cell_idx in range(min(n_cells, len(VALID_BARCODES))):
        cell_barcode = VALID_BARCODES[cell_idx]
        cell_barcodes_used.append(cell_barcode)

        # Select abundance profile based on cell index
        profile = abundance_profiles[cell_idx % len(abundance_profiles)]

        # Create weighted template list for this cell
        weighted_templates = []
        for i, (template, name) in enumerate(microbe_templates):
            weighted_templates.extend([(template, name)] * profile[i])

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

        # Microbial reads for this cell (using weighted templates)
        for _ in range(n_microbe_reads_per_cell):
            read_id += 1
            umi = _random_umi(12)

            # Pick microbe based on weighted abundance profile
            template, species_name = random.choice(weighted_templates)

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
