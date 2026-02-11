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
# CLI entry point for generating test data
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import sys

    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("test_data")
    paths = generate_test_fastq(out)
    print(f"Generated test data:")
    for k, v in paths.items():
        print(f"  {k}: {v}  ({v.stat().st_size} bytes)")
