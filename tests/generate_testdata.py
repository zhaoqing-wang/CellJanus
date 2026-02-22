#!/usr/bin/env python3
"""
Generate realistic test FASTQ files for CellJanus.

Modelled after MetaScope's extdata/ — ships small but real-format
FASTQ files so users and CI can run `celljanus qc` immediately.

Output files:
  testdata/
  ├── reads_R1.fastq.gz          # Paired-end R1 (host + microbe + low-Q)
  ├── reads_R2.fastq.gz          # Paired-end R2
  ├── reads_SE.fastq.gz          # Single-end version
  └── README.md                  # Describes the test data
"""

from __future__ import annotations

import gzip
import random
from pathlib import Path

random.seed(2026)

# ── Reference fragments ──────────────────────────────────────────────
# Human β-globin (HBB) — NM_000518.5, exon 1-2 junction
HUMAN_HBB = (
    "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATG"
    "AAGTTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTC"
    "CTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAG"
    "TGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCTCACTGCAGT"
)

# Human ACTB (actin beta) — NM_001101.5
HUMAN_ACTB = (
    "ATGGATGATGATATCGCCGCGCTCGTCGTCGACAACGGCTCCGGCATGTGCAAAGCCGGCTTCGCGGG"
    "CGACGATGCCCCCCGGGCCGTCTTCCCCTCCATCGTGGGGCGCCCCAGGCACCAGGGCGTGATGGTGGG"
    "CATGGGTCAGAAGGATTCCTATGTGGGCGACGAGGCCCAGAGCAAGAGAGGCATCCTCACCCTGAAGTA"
)

# S. aureus 16S rRNA — NR_037007.1
STAPH_16S = (
    "TGGAGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAAC"
    "AGATAAGGAGCTTGCTCCTTTGACGTTAGCGGCGGACGGGTGAGTAACACGTGGATAACCTACCTATA"
    "AGACTGGGATAACTTCGGGAAACCGGAGCTAATACCGGATAATATATTGAACCGCATGGTTCAATAGTG"
    "AAAGACGGTTTTGCTGTCACTTATAGATGGATCCGCGCCGTATTAGCTAGTTGGTAAGGTAACGGCTTA"
)

# E. coli 16S rRNA — NR_024570.1
ECOLI_16S = (
    "AAATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACG"
    "GTAACAGGAAGCAGCTTGCTGCTTCGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCC"
    "TGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGA"
    "CCTTCGGGCCTCTTGCCATCAGATGTGCCCAGATGGGATTAGCTAGTAGGTGGGGTAACGGCTCACCTA"
)

# K. pneumoniae 16S rRNA — NR_036794.1 fragment
KPNEU_16S = (
    "AGAGTTTGATCCTGGCTCAGATTGAACGCTGGCGGCAGGCTTAACACATGCAAGTCGAGCGGTAGCAC"
    "AGGGAGCTTGCTCCCTGGGTGACGAGCGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAG"
    "GGGGATAACTACTGGAAACGGTAGCTAATACCGCATAATGTCGCAAGACCAAAGAGGGGGACCTTAGGG"
)


def _mutate(seq: str, rate: float = 0.02) -> str:
    """Introduce random SNPs at the given rate."""
    out = list(seq)
    for i in range(len(out)):
        if random.random() < rate:
            out[i] = random.choice("ACGT")
    return "".join(out)


def _qual_string(length: int, min_q: int = 20, max_q: int = 40) -> str:
    """Illumina-style Phred+33 quality string."""
    return "".join(chr(random.randint(min_q + 33, max_q + 33)) for _ in range(length))


def _instrument_header(read_id: int, pair: int = 1) -> str:
    """Realistic Illumina read header."""
    lane = random.randint(1, 4)
    tile = random.randint(1101, 2228)
    x = random.randint(1000, 30000)
    y = random.randint(1000, 30000)
    return f"@E00489:42:HVNCCCCXY:{lane}:{tile}:{x}:{y} {pair}:N:0:ATCACG+TTAGGC"


def _make_read(
    template: str,
    read_len: int,
    read_id: int,
    pair: int,
    min_q: int = 20,
    max_q: int = 40,
    mutation_rate: float = 0.02,
) -> str:
    """Create a single FASTQ record."""
    start = random.randint(0, max(0, len(template) - read_len))
    seq = template[start : start + read_len]
    if len(seq) < read_len:
        # Pad with random bases if template too short from this start
        seq += "".join(random.choices("ACGT", k=read_len - len(seq)))
    seq = _mutate(seq, mutation_rate)
    qual = _qual_string(len(seq), min_q, max_q)
    header = _instrument_header(read_id, pair)
    return f"{header}\n{seq}\n+\n{qual}\n"


def generate(output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)

    r1_records = []
    r2_records = []
    se_records = []

    rid = 0

    # ── Host reads (human — 60% of total) ────────────────────────
    host_templates = [HUMAN_HBB, HUMAN_ACTB]
    for _ in range(600):
        rid += 1
        tmpl = random.choice(host_templates)
        r1_records.append(_make_read(tmpl, 150, rid, 1, 25, 38, 0.02))
        r2_records.append(_make_read(tmpl, 150, rid, 2, 23, 37, 0.02))
        se_records.append(_make_read(tmpl, 150, rid, 1, 25, 38, 0.02))

    # ── Microbial reads (bacteria — 30% of total) ────────────────
    microbe_templates = [STAPH_16S, ECOLI_16S, KPNEU_16S]
    for _ in range(300):
        rid += 1
        tmpl = random.choice(microbe_templates)
        r1_records.append(_make_read(tmpl, 150, rid, 1, 20, 35, 0.01))
        r2_records.append(_make_read(tmpl, 150, rid, 2, 18, 34, 0.01))
        se_records.append(_make_read(tmpl, 150, rid, 1, 20, 35, 0.01))

    # ── Low-quality / adapter reads (10%) ────────────────────────
    ILLUMINA_ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
    for _ in range(100):
        rid += 1
        # Short insert → read-through into adapter
        insert = "".join(random.choices("ACGT", k=random.randint(20, 60)))
        seq = (insert + ILLUMINA_ADAPTER)[:150]
        if len(seq) < 150:
            seq += "".join(random.choices("ACGT", k=150 - len(seq)))
        qual = _qual_string(150, 2, 15)  # very low quality
        hdr = _instrument_header(rid, 1)
        r1_records.append(f"{hdr}\n{seq}\n+\n{qual}\n")
        # R2 is pure low quality
        seq2 = "".join(random.choices("ACGT", k=150))
        qual2 = _qual_string(150, 2, 12)
        hdr2 = _instrument_header(rid, 2)
        r2_records.append(f"{hdr2}\n{seq2}\n+\n{qual2}\n")
        se_records.append(f"{hdr}\n{seq}\n+\n{qual}\n")

    # Shuffle
    indices = list(range(len(r1_records)))
    random.shuffle(indices)

    # Write paired-end
    with gzip.open(output_dir / "reads_R1.fastq.gz", "wt") as f:
        for i in indices:
            f.write(r1_records[i])
    with gzip.open(output_dir / "reads_R2.fastq.gz", "wt") as f:
        for i in indices:
            f.write(r2_records[i])
    # Write single-end
    with gzip.open(output_dir / "reads_SE.fastq.gz", "wt") as f:
        for i in indices:
            f.write(se_records[i])

    print(f"Generated {len(r1_records)} reads (PE + SE) in {output_dir}/")
    for name in ["reads_R1.fastq.gz", "reads_R2.fastq.gz", "reads_SE.fastq.gz"]:
        p = output_dir / name
        print(f"  {name}: {p.stat().st_size:,} bytes")


if __name__ == "__main__":
    generate(Path(__file__).parent.parent / "testdata")
