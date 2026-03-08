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
# Human genomic DNA fragments extracted from hg38.
# These are contiguous genomic regions (NOT mRNA/cDNA) so they align
# correctly to both the test host genome and real hg38 with Bowtie2.

# TP53 region — chr17:7,674,800-7,675,300 (500bp, GC=57.4%)
HUMAN_TP53 = (
    "GGGCCACTGACAACCACCCTTAACCCCTCCTCCCAGAGACCCCAGTTGCAAACCAGACCTCAGGCGGC"
    "TCATAGGGCACCACCACACTATGTCGAAAAGTGTTTCTGTCATCCAAATACTCCACACGCAAATTTCC"
    "TTCCACTCGGATAAGATGCTGAGGAGGGGCCAGACCTAAGAGCAATCAGTGAGGAATCAGAGGCCTG"
    "GGGACCCTGGGCAACCAGCCCTGTCGTCTCTCCAGCCCCAGCTGCTCACCATCGCTATCTGAGCAGCG"
    "CTCATGGTGGGGGCAGCGCCTCACAACCTCCGTCATGTGCTGTGACTGCTTGTAGATGGCCATGGCGC"
    "GGACGCGGGTGCCGGGCGGGGGTGTGGAATCAACCCACAGCTGCACAGGGCAGGTCTTGGCCAGTTGG"
    "CAAAACATCTTGTTGAGGGCAGGGGAGTACTGTAGGAAGAGGAAGGAGACAGAGTTGAAAGTCAGGGC"
    "ACAAGTGAACAGATAAAGCAACTGG"
)

# MSH2 region — chr2:47,596,500-47,597,000 (500bp, GC=51.0%)
HUMAN_MSH2 = (
    "AAGAGGGCTTGGATTCAAGTGTGAGAGAAATGGCAGCTCCCCATGAGGGAAGGGGCATGTTTTCTCTT"
    "TCTTACCAGCTGTATTCTCCTGATCCAAAATCATTATTCATCATCAAAATAAAATTCCAGAGGGCCAGC"
    "TTGTTCTGGCATGGCAAGGACAGTCAGATGACCTGAGGGGGCATAATGAGGGGAAGGATCATACAGTT"
    "AGTGGGGCTGAGACACTGGTGGCATTGGGACAGATGAGTGGAGATCACGCCTATGCCAAGGAACTCAT"
    "CTTCCCCTCGTGCTGCAGCCTTCAAAAGTATCTCAGCACAAGAAGAGGGAGGGAGAAAGAAAAGCTACA"
    "CTTTCCAGAGTCCTGGAGTGGACATTGGCTATTTCTTCCCACCCAGACTCCATTGCCCTTGAGGAACAC"
    "CTCCTCCCCTTTCTCAGCCCCTCTGGTTTGGATGAGGCAGATCCTGTGATGCTCACCCCCAACCCTTGC"
    "CACAGGGTCAGACCAGGACC"
)

# Validated genome fragments from RefSeq reference genomes.
# Each classifies at species level with Kraken2 standard_8.

# E. coli K-12 MG1655 genome fragment (taxid 562)
ECOLI_GENOME = (
    "TGCCGGTCCAGAGTGCGCCAAAACACGTCCTGTTCTGGAGGGGAGTTTCAGACACTCCGTTTCTTGCCTG"
    "AAAGTCGATCCGCTTTAAAAACAATAAGGGCTGACAGTTGTCAGCCCTTTTTCACGCTAAAAGCGATTAT"
    "TTATTCCCGCCAGATGATATGGCAAAGTTTGTGATCTTTTTCGCGGCATAACAGAATGCGGGCAAAAACA"
    "TCGTTGATTTCACCATCTTCACTGTCCGCCAGACCAATCACCACTTCGGCAAAAAAGTCCGGGTTCAAAT"
    "CGAAATCAACATGTTCCTGCCAGTCTTCCGCCGGATCGAATAACTCCGCGCCGCCGCGCTCTTCAAACTG"
    "AAGATTGAACAGCAGAACATCAGCTGGATCGAGATTGTCGGCAGCCAGTTCGAGAAAAATATCGTAGGCC"
    "TGTTCAAGCGTTTCATCTTCAGTCAGGCGATTGTTTAGATCCATATCCATAGTTACTACCTGTTTAACCT"
    "CTGTTGGCGA"
)

# S. aureus NCTC 8325 genome fragment (taxid 1280)
SAUREUS_GENOME = (
    "TAATTTTTTAAACAAATTTAAGAATACATAGTAATATAACAATCTAAACATAAAAACTTTTAACACAACA"
    "CTTAAACCAATGCTTTAATTTTCAATACGTAGCTATAATTTGTTGTAAAATCAAAAAGGTTAAAATGTTA"
    "ATTTTCAAAAAAAGGCTCAAAATATGTTTGATTTAGTTATTAAATGTTAAGATATATAAGACTACTATTT"
    "CTTTGTAAAAATGAATCCGATTTACGAGTGAGTAATAGTGAAGGCAGTTTTAAGTTGAAGAAGGCAAAAA"
    "GAGTAAATGTTTATTAATATTTGTAGAAACTAGGTAAGCAAATTAGTTGTGAAAATGTTAATGGTTGCGT"
    "GATAATTTCTATATTTAAATTAGTTTGAAGTGAGGGAGAGTATGTCGAATCAAAATTACGACTACAATAA"
    "AAATGAAGATGGAAGTAAGAAGAAAATGAGTACAACAGCGAAAGTAGTTAGCATTGCGACGGTATTGCTA"
    "TTACTCGGAG"
)

# K. pneumoniae HS11286 genome fragment (taxid 573)
KPNEU_GENOME = (
    "CTGTCTTCAGGGGGGGTTAACGGGTAAAAGCGGCGTGTAGTGTCTCCACGCTACGTGCCGCACGCAGCGT"
    "GCGTCCGGGATTGCTCTGGCTGGCCAGCAGGCTCTGCACCTTGCTGACGATATGCGCGCCGATCGGGATA"
    "GCGGAGGTGGCCGCTGGCGATGGCGCATTGCAGGTATGTATCGAACGCGCGGTGGTAACGAAGAGAAAAT"
    "CGTCGATCAGCTTGCCCTGCGGCGACACGGCCTGCGCCCGCACGCCGGCCGGCCACGGGCGTAAATCGCT"
    "CAGCGTCAGGCTCGGGCAGTATTTTTGCACCAGCCGCAGATAGCCGCTTCTGCACAGCGAATTCTTCATT"
    "TCGCCCAGTCCGGAGCGCAGGTTGTTTTGCAGTACCCGGCGGATGCCGGGCGAGGTGAGGATCTCCAGGG"
    "TATCGGCCAGCGAGATATCGCGCTTGCGGTAGCCTTCGCGCTTCAGAGCCAGCACGGCGTTGGGGCCGAC"
    "GGTGACGCTG"
)


def _mutate(seq: str, rate: float = 0.02) -> str:
    """Introduce random SNPs at the given rate."""
    out = list(seq)
    for i in range(len(out)):
        if random.random() < rate:
            out[i] = random.choice("ACGT")
    return "".join(out)


_COMP = str.maketrans("ACGTacgt", "TGCAtgca")

def _revcomp(seq: str) -> str:
    """Reverse complement of a DNA sequence."""
    return seq.translate(_COMP)[::-1]


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
    """Create a single FASTQ record (for single-end)."""
    start = random.randint(0, max(0, len(template) - read_len))
    seq = template[start : start + read_len]
    if len(seq) < read_len:
        seq += "".join(random.choices("ACGT", k=read_len - len(seq)))
    seq = _mutate(seq, mutation_rate)
    qual = _qual_string(len(seq), min_q, max_q)
    header = _instrument_header(read_id, pair)
    return f"{header}\n{seq}\n+\n{qual}\n"


def _make_pe_reads(
    template: str,
    read_len: int,
    read_id: int,
    min_q_r1: int = 20,
    max_q_r1: int = 40,
    min_q_r2: int = 18,
    max_q_r2: int = 38,
    mutation_rate: float = 0.02,
) -> tuple[str, str]:
    """
    Create a proper FR paired-end read pair from a template.

    R1 reads forward from the fragment start.
    R2 reads the reverse complement from the fragment end.
    This produces concordant Bowtie2 alignments.
    """
    # Fragment (insert) size: between read_len and template length
    max_insert = min(len(template), read_len * 3)
    min_insert = read_len
    insert_size = random.randint(min_insert, max(min_insert, max_insert))

    # Fragment start position
    max_start = max(0, len(template) - insert_size)
    frag_start = random.randint(0, max_start)
    frag_end = frag_start + insert_size

    # R1: forward read from fragment start
    r1_seq = template[frag_start : frag_start + read_len]
    if len(r1_seq) < read_len:
        r1_seq += "".join(random.choices("ACGT", k=read_len - len(r1_seq)))
    r1_seq = _mutate(r1_seq, mutation_rate)

    # R2: reverse complement from fragment end
    r2_start = max(frag_start, frag_end - read_len)
    r2_seq = _revcomp(template[r2_start : r2_start + read_len])
    if len(r2_seq) < read_len:
        r2_seq += "".join(random.choices("ACGT", k=read_len - len(r2_seq)))
    r2_seq = _mutate(r2_seq, mutation_rate)

    r1_qual = _qual_string(read_len, min_q_r1, max_q_r1)
    r2_qual = _qual_string(read_len, min_q_r2, max_q_r2)
    r1_header = _instrument_header(read_id, 1)
    r2_header = _instrument_header(read_id, 2)

    r1_rec = f"{r1_header}\n{r1_seq}\n+\n{r1_qual}\n"
    r2_rec = f"{r2_header}\n{r2_seq}\n+\n{r2_qual}\n"
    return r1_rec, r2_rec


def generate(output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)

    r1_records = []
    r2_records = []
    se_records = []

    rid = 0

    # ── Host reads (human — 60% of total) ────────────────────────
    host_templates = [HUMAN_TP53, HUMAN_MSH2]
    for _ in range(600):
        rid += 1
        tmpl = random.choice(host_templates)
        r1_rec, r2_rec = _make_pe_reads(tmpl, 150, rid, 25, 38, 23, 37, 0.02)
        r1_records.append(r1_rec)
        r2_records.append(r2_rec)
        se_records.append(_make_read(tmpl, 150, rid, 1, 25, 38, 0.02))

    # ── Microbial reads (bacteria — 30% of total) ────────────────
    microbe_templates = [SAUREUS_GENOME, ECOLI_GENOME, KPNEU_GENOME]
    for _ in range(300):
        rid += 1
        tmpl = random.choice(microbe_templates)
        r1_rec, r2_rec = _make_pe_reads(tmpl, 150, rid, 20, 35, 18, 34, 0.01)
        r1_records.append(r1_rec)
        r2_records.append(r2_rec)
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
