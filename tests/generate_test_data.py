"""
Test data generator for CellJanus.

Creates minimal synthetic FASTQ files with reads from:
  • Human genome (gene fragments) — should be filtered as host
  • Bacterial genome sequences — should be classified as microbes
  • Low-quality reads — should be trimmed by QC

Uses validated genome fragments from NCBI RefSeq reference genomes.
Each microbial fragment classifies at species level with both the
minimal kraken2_testdb and the full Kraken2 standard_8 database.
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


_COMP = str.maketrans("ACGTacgt", "TGCAtgca")

def _revcomp(seq: str) -> str:
    """Reverse complement of a DNA sequence."""
    return seq.translate(_COMP)[::-1]


def _make_pe_reads(
    template: str,
    read_length: int,
    mutation_rate: float = 0.02,
    min_q: int = 20,
    max_q: int = 38,
) -> tuple[str, str]:
    """
    Generate proper FR paired-end sequences from a template.

    R1 reads forward from the fragment start.
    R2 reads the reverse complement from the fragment end.
    Returns (r1_seq, r2_seq) after mutation.
    """
    max_insert = min(len(template), read_length * 3)
    insert_size = random.randint(read_length, max(read_length, max_insert))

    max_start = max(0, len(template) - insert_size)
    frag_start = random.randint(0, max_start)
    frag_end = frag_start + insert_size

    # R1: forward from fragment start
    r1_seq = template[frag_start : frag_start + read_length]
    if len(r1_seq) < read_length:
        r1_seq += _random_seq(read_length - len(r1_seq))

    # R2: reverse complement from fragment end
    r2_start = max(frag_start, frag_end - read_length)
    r2_seq = _revcomp(template[r2_start : r2_start + read_length])
    if len(r2_seq) < read_length:
        r2_seq += _random_seq(read_length - len(r2_seq))

    # Apply mutations
    def _mutate(seq: str) -> str:
        out = list(seq)
        for j in range(len(out)):
            if random.random() < mutation_rate:
                out[j] = random.choice("ACGT")
        return "".join(out)

    return _mutate(r1_seq), _mutate(r2_seq)


# ---------------------------------------------------------------------------
# Known test sequences — validated genome fragments from RefSeq
# ---------------------------------------------------------------------------
# Each 500bp fragment was validated to classify at species level with
# Kraken2 standard_8 (k2_standard_08gb). The first 91bp (10x R2 read
# length) and first 150bp (bulk read length) both resolve to the
# correct species taxid.

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

# Escherichia coli K-12 MG1655 genome fragment (taxid 562)
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

# Staphylococcus aureus NCTC 8325 genome fragment (taxid 1280)
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

# Klebsiella pneumoniae HS11286 genome fragment (taxid 573)
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

# Pseudomonas aeruginosa PAO1 genome fragment (taxid 287)
PAERUG_GENOME = (
    "TGCCAGGCCGTCGCGCGCGGCGGCCTCGAACGCTTCACCGCCGGGACCGCGGTACAGCGCCAGCGGGCCG"
    "TCGCGCTGGAGCAGCTCGGCTTCGATGAAATTGGTCTCGCCAATGAAGTCGGCGACAAAGCGGTTGCCGG"
    "GGCGCTCGTAGATCTCTTCCGGCGGTCCGACCTGTTGCACCTGGCCCTCGGAGAGCACGGCGATGCGATC"
    "GGACATGGTCAGCGCCTCCTCCTGGTCGTGGGTGACGAAAATGAAGGTGATCCCGGTCTTCGCTTGGATC"
    "GCCTTGAGTTCCTCGCGCATCGCCTGGCGCAGCTTCAGGTCGAGGGCCGAGAGCGGTTCGTCGAGCAGCA"
    "GCACCCGCGGTTGCGGCGCCAGCGCACGGGCCAGTGCCACGCGTTGCTGCTGGCCTCCTGAAAGCTGCCC"
    "GGGACGGCGATCGGCGAAACGTTCCATCTGCACCAGGGCGAGCATCTCGCGGACCCGCTCGGCGATCTCC"
    "GCCTTGCCCA"
)

# Bacillus subtilis 168 genome fragment (taxid 1423)
BSUB_GENOME = (
    "CAACCGTTCCCCTGTCGATCACGGAAAGATCAGATTAACGCGTTTTCCAGACTGCGCTCGACTAAGATAA"
    "CAAATGTCTTTTTTTCTGGAAATTTCTATAATAACTTGTCTGACTAGTTCGACAATTAAAATCTGAATTC"
    "CTTCCAATTGTTTATTAAGTAAAATATTTTTCAGATTTACCATCAGCGCTCTGGTGATTGCAAGCTGCTG"
    "CTGTTGACCGCCTTTTCGTTTTGAACAAACGTTTGATTAAAACAAATAGCCCCCTCTGTTTAGCCGGAGT"
    "CGGCCTAATACAGCGACAATCTGCCCCGGGCGAAGCTAAGGCTGTGCAAAACCATCGATTGATCATAACC"
    "AGACGTTACATTCCTGATGGCAGCCCGTTGTTTTGCTCCGCAACCGGAAGCATGGCATGGAGGCCAGTCT"
    "GAGGCGCCGACCTCCACGATACCGAGAACCCGACTATTGCCATTAACTCCGCATCATGGACGGTCACTTC"
    "GCTGATCGCC"
)

# Salmonella enterica LT2 genome fragment (taxid 28901)
SENTER_GENOME = (
    "TTCAGCATCAGCATTAGCGGTGGTCGTAAAGGCTTTATTGATTGCCTGGTTGCATCACCTGACGGCCAAC"
    "AGTGTGCCATCGAAGTTGACAAACGCACACCTCGCAGCCGCTCGCTGTTGAAACTTAGCGAACTGCCTGA"
    "AGGGATATCTGGTTTCGTTCTTCTCAAAGACGGTAAGCATCCTTTGCGTTATAGCGAAGGCGGCATCGAC"
    "GTTATCCGGGCGACGAAATTTAAGTGAATTGAATTAGAAGGGTGGCTGGCAGCTTTTGGGGAGGCCACCA"
    "GCCATGTGAGGGGGAATCCATGAAAACCACATCACAAAATCATTATCGCATCGATATGGGGGCTGGACAA"
    "TGCTGACCATCACGCCAAATTTTGCCCAGGAACGCGCGCTGAACATGCTGCGCCGCGACTGGAAGTCGCA"
    "TAATACTTTTATGGTGTATGCACCCACTGGCAGCGGCAAAACGGGGTTAGCAGCGTTTATCGTTGACGGT"
    "TTCGTCAGTC"
)

# Streptococcus pneumoniae TIGR4 genome fragment (taxid 1313)
SPNEU_GENOME = (
    "AAAGTCAAAAATGAATGACGAAAAATAACCCTGAGAGAGGCTGGAGCCTCTCTTTTTTGTGCAGTTTAGG"
    "AGCTAAAGGGAACAGAATGGAGAAAATGGAACAAATGTGTTTTCTAATCTGTTAGACTGTATCTAGAAAG"
    "GGGAAAATTATGATTAAAGAATTGTATGAAGAAGTCCAAGGGACTGTGTATAAGTGTAGAAATGAATATT"
    "ACCTTCATTTATGGGAATTGTCGGATTGGGAGCAAGAAGGCATGCTCTGCTTACATGAATTGATTAGTAG"
    "AGAAGAAGGACTGGTAGACGATATTCCACGTTTAAGGAAATATTTCAAGACCAAGTTTCGAAATCGAATT"
    "TTAGACTATATCCGTAAACAGGAAAGTCAGAAGCGTAGATACGATAAAGAACCCTATGAAGAAGTGGGTG"
    "AGATCAGTCATCGTATAAGTGAGGGGGGTCTCTGGCTAGATGATTATTATCTCTTTCATGAAACACTAAG"
    "AGATTATAGA"
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

    # --- Host reads (human-like, proper FR paired-end) ---
    host_templates = [HUMAN_TP53, HUMAN_MSH2]
    for i in range(n_host_reads):
        read_id += 1
        template = random.choice(host_templates)
        r1_seq, r2_seq = _make_pe_reads(template, read_length, 0.02, 25, 38)

        qual = _random_quality(read_length, min_q=25, max_q=38)
        records_r1.append(f"@HOST_READ_{read_id}/1\n{r1_seq}\n+\n{qual}\n")
        if paired:
            qual2 = _random_quality(read_length, min_q=25, max_q=38)
            records_r2.append(f"@HOST_READ_{read_id}/2\n{r2_seq}\n+\n{qual2}\n")

    # --- Microbial reads (7 species with varied abundances) ---
    microbe_templates = [
        (SAUREUS_GENOME, "SAUREUS", 30),   # 30% Staphylococcus aureus
        (ECOLI_GENOME, "ECOLI", 25),       # 25% Escherichia coli
        (KPNEU_GENOME, "KPNEU", 15),      # 15% Klebsiella pneumoniae
        (PAERUG_GENOME, "PAERUG", 12),     # 12% Pseudomonas aeruginosa
        (BSUB_GENOME, "BSUB", 8),         #  8% Bacillus subtilis
        (SENTER_GENOME, "SENTER", 6),     #  6% Salmonella enterica
        (SPNEU_GENOME, "SPNEU", 4),       #  4% Streptococcus pneumoniae
    ]
    # Create weighted list for random selection
    weighted_templates = []
    for template, name, weight in microbe_templates:
        weighted_templates.extend([(template, name)] * weight)

    for i in range(n_microbe_reads):
        read_id += 1
        template, species = random.choice(weighted_templates)
        r1_seq, r2_seq = _make_pe_reads(template, read_length, 0.01, 20, 35)

        qual = _random_quality(read_length, min_q=20, max_q=35)
        records_r1.append(f"@MICROBE_{species}_{read_id}/1\n{r1_seq}\n+\n{qual}\n")
        if paired:
            qual2 = _random_quality(read_length, min_q=20, max_q=35)
            records_r2.append(f"@MICROBE_{species}_{read_id}/2\n{r2_seq}\n+\n{qual2}\n")

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
    # Uses validated genome fragments that classify at species level
    microbe_templates = [
        (SAUREUS_GENOME, "Staphylococcus_aureus"),
        (ECOLI_GENOME, "Escherichia_coli"),
        (KPNEU_GENOME, "Klebsiella_pneumoniae"),
        (PAERUG_GENOME, "Pseudomonas_aeruginosa"),
        (BSUB_GENOME, "Bacillus_subtilis"),
        (SENTER_GENOME, "Salmonella_enterica"),
        (SPNEU_GENOME, "Streptococcus_pneumoniae"),
    ]

    # Define abundance profiles for different cell groups
    # Each profile is (weights for each species)
    abundance_profiles = [
        [30, 25, 15, 12, 8, 6, 4],  # Profile A: S.aureus/E.coli dominant
        [10, 35, 20, 15, 10, 5, 5],  # Profile B: E.coli dominant
        [15, 15, 30, 20, 10, 5, 5],  # Profile C: K.pneumoniae dominant
        [20, 10, 10, 25, 20, 10, 5],  # Profile D: P.aeruginosa/B.subtilis
        [5, 10, 15, 15, 20, 25, 10],  # Profile E: S.enterica dominant
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

            # R2: host-like sequence (human genomic DNA fragment)
            host_template = random.choice([HUMAN_TP53, HUMAN_MSH2])
            start = random.randint(0, max(0, len(host_template) - read_length))
            r2_seq = host_template[start : start + read_length]
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
