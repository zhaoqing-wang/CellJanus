# CellJanus Test Data

Synthetic paired-end and single-end Illumina FASTQ files for quick testing.

## Files

### Bulk RNA-seq Test Data

| File | Description | Reads | Size |
|------|-------------|------:|-----:|
| `reads_R1.fastq.gz` | Paired-end R1 | 1,000 | ~127 KB |
| `reads_R2.fastq.gz` | Paired-end R2 | 1,000 | ~128 KB |
| `reads_SE.fastq.gz` | Single-end | 1,000 | ~127 KB |

### scRNA-seq Test Data (v0.1.5+)

| File | Description | Reads | Cells |
|------|-------------|------:|------:|
| `scrnaseq/scrna_R1.fastq.gz` | 10x-style R1 with CB+UMI | 500 | 10 |
| `scrnaseq/scrna_R2.fastq.gz` | 10x-style R2 with cDNA | 500 | 10 |
| `scrnaseq/barcodes.txt` | Cell barcode whitelist | — | 10 |

## Read Composition

### Bulk Test Data

| Category | Count | Template Source | Expected Behaviour |
|----------|------:|-----------------|-------------------|
| Host (human) | 600 | HBB (β-globin), ACTB (actin) | Aligned by Bowtie2 → host BAM |
| Microbial | 300 | *S. aureus*, *E. coli*, *K. pneumoniae* 16S rRNA | Unmapped → Kraken2 classification |
| Low-quality | 100 | Random + adapter read-through | Filtered by fastp (Q < 15) |

### scRNA-seq Test Data

| Category | Count | Description |
|----------|------:|-------------|
| Host reads | 450 (45/cell × 10 cells) | Human gene fragments with cell barcodes |
| Microbial reads | 50 (5/cell × 10 cells) | *S. aureus*, *E. coli* 16S with cell barcodes |

Headers include 10x-style tags: `CB:Z:AAACCTGAGCGATGAC UB:Z:ACTGACTGACTG`

## Format

- Illumina-style headers: `@E00489:42:HVNCCCCXY:lane:tile:x:y pair:N:0:ATCACG+TTAGGC`
- Read length: 150 bp (bulk), 91 bp (scRNA-seq)
- Quality encoding: Phred+33
- Compression: gzip

## Regenerate

```bash
# Bulk test data
python tests/generate_testdata.py

# scRNA-seq test data
python -c "from tests.generate_test_data import generate_scrnaseq_fastq; from pathlib import Path; generate_scrnaseq_fastq(Path('testdata/scrnaseq'))"
```
