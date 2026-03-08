# CellJanus Test Data

Synthetic paired-end and single-end Illumina FASTQ files for quick testing.
All microbial reads are derived from validated NCBI RefSeq genome fragments
that classify at species level with both the minimal `kraken2_testdb` and
the full Kraken2 `standard_8` database.

## Files

### Bulk RNA-seq Test Data

| File | Description | Reads | Size |
|------|-------------|------:|-----:|
| `reads_R1.fastq.gz` | Paired-end R1 | 1,000 | ~126 KB |
| `reads_R2.fastq.gz` | Paired-end R2 | 1,000 | ~127 KB |
| `reads_SE.fastq.gz` | Single-end | 1,000 | ~126 KB |

### scRNA-seq Test Data

| File | Description | Reads | Cells |
|------|-------------|------:|------:|
| `scrnaseq/scrna_R1.fastq.gz` | 10x-style R1 with CB+UMI | 15,000 | 300 |
| `scrnaseq/scrna_R2.fastq.gz` | 10x-style R2 with cDNA | 15,000 | 300 |
| `scrnaseq/barcodes.txt` | Cell barcode whitelist | — | 300 |

## Read Composition

### Bulk Test Data

| Category | Count | Template Source | Expected Behaviour |
|----------|------:|-----------------|-------------------|
| Host (human) | 600 | HBB (β-globin), ACTB (actin) | Aligned by Bowtie2 → host BAM |
| Microbial | 300 | *S. aureus*, *E. coli*, *K. pneumoniae* genome fragments | Unmapped → Kraken2 classification |
| Low-quality | 100 | Random + adapter read-through | Filtered by fastp (Q < 15) |

### scRNA-seq Test Data

| Category | Count | Description |
|----------|------:|-------------|
| Host reads | 12,600 (42/cell × 300 cells) | Human β-globin fragments with cell barcodes |
| Microbial reads | 2,400 (8/cell × 300 cells) | 7 species genome fragments with cell barcodes |

#### Microbial Species (7)

| Species | Taxid | Source Genome |
|---------|------:|---------------|
| *Staphylococcus aureus* | 1280 | NCTC 8325 |
| *Escherichia coli* | 562 | K-12 MG1655 |
| *Klebsiella pneumoniae* | 573 | HS11286 |
| *Pseudomonas aeruginosa* | 287 | PAO1 |
| *Bacillus subtilis* | 1423 | 168 |
| *Salmonella enterica* | 28901 | LT2 |
| *Streptococcus pneumoniae* | 1313 | TIGR4 |

Headers include 10x-style tags: `CB:Z:<16bp_barcode> UB:Z:<12bp_UMI>`

## Format

- Bulk headers: `@HOST_READ_1/1`, `@MICROBE_ECOLI_601/1`, `@LOWQ_READ_951/1`
- scRNA-seq headers: `@A00123:456:ABCDEFGHI:1:1101:x:y CB:Z:<barcode> UB:Z:<umi>`
- Read length: 150 bp (bulk), 91 bp (scRNA-seq R2)
- Quality encoding: Phred+33
- Compression: gzip

## Regenerate

```bash
# Bulk test data (generates reads_R1, reads_R2, reads_SE)
python tests/generate_testdata.py

# All test data (bulk + scRNA-seq)
python tests/generate_test_data.py testdata
```
