# CellJanus Test Data

Synthetic paired-end and single-end Illumina FASTQ files for quick testing.

## Files

| File | Description | Reads | Size |
|------|-------------|------:|-----:|
| `reads_R1.fastq.gz` | Paired-end R1 | 1,000 | ~127 KB |
| `reads_R2.fastq.gz` | Paired-end R2 | 1,000 | ~128 KB |
| `reads_SE.fastq.gz` | Single-end | 1,000 | ~127 KB |

## Read Composition

| Category | Count | Template Source | Expected Behaviour |
|----------|------:|-----------------|-------------------|
| Host (human) | 600 | HBB (β-globin), ACTB (actin) | Aligned by Bowtie2 → host BAM |
| Microbial | 300 | *S. aureus*, *E. coli*, *K. pneumoniae* 16S rRNA | Unmapped → Kraken2 classification |
| Low-quality | 100 | Random + adapter read-through | Filtered by fastp (Q < 15) |

## Format

- Illumina-style headers: `@E00489:42:HVNCCCCXY:lane:tile:x:y pair:N:0:ATCACG+TTAGGC`
- Read length: 150 bp
- Quality encoding: Phred+33
- Compression: gzip

## Regenerate

```bash
python tests/generate_testdata.py
```
