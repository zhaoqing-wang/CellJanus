"""
CellJanus: A Dual-Perspective Tool for Deconvolving Host Single-Cell
and Microbial Transcriptomes.

Pipeline: FASTQ → QC (fastp) → Align hg38 (Bowtie2/STAR) → Extract unmapped
(samtools) → Classify (Kraken2) → Quantify (Bracken) → Visualize

Supports:
  - Bulk RNA-seq (standard pipeline)
  - scRNA-seq with per-cell barcode tracking (10x Genomics, Parse Bio)
"""

__version__ = "0.2.1"
__author__ = "CellJanus Team"
