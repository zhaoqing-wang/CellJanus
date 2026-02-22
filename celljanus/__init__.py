"""
CellJanus: A Dual-Perspective Tool for Deconvolving Host Single-Cell
and Microbial Transcriptomes.

Pipeline: FASTQ → QC (fastp) → Align hg38 (Bowtie2) → Extract unmapped
(samtools) → Classify (Kraken2) → Quantify (Bracken) → Visualize
"""

__version__ = "0.1.1"
__author__ = "CellJanus Team"
