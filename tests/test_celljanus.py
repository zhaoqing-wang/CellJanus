"""Tests for CellJanus core modules."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from celljanus.config import CellJanusConfig, default_threads, find_tools, total_memory_gb


# ──────────────────────────────────────────────────────────────────────
# Config tests
# ──────────────────────────────────────────────────────────────────────


class TestConfig:
    def test_default_config(self):
        cfg = CellJanusConfig()
        assert cfg.threads >= 1
        assert cfg.max_memory_gb > 0
        assert cfg.qc_min_quality == 15
        assert cfg.qc_min_length == 36

    def test_ensure_dirs(self, tmp_path):
        cfg = CellJanusConfig(output_dir=tmp_path / "out")
        cfg.ensure_dirs()
        assert cfg.output_dir.exists()
        assert cfg.temp_dir.exists()

    def test_system_summary(self):
        cfg = CellJanusConfig()
        s = cfg.system_summary
        assert "CPUs=" in s
        assert "RAM=" in s

    def test_find_tools_returns_dict(self):
        tools = find_tools()
        assert isinstance(tools, dict)
        assert "fastp" in tools
        assert "bowtie2" in tools
        assert "samtools" in tools
        assert "kraken2" in tools

    def test_default_threads(self):
        t = default_threads()
        assert t >= 1

    def test_total_memory(self):
        mem = total_memory_gb()
        assert mem > 0


# ──────────────────────────────────────────────────────────────────────
# Utility tests
# ──────────────────────────────────────────────────────────────────────


class TestUtils:
    def test_fmt_elapsed(self):
        from celljanus.utils import fmt_elapsed

        assert fmt_elapsed(65) == "1:05"
        assert fmt_elapsed(3661) == "1:01:01"
        assert fmt_elapsed(5) == "0:05"

    def test_file_size_human(self, tmp_path):
        from celljanus.utils import file_size_human

        f = tmp_path / "test.txt"
        f.write_text("x" * 1024)
        s = file_size_human(f)
        assert "KB" in s or "B" in s

    def test_ensure_parent(self, tmp_path):
        from celljanus.utils import ensure_parent

        p = tmp_path / "a" / "b" / "c.txt"
        ensure_parent(p)
        assert p.parent.exists()


# ──────────────────────────────────────────────────────────────────────
# Test data generation tests
# ──────────────────────────────────────────────────────────────────────


class TestTestData:
    def test_generate_test_fastq(self, tmp_path):
        from tests.generate_test_data import generate_test_fastq

        result = generate_test_fastq(
            tmp_path,
            n_host_reads=100,
            n_microbe_reads=50,
            n_low_quality=10,
        )
        assert result["read1"].exists()
        assert result["read2"].exists()

        # Check file is valid gzipped FASTQ
        with gzip.open(result["read1"], "rt") as fh:
            lines = fh.readlines()
        # 160 reads × 4 lines each
        assert len(lines) == (100 + 50 + 10) * 4

    def test_generate_single_end(self, tmp_path):
        from tests.generate_test_data import generate_test_fastq

        result = generate_test_fastq(
            tmp_path / "se", paired=False, n_host_reads=20, n_microbe_reads=10, n_low_quality=5
        )
        assert result["read1"].exists()
        assert "read2" not in result


# ──────────────────────────────────────────────────────────────────────
# QC module tests (require fastp)
# ──────────────────────────────────────────────────────────────────────


class TestQC:
    @pytest.mark.skipif(find_tools().get("fastp") is None, reason="fastp not installed")
    def test_run_qc_single_end(self, test_read1, tmp_path):
        from celljanus.qc import run_qc

        out_r1, out_r2, report = run_qc(
            test_read1,
            tmp_path / "qc_out",
            cfg=CellJanusConfig(output_dir=tmp_path),
        )
        assert out_r1.exists()
        assert out_r2 is None
        assert report.total_reads_before > 0
        assert report.total_reads_after > 0
        assert report.reads_retained_pct > 0

    @pytest.mark.skipif(find_tools().get("fastp") is None, reason="fastp not installed")
    def test_run_qc_paired_end(self, test_read1, test_read2, tmp_path):
        from celljanus.qc import run_qc

        out_r1, out_r2, report = run_qc(
            test_read1,
            tmp_path / "qc_pe",
            read2=test_read2,
            cfg=CellJanusConfig(output_dir=tmp_path),
        )
        assert out_r1.exists()
        assert out_r2 is not None
        assert out_r2.exists()


# ──────────────────────────────────────────────────────────────────────
# CLI tests
# ──────────────────────────────────────────────────────────────────────


class TestCLI:
    def test_cli_help(self):
        from click.testing import CliRunner
        from celljanus.cli import main

        runner = CliRunner()
        result = runner.invoke(main, ["--help"])
        assert result.exit_code == 0
        assert "CellJanus" in result.output

    def test_cli_version(self):
        from click.testing import CliRunner
        from celljanus.cli import main

        runner = CliRunner()
        result = runner.invoke(main, ["--version"])
        assert result.exit_code == 0
        assert "0.1.5" in result.output

    def test_cli_check(self):
        from click.testing import CliRunner
        from celljanus.cli import main

        runner = CliRunner()
        result = runner.invoke(main, ["check"])
        assert result.exit_code == 0

    def test_cli_qc_help(self):
        from click.testing import CliRunner
        from celljanus.cli import main

        runner = CliRunner()
        result = runner.invoke(main, ["qc", "--help"])
        assert result.exit_code == 0
        assert "--read1" in result.output

    def test_cli_run_help(self):
        from click.testing import CliRunner
        from celljanus.cli import main

        runner = CliRunner()
        result = runner.invoke(main, ["run", "--help"])
        assert result.exit_code == 0
        assert "--host-index" in result.output
        assert "--kraken2-db" in result.output

    def test_cli_scrnaseq_help(self):
        from click.testing import CliRunner
        from celljanus.cli import main

        runner = CliRunner()
        result = runner.invoke(main, ["scrnaseq", "--help"])
        assert result.exit_code == 0
        assert "--barcode-mode" in result.output
        assert "--kraken2-db" in result.output


# ──────────────────────────────────────────────────────────────────────
# scRNA-seq module tests
# ──────────────────────────────────────────────────────────────────────


class TestScRNASeq:
    def test_extract_barcode_10x(self):
        from celljanus.scrnaseq import extract_barcode_10x

        # Test 10x-style header
        header = "@A00123:456:ABCDEF:1:1101:1000:2000 CB:Z:AAACCTGAGCGATGAC UB:Z:ACTGACTGACTG"
        cb, umi = extract_barcode_10x(header)
        assert cb == "AAACCTGAGCGATGAC"
        assert umi == "ACTGACTGACTG"

    def test_extract_barcode_missing(self):
        from celljanus.scrnaseq import extract_barcode_10x

        header = "@read1/1"
        cb, umi = extract_barcode_10x(header)
        assert cb is None
        assert umi is None

    def test_barcode_config_defaults(self):
        from celljanus.scrnaseq import BarcodeConfig

        cfg = BarcodeConfig()
        assert cfg.mode == "10x"
        assert cfg.cb_length == 16
        assert cfg.umi_length == 12
        assert cfg.min_reads_per_cell == 1

    def test_cell_microbial_abundance(self):
        from celljanus.scrnaseq import CellMicrobialAbundance

        abundance = CellMicrobialAbundance()

        # Add some classifications
        abundance.add_classification("CELL1", "E. coli", umi="UMI1")
        abundance.add_classification("CELL1", "E. coli", umi="UMI2")
        abundance.add_classification("CELL1", "S. aureus", umi="UMI3")
        abundance.add_classification("CELL2", "E. coli", umi="UMI4")

        # Test summary
        summary = abundance.summary()
        assert summary["total_cells"] == 2
        assert summary["species_detected"] == 2
        assert summary["total_microbial_reads"] == 4

    def test_cell_abundance_to_matrix(self):
        from celljanus.scrnaseq import CellMicrobialAbundance

        abundance = CellMicrobialAbundance()
        abundance.add_classification("CELL1", "E. coli")
        abundance.add_classification("CELL1", "E. coli")
        abundance.add_classification("CELL2", "S. aureus")

        matrix = abundance.to_matrix(min_reads=1)
        assert len(matrix) == 2
        assert "E. coli" in matrix.columns or "S. aureus" in matrix.columns

    def test_wsl2_detection(self):
        from celljanus.scrnaseq import detect_wsl2, is_cross_filesystem_path
        import platform

        # These are functional tests - results depend on environment
        wsl2 = detect_wsl2()
        assert isinstance(wsl2, bool)

        # Test cross-filesystem detection
        # On native Windows, these functions return False since they're designed for WSL2
        from pathlib import Path

        if platform.system() == "Windows":
            # On Windows, is_cross_filesystem_path always returns False
            assert is_cross_filesystem_path(Path("/mnt/c/Users/test")) == False
            assert is_cross_filesystem_path(Path("/home/user/data")) == False
        else:
            # On Linux/WSL2, check actual path patterns
            assert is_cross_filesystem_path(Path("/mnt/c/Users/test")) == True
            assert is_cross_filesystem_path(Path("/home/user/data")) == False

    def test_generate_scrnaseq_testdata(self, tmp_path):
        from tests.generate_test_data import generate_scrnaseq_fastq

        result = generate_scrnaseq_fastq(
            tmp_path / "scrnaseq",
            n_cells=5,
            reads_per_cell=20,
            n_microbe_reads_per_cell=3,
        )

        assert result["read1"].exists()
        assert result["read2"].exists()
        assert result["whitelist"].exists()
        assert result["n_cells"] == 5
        assert result["total_reads"] == 5 * 20
