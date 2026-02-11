"""Pytest fixtures for CellJanus tests."""

from __future__ import annotations

from pathlib import Path

import pytest


@pytest.fixture(scope="session")
def test_data_dir(tmp_path_factory) -> Path:
    """Generate test FASTQ files once per session."""
    from tests.generate_test_data import generate_test_fastq

    d = tmp_path_factory.mktemp("celljanus_test")
    generate_test_fastq(d)
    return d


@pytest.fixture
def test_read1(test_data_dir) -> Path:
    return test_data_dir / "test_R1.fastq.gz"


@pytest.fixture
def test_read2(test_data_dir) -> Path:
    return test_data_dir / "test_R2.fastq.gz"
