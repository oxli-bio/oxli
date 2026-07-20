"""Shared fixtures for the oxli pytest-codspeed benchmark suite.

These benchmarks live outside ``testpaths`` (``src/python/tests``) so the normal
test run never collects them. Run them with::

    pytest src/python/benchmarks --codspeed

Realistic workloads use the committed ``doc/example.fa`` genome fragment
(~347 KB, ~350k k-mers); per-call micro-benchmarks use small deterministic
synthetic sequences so their measured work is stable across runs.
"""

import gzip
import random
from pathlib import Path

import pytest

import oxli

# doc/example.fa lives at the repo root: benchmarks -> python -> src -> root.
DOC_DIR = Path(__file__).resolve().parents[3] / 'doc'
EXAMPLE_FA = DOC_DIR / 'example.fa'

# Default k-mer size for the shared tables. 21 matches the documented examples.
KSIZE = 21


def read_fasta_seq(path):
    """Concatenate every sequence line of a FASTA file into a single string.

    Parameters
    ----------
    path : pathlib.Path or str
        Path to a plain (uncompressed) FASTA file.

    Returns
    -------
    str
        All non-header lines concatenated, with newlines stripped.
    """
    parts = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith('>'):
                parts.append(line.strip())
    return ''.join(parts)


def random_dna(length, seed=42):
    """Generate a deterministic random DNA sequence.

    Parameters
    ----------
    length : int
        Number of bases to generate.
    seed : int, optional
        Seed for the pseudo-random generator so benchmark work is reproducible.

    Returns
    -------
    str
        A DNA string of the requested length over the alphabet ``ACGT``.
    """
    rng = random.Random(seed)
    return ''.join(rng.choice('ACGT') for _ in range(length))


@pytest.fixture(scope='session')
def example_seq():
    """The concatenated sequence of ``doc/example.fa`` (skips if absent)."""
    if not EXAMPLE_FA.exists():
        pytest.skip(f'benchmark data not found: {EXAMPLE_FA}')
    return read_fasta_seq(EXAMPLE_FA)


@pytest.fixture(scope='session')
def example_fa_path():
    """Filesystem path to ``doc/example.fa`` (skips if absent)."""
    if not EXAMPLE_FA.exists():
        pytest.skip(f'benchmark data not found: {EXAMPLE_FA}')
    return str(EXAMPLE_FA)


@pytest.fixture(scope='session')
def example_gz_path(tmp_path_factory, example_seq):
    """A gzip-compressed copy of the example genome, written once per session."""
    data_dir = tmp_path_factory.mktemp('bench_data')
    gz_path = data_dir / 'example.fa.gz'
    with gzip.open(gz_path, 'wt') as fh:
        fh.write('>example\n')
        fh.write(example_seq)
        fh.write('\n')
    return str(gz_path)


@pytest.fixture(scope='session')
def populated_table(example_seq):
    """A read-only k=21 table populated from the example genome.

    Shared across read-only benchmarks (retrieval, set ops, metrics,
    serialization, histograms). Benchmarks MUST NOT mutate this table.
    """
    table = oxli.KmerCountTable(KSIZE)
    table.consume(example_seq)
    return table


@pytest.fixture(scope='session')
def second_table():
    """A second read-only k=21 table from independent synthetic DNA.

    Provides partial overlap with ``populated_table`` for set operations and
    similarity metrics.
    """
    table = oxli.KmerCountTable(KSIZE)
    table.consume(random_dna(200_000, seed=7))
    return table


@pytest.fixture(scope='session')
def stored_table(example_seq):
    """A read-only k=21 table built with ``store_kmers=True``.

    Used by ``unhash`` and ``dump_kmers`` benchmarks that require the reverse
    hash-to-k-mer mapping.
    """
    table = oxli.KmerCountTable(KSIZE, store_kmers=True)
    table.consume(example_seq)
    return table


@pytest.fixture(scope='session')
def sample_hashes(populated_table):
    """A list of real hash keys drawn from ``populated_table``."""
    return list(populated_table.hashes)


@pytest.fixture(scope='session')
def count_kmers():
    """A fixed list of 2000 distinct k-mers for ``count``/``get`` benchmarks."""
    return [random_dna(KSIZE, seed=i) for i in range(2000)]


@pytest.fixture(scope='session')
def mutation_seq():
    """A ~50k-base synthetic sequence for rebuilding tables in mutation benchmarks.

    Small enough that rebuilding a fresh table each call is cheap relative to the
    mutating operation under test, but large enough to be representative.
    """
    return random_dna(50_000, seed=3)
