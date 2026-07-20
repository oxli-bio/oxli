"""Benchmarks for mutating operations.

Each benchmark rebuilds fresh table state inside the measured callable because
these operations mutate the table (running them repeatedly on a shared table
would not be representative). The rebuild is kept cheap relative to the
operation under test.
"""

import oxli


def test_drop(benchmark, count_kmers):
    """drop() a batch of k-mers from a freshly counted table."""

    def run():
        table = oxli.KmerCountTable(21)
        for kmer in count_kmers:
            table.count(kmer)
        for kmer in count_kmers:
            table.drop(kmer)

    benchmark(run)


def test_drop_hash(benchmark, sample_hashes):
    """drop_hash() a batch of hashes from a freshly built table."""
    keys = sample_hashes[:5000]

    def run():
        table = oxli.KmerCountTable(21)
        for h in keys:
            table.count_hash(h)
        for h in keys:
            table.drop_hash(h)

    benchmark(run)


def test_mincut(benchmark, mutation_seq):
    """mincut() a full-table scan removing low-count k-mers."""

    def run():
        table = oxli.KmerCountTable(21)
        table.consume(mutation_seq)
        table.mincut(2)

    benchmark(run)


def test_maxcut(benchmark, mutation_seq):
    """maxcut() a full-table scan removing high-count k-mers."""

    def run():
        table = oxli.KmerCountTable(21)
        table.consume(mutation_seq)
        table.maxcut(1)

    benchmark(run)


def test_add(benchmark, mutation_seq, second_table):
    """add() merges another table's counts into a freshly built table."""

    def run():
        table = oxli.KmerCountTable(21)
        table.consume(mutation_seq)
        table.add(second_table)

    benchmark(run)
