"""Benchmarks for counting individual k-mers / hashes."""

import oxli


def test_count(benchmark, count_kmers):
    """count() a batch of k-mers into a fresh table (hash + increment)."""

    def run():
        table = oxli.KmerCountTable(21)
        for kmer in count_kmers:
            table.count(kmer)

    benchmark(run)


def test_count_hash(benchmark, sample_hashes):
    """count_hash() a batch of pre-hashed values (increment only, no hashing)."""
    keys = sample_hashes[:5000]

    def run():
        table = oxli.KmerCountTable(21)
        for h in keys:
            table.count_hash(h)

    benchmark(run)
