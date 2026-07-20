"""Benchmarks for hashing and canonicalization."""

import oxli


def test_hash_kmer(benchmark, count_kmers):
    """hash_kmer() over a batch of k-mers (per-call PyO3 + MurmurHash cost)."""
    table = oxli.KmerCountTable(21)

    def run():
        for kmer in count_kmers:
            table.hash_kmer(kmer)

    benchmark(run)


def test_canon(benchmark, count_kmers):
    """canon() over a batch of k-mers (reverse-complement + compare)."""
    table = oxli.KmerCountTable(21)

    def run():
        for kmer in count_kmers:
            table.canon(kmer)

    benchmark(run)


def test_kmers_and_hashes(benchmark, example_seq):
    """kmers_and_hashes() over the whole example genome (allocation-heavy)."""
    table = oxli.KmerCountTable(21)

    def run():
        table.kmers_and_hashes(example_seq)

    benchmark(run)


def test_unhash(benchmark, stored_table, sample_hashes):
    """unhash() over a batch of stored hashes (reverse hash->k-mer lookup)."""
    keys = sample_hashes[:5000]

    def run():
        for h in keys:
            stored_table.unhash(h)

    benchmark(run)
