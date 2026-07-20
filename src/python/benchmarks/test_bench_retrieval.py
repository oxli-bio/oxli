"""Benchmarks for count retrieval (single and batch)."""


def test_get(benchmark, populated_table, count_kmers):
    """get() a batch of k-mers (hash + lookup) from a populated table."""

    def run():
        for kmer in count_kmers:
            populated_table.get(kmer)

    benchmark(run)


def test_get_hash(benchmark, populated_table, sample_hashes):
    """get_hash() a batch of hashes via a Python loop."""
    keys = sample_hashes[:5000]

    def run():
        for h in keys:
            populated_table.get_hash(h)

    benchmark(run)


def test_get_hash_array(benchmark, populated_table, sample_hashes):
    """get_hash_array() the same batch in one call (native batch lookup).

    Compare against ``test_get_hash`` to see the batching win over a Python loop.
    """
    keys = sample_hashes[:5000]

    def run():
        populated_table.get_hash_array(keys)

    benchmark(run)
