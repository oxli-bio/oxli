"""Benchmarks for similarity metrics between two tables."""


def test_jaccard(benchmark, populated_table, second_table):
    """jaccard() similarity (intersection / union of hash sets)."""
    benchmark(lambda: populated_table.jaccard(second_table))


def test_cosine(benchmark, populated_table, second_table):
    """cosine() similarity (Rayon-parallel dot product + magnitudes)."""
    benchmark(lambda: populated_table.cosine(second_table))
