"""Benchmarks for set operations between two tables."""


def test_union(benchmark, populated_table, second_table):
    """union() of two k-mer hash sets."""
    benchmark(lambda: populated_table.union(second_table))


def test_intersection(benchmark, populated_table, second_table):
    """intersection() of two k-mer hash sets."""
    benchmark(lambda: populated_table.intersection(second_table))


def test_difference(benchmark, populated_table, second_table):
    """difference() of two k-mer hash sets."""
    benchmark(lambda: populated_table.difference(second_table))


def test_symmetric_difference(benchmark, populated_table, second_table):
    """symmetric_difference() of two k-mer hash sets."""
    benchmark(lambda: populated_table.symmetric_difference(second_table))
