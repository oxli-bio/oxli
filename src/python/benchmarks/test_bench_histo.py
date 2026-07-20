"""Benchmarks for histograms and whole-table property scans."""


def test_histo_zero(benchmark, populated_table):
    """histo(zero=True) — fills every frequency from 0 to the max count."""
    benchmark(lambda: populated_table.histo(True))


def test_histo_nonzero(benchmark, populated_table):
    """histo(zero=False) — only observed frequencies."""
    benchmark(lambda: populated_table.histo(False))


def test_hashes_property(benchmark, populated_table):
    """hashes property — collect all hash keys into a list."""
    benchmark(lambda: populated_table.hashes)


def test_sum_counts_property(benchmark, populated_table):
    """sum_counts property — full-table count summation."""
    benchmark(lambda: populated_table.sum_counts)
