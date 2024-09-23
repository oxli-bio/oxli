import pytest

import oxli


@pytest.fixture
def kmer_count_table():
    """Fixture to create a KmerCountTable with ksize=4."""
    table = oxli.KmerCountTable(ksize=4)
    return table


def test_min_empty_table(kmer_count_table):
    """Test min on an empty KmerCountTable.

    Edge case: When the table is empty, min should return 0.
    """
    assert kmer_count_table.min == 0, "min should return 0 for an empty table"


def test_max_empty_table(kmer_count_table):
    """Test max on an empty KmerCountTable.

    Edge case: When the table is empty, max should return 0.
    """
    assert kmer_count_table.max == 0, "max should return 0 for an empty table"


def test_min_non_empty_table(kmer_count_table):
    """Test min on a non-empty KmerCountTable."""
    kmer_count_table.count("AAAA")  # Adding 1 k-mer
    kmer_count_table.count("TTTT")  # Another k-mer with same hash (canonical k-mer)
    kmer_count_table.consume("CCCCCC")  # Count "CCCC" 3 times

    assert (
        kmer_count_table.min == 2
    ), "min should return the minimum count value, in this case 2"


def test_max_non_empty_table(kmer_count_table):
    """Test max on a non-empty KmerCountTable."""
    kmer_count_table.count("AAAA")  # Adding k-mers
    kmer_count_table.count("TTTT")  # Another k-mer with same hash (canonical k-mer)
    kmer_count_table.count("CCCC")  # Another distinct k-mer

    assert (
        kmer_count_table.max == 2
    ), "max should return the maximum count value, in this case 2"


def test_histo_zero_false_empty_table(kmer_count_table):
    """Test histo(zero=False) on an empty KmerCountTable.

    Edge case: When the table is empty, histo() should return an empty list.
    """
    assert (
        kmer_count_table.histo(zero=False) == []
    ), "histo() should return an empty list for an empty table"


def test_histo_zero_true_empty_table(kmer_count_table):
    """Test histo(zero=True) on an empty KmerCountTable.

    Edge case: When the table is empty, histo() should return [(0, 0)].
    """
    assert kmer_count_table.histo(zero=True) == [
        (0, 0)
    ], "histo(zero=True) should return [(0, 0)] for an empty table"


def test_histo_zero_false_non_empty_table(kmer_count_table):
    """
    Test histo(zero=False) on a non-empty KmerCountTable.
    Only observed frequencies should be included in the histogram.
    """
    kmer_count_table.count("AAAA")  # Add k-mer, counts=1
    kmer_count_table.count("AAAA")  # Add k-mer, counts=2
    kmer_count_table.count("TTTT")  # Add another k-mer, canonical hash same, counts=3
    kmer_count_table.count("CCCC")  # Add distinct k-mer, counts=1

    expected_histo = [(1, 1), (3, 1)]  # 1 k-mer observed once, 1 observed thrice
    assert (
        kmer_count_table.histo(zero=False) == expected_histo
    ), "histo(zero=False) should only return observed frequencies"


def test_histo_zero_true_non_empty_table(kmer_count_table):
    """
    Test histo(zero=True) on a non-empty KmerCountTable.
    All frequencies up to the maximum count should be included, including zero frequencies.
    """
    kmer_count_table.count("AAAA")  # Add k-mer, counts=1
    kmer_count_table.count("AAAA")  # Add k-mer, counts=2
    kmer_count_table.count("TTTT")  # Add another k-mer, canonical hash same, counts=3
    kmer_count_table.count("CCCC")  # Add distinct k-mer, counts=1

    expected_histo = [
        (0, 0),
        (1, 1),
        (2, 0),
        (3, 1),
    ]  # Include 0 frequency, 1 k-mer observed once, 0 observed twice, 1 observed thrice
    assert (
        kmer_count_table.histo(zero=True) == expected_histo
    ), "histo(zero=True) should include all frequencies up to max"


def test_histo_with_large_max_count(kmer_count_table):
    """Test histo() when there is a large maximum count in the table.

    Edge case: The histogram should correctly account for large frequency values.
    """
    for _ in range(5):
        kmer_count_table.count("AAAA")  # Add the same k-mer 100 times

    expected_histo = [
        (0, 0),
        (1, 0),
        (2, 0),
        (3, 0),
        (4, 0),
        (5, 1),
    ]  # 1 k-mer with count 100
    assert (
        kmer_count_table.histo(zero=True) == expected_histo
    ), "histo() include all zero counts up to max observed count."
