import pytest

from oxli import KmerCountTable


def test_add_basic():
    """Test basic addition of two KmerCountTables."""
    table1 = KmerCountTable(5)
    table2 = KmerCountTable(5)

    table1.consume("ATGCATGCA")
    table2.consume("ATGCATGCA")

    counts_added, new_keys = table1.add(table2)

    assert counts_added == 5  # 5 kmers in "ATGCATGCA"
    assert new_keys == 0  # No new keys added
    assert table1.sum_counts == 10  # Each kmer counted twice


def test_add_different_content():
    """Test addition of KmerCountTables with different content."""
    table1 = KmerCountTable(5)
    table2 = KmerCountTable(5)

    table1.consume("ATGCATGCA")
    table2.consume("TGCATGCATGG")

    counts_added, new_keys = table1.add(table2)

    assert len(table1) == 3
    assert counts_added == 7
    assert new_keys == 1  # One new kmer: "CATGG" "3442404512935954368"
    assert table1.sum_counts == 12


def test_add_different_ksize():
    """Test addition of KmerCountTables with different ksizes."""
    table1 = KmerCountTable(5)
    table2 = KmerCountTable(6)

    with pytest.raises(ValueError):
        table1.add(table2)


def test_add_empty_tables():
    """Test addition of empty KmerCountTables."""
    table1 = KmerCountTable(5)
    table2 = KmerCountTable(5)

    counts_added, new_keys = table1.add(table2)

    assert counts_added == 0
    assert new_keys == 0
    assert table1.sum_counts == 0


def test_add_to_empty_table():
    """Test addition to an empty KmerCountTable."""
    table1 = KmerCountTable(5)
    table2 = KmerCountTable(5)

    table2.consume("ATGCATGCA")

    counts_added, new_keys = table1.add(table2)

    assert len(table1) == 2
    assert counts_added == 5
    assert new_keys == 2  # All keys are new
    assert table1.sum_counts == 5


def test_add_consumed_attribute():
    """Test that the 'consumed' attribute is correctly updated."""
    table1 = KmerCountTable(5)
    table2 = KmerCountTable(5)

    table1.consume("ATGCA")
    table2.consume("TGCAT")

    initial_consumed = table1.consumed
    table1.add(table2)

    assert table1.consumed == initial_consumed + table2.consumed


@pytest.mark.parametrize(
    "store_kmers1,store_kmers2",
    [(True, True), (True, False), (False, True), (False, False)],
)
def test_add_store_kmers_combinations(store_kmers1, store_kmers2, capfd):
    """Test addition with different combinations of store_kmers option."""
    table1 = KmerCountTable(5, store_kmers=store_kmers1)
    table2 = KmerCountTable(5, store_kmers=store_kmers2)

    table1.consume("ATGCA")
    table2.consume("GGCAT")

    counts_added, new_keys = table1.add(table2)

    assert counts_added == 1
    assert new_keys == 1

    captured = capfd.readouterr()
    if store_kmers1 and not store_kmers2:
        assert "Warning: Incoming table does not store k-mers" in captured.err

    if store_kmers1 and store_kmers2:
        assert table1.dump_kmers(sortkeys=True) == [("ATGCA", 1), ("ATGCC", 1)]


def test_add_large_tables():
    """Test addition of large KmerCountTables to check performance."""
    table1 = KmerCountTable(5)
    table2 = KmerCountTable(5)

    long_seq = "ATGC" * 100000  # 400,000 base pairs
    table1.consume(long_seq)
    table2.consume(long_seq)

    counts_added, new_keys = table1.add(table2)

    assert counts_added == 399996  # (400000 - 5 + 1) kmers
    assert new_keys == 0
    assert table1.sum_counts == 799992  # Each kmer counted twice


def test_add_multiple_times():
    """Test adding multiple KmerCountTables."""
    table1 = KmerCountTable(5)
    table2 = KmerCountTable(5)
    table3 = KmerCountTable(5)

    table1.consume("ATGCA")
    table2.consume("TGCAT")
    table3.consume("GCATG")

    table1.add(table2)
    counts_added, new_keys = table1.add(table3)

    assert counts_added == 1
    assert new_keys == 1
    assert table1.sum_counts == 3


# def test_add_self():
#    """Test adding a KmerCountTable to itself."""
#    table = KmerCountTable(5)
#    table.consume("ATGCATGCA")
#
#    with pytest.raises(ValueError, match="Cannot add KmerCountTable to itself."):
#        table.add(table)
#
#    # Ensure the original table is unchanged
#    assert table.sum_counts() == 5  # 5 kmers in "ATGCATGCA"


# Run the tests
if __name__ == "__main__":
    pytest.main([__file__])
