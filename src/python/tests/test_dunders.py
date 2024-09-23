import pytest

import oxli
from test_basic import create_sample_kmer_table


# Test __len__
def test_len_initial():
    kmer_table = oxli.KmerCountTable(ksize=16)
    assert len(kmer_table) == 0, "Initial length should be 0"


def test_len_after_count():
    kmer_table = oxli.KmerCountTable(ksize=16)
    kmer_table.count("ACGTACGTACGTACGT")  # Adds 1 unique k-mer
    assert len(kmer_table) == 1, "Length should be 1 after adding one unique k-mer"


def test_len_after_multiple_counts():
    kmer_table = oxli.KmerCountTable(ksize=16)
    kmer_table.count("ACGTACGTACGTACGT")  # Adds 1 unique k-mer
    kmer_table.count("ACGTACGTACGTACGT")  # Adds 1 repeat k-mer
    kmer_table.count("CCCCCCCCCCCCCCCC")  # Adds 1 unique k-mer
    kmer_table.consume("GCTAGCTAGCTA")  # Adds 0 k-mers
    assert len(kmer_table) == 2, "Length should be 2 after adding two unique k-mers"


# Test iter methods


def test_iterable():
    kmer_table = create_sample_kmer_table(3, ["AAA", "TTT", "AAC"])
    hash_aaa = kmer_table.hash_kmer("AAA")  # 10679328328772601858
    hash_ttt = kmer_table.hash_kmer("TTT")  # 10679328328772601858
    hash_aac = kmer_table.hash_kmer("AAC")  # 6579496673972597301

    # Collect items from the iterator
    items = list(kmer_table)

    # Check if the items contain the expected tuples
    assert 2 in [
        count for _, count in items
    ], "Counts should be present in the iterated items"
    assert 6579496673972597301 in [
        key for key, _ in items
    ], "keys should be present in the iterated items"
    assert len(items) == 2, "There should be 2 k-mers in the table"


def test_iter_empty():
    kmer_table = oxli.KmerCountTable(ksize=16)

    # Collect items from an empty iterator
    items = list(kmer_table)

    # Ensure that no items are returned from an empty table
    assert items == [], "Iterator should be empty for an empty KmerCountTable"


# Test __set__ and __get__


def test_setitem():
    """Set values using the indexing syntax (obj[key] = value)"""
    kmer_table = oxli.KmerCountTable(ksize=16)
    kmer_table["ACGTACGTACGTACGT"] = 5  # Set count directly
    assert (
        kmer_table["ACGTACGTACGTACGT"] == 5
    ), "Value should be 5 after setting with __setitem__"


def test_getitem():
    """Query an object to using the indexing syntax (obj[key])"""
    kmer_table = oxli.KmerCountTable(ksize=16)
    kmer_table["ACGTACGTACGTACGT"] = 5
    assert (
        kmer_table["ACGTACGTACGTACGT"] == 5
    ), "Value should be 5 after setting with __setitem__"
    assert kmer_table["ACGTACGTACGTACGT"] == kmer_table.get(
        "ACGTACGTACGTACGT"
    ), "Behaviour should be same as .get()"

    # Check for a k-mer that does not exist
    assert (
        kmer_table["CCCCCCCCCCCCCCCC"] == 0
    ), "Default value for non-existent k-mer should be 0"


def test_setitem_update():
    kmer_table = oxli.KmerCountTable(ksize=16)
    kmer_table.count("ACGTACGTACGTACGT")  # Set count to 1
    kmer_table["ACGTACGTACGTACGT"] = 5  # Update count to 5
    assert kmer_table.get("ACGTACGTACGTACGT") == 5
    kmer_table["ACGTACGTACGTACGT"] = 10  # Update the count
    assert (
        kmer_table["ACGTACGTACGTACGT"] == 10
    ), "Value should be updated to 10 after setting with __setitem__"
