from os import remove
import tempfile

import pytest

from oxli import KmerCountTable


@pytest.fixture
def kmer_count_table():
    """Fixture to set up a KmerCountTable instance with sample data."""
    kct = KmerCountTable(ksize=4, store_kmers=True)
    kct.count("AAAA")  # 17832910516274425539
    kct.count("TTTT")  # 17832910516274425539
    kct.count("AATT")  # 382727017318141683
    kct.count("GGGG")  # 73459868045630124
    kct.count("GGGG")  # 73459868045630124
    return kct


@pytest.fixture
def empty_kmer_count_table():
    """Fixture to set up an empty KmerCountTable instance."""
    return KmerCountTable(ksize=4, store_kmers=True)


def test_dump_conflicting_sort_options(kmer_count_table):
    """Test that passing both sortcounts=True and sortkeys=True raises a ValueError."""
    with pytest.raises(
        ValueError, match="Cannot sort by both counts and keys at the same time."
    ):
        kmer_count_table.dump(file=None, sortcounts=True, sortkeys=True)


def test_dump_no_sorting(kmer_count_table):
    """Test the dump function with no sorting (both sortcounts and sortkeys are False)."""
    result = kmer_count_table.dump(file=None, sortcounts=False, sortkeys=False)

    # Expected output same order as for iterator
    expected = list(kmer_count_table)
    # [(17832910516274425539, 2), (382727017318141683, 1), (73459868045630124, 2)]

    assert result == expected, f"Expected {expected}, but got {result}"


def test_dump_sortcounts_with_ties(kmer_count_table):
    """Test the dump function with sortcounts=True, ensuring it handles ties in counts."""
    result = kmer_count_table.dump(file=None, sortcounts=True, sortkeys=False)

    # Expected output sorted by count, with secondary sorting by hash for ties
    expected = [
        (382727017318141683, 1),  # 'AATT'
        (73459868045630124, 2),  # 'GGGG' (lower hash than 'AAAA')
        (17832910516274425539, 2),  # 'AAAA'/'TTTT'
    ]

    assert result == expected, f"Expected {expected}, but got {result}"


def test_dump_single_kmer():
    """Test the dump function with only a single k-mer counted."""
    kct = KmerCountTable(ksize=4)
    kct.count("AAAA")  # Hash for 'AAAA'/'TTTT'

    result = kct.dump(file=None, sortcounts=True, sortkeys=False)

    expected = [
        (17832910516274425539, 1)  # 'AAAA'/'TTTT'
    ]

    assert result == expected, f"Expected {expected}, but got {result}"


def test_dump_write_to_file(kmer_count_table):
    """Test the dump function when writing to a file.

    This test checks if the function correctly writes the hash:count pairs to a file.
    """
    with tempfile.NamedTemporaryFile(delete=False) as temp_file:
        temp_file_path = temp_file.name

    kmer_count_table.dump(file=temp_file_path, sortcounts=True, sortkeys=False)

    with open(temp_file_path, "r") as f:
        lines = f.readlines()

    # Expected output sorted by count then hash (default behavior)
    expected_lines = [
        f"{382727017318141683}\t1\n",  # 'AATT'
        f"{73459868045630124}\t2\n",  # 'GGGG'
        f"{17832910516274425539}\t2\n",  # 'AAAA'/'TTTT'
    ]

    assert lines == expected_lines, f"Expected {expected_lines}, but got {lines}"

    # Cleanup
    remove(temp_file_path)


def test_dump_write_to_file_sortkeys(kmer_count_table):
    """Test the dump function with sortkeys=True when writing to a file."""
    with tempfile.NamedTemporaryFile(delete=False) as temp_file:
        temp_file_path = temp_file.name

    kmer_count_table.dump(file=temp_file_path, sortkeys=True)

    with open(temp_file_path, "r") as f:
        lines = f.readlines()

    # Expected output sorted by hash keys
    expected_lines = [
        f"{73459868045630124}\t2\n",  # 'GGGG'
        f"{382727017318141683}\t1\n",  # 'AATT'
        f"{17832910516274425539}\t2\n",  # 'AAAA'/'TTTT'
    ]

    assert lines == expected_lines, f"Expected {expected_lines}, but got {lines}"

    # Cleanup
    remove(temp_file_path)


def test_dump_sortkeys(kmer_count_table):
    """Test the dump function with sortkeys=True.

    This test verifies if the function sorts by hash keys when `sortkeys` is set to True.
    """
    result = kmer_count_table.dump(file=None, sortkeys=True)

    # Expected output sorted by hash key
    expected = [
        (73459868045630124, 2),  # 'GGGG'
        (382727017318141683, 1),  # 'AATT'
        (17832910516274425539, 2),  # 'AAAA'/'TTTT'
    ]

    assert result == expected, f"Expected {expected}, but got {result}"


def test_dump_invalid_file_path(kmer_count_table):
    """Test that passing an invalid file path raises an error."""
    with pytest.raises(OSError):
        kmer_count_table.dump(file="", sortkeys=True)


def test_dump_hash_empty_table(empty_kmer_count_table):
    """Test the dump function on an empty KmerCountTable.

    This test checks that the function handles an empty table correctly.
    """
    # Test that calling dump without file returns an empty list
    result = empty_kmer_count_table.dump(file=None, sortkeys=False)
    assert result == [], "Expected an empty list from an empty KmerCountTable"

    # Test that calling dump with a file writes nothing to the file
    with tempfile.NamedTemporaryFile(delete=False) as temp_file:
        temp_file_path = temp_file.name

    empty_kmer_count_table.dump(file=temp_file_path, sortkeys=False)

    with open(temp_file_path, "r") as f:
        lines = f.readlines()

    assert lines == [], "Expected an empty file for an empty KmerCountTable"

    # Cleanup
    remove(temp_file_path)


# Tests for dump_kmers()


def test_dump_kmers_conflicting_sort_options(kmer_count_table):
    """Test that passing both sortcounts=True and sortkeys=True raises a ValueError."""
    with pytest.raises(
        ValueError, match="Cannot sort by both counts and kmers at the same time."
    ):
        kmer_count_table.dump_kmers(file=None, sortcounts=True, sortkeys=True)


def test_dump_kmers_sortcounts_with_ties(kmer_count_table):
    """Test the dump_kmers function with sortcounts=True, ensuring it handles ties in counts."""
    result = kmer_count_table.dump_kmers(file=None, sortcounts=True, sortkeys=False)

    # Expected output sorted by count, with secondary sorting by kmer for ties
    expected = [
        ("AATT", 1),
        (
            "AAAA",
            2,
        ),  # 'AAAA'/'TTTT' is tied with 'GGGG / CCCC' on counts, 'AAAA' is lexicographically smaller
        ("CCCC", 2),
    ]

    assert result == expected, f"Expected {expected}, but got {result}"


def test_dump_kmers_single_kmer():
    """Test the dump_kmers function with only a single k-mer counted."""
    kct = KmerCountTable(ksize=4, store_kmers=True)
    kct.count("AAAA")  # Canonical kmer: 'AAAA'

    result = kct.dump_kmers(file=None, sortcounts=True, sortkeys=False)

    expected = [("AAAA", 1)]

    assert result == expected, f"Expected {expected}, but got {result}"


def test_dump_kmers_write_to_file(kmer_count_table):
    """Test the dump_kmers function when writing to a file.

    This test checks if the function correctly writes the kmer:count pairs to a file.
    """
    with tempfile.NamedTemporaryFile(delete=False) as temp_file:
        temp_file_path = temp_file.name

    kmer_count_table.dump_kmers(file=temp_file_path, sortcounts=True, sortkeys=False)

    with open(temp_file_path, "r") as f:
        lines = f.readlines()

    # Expected output sorted by count then kmer (default behavior)
    expected_lines = [
        f"AATT\t1\n",
        f"AAAA\t2\n",  # 'AAAA'/'TTTT'
        f"CCCC\t2\n",
    ]

    assert lines == expected_lines, f"Expected {expected_lines}, but got {lines}"

    # Cleanup
    remove(temp_file_path)


def test_dump_kmers_write_to_file_sortkeys(kmer_count_table):
    """Test the dump_kmers function with sortkeys=True when writing to a file."""
    with tempfile.NamedTemporaryFile(delete=False) as temp_file:
        temp_file_path = temp_file.name

    kmer_count_table.dump_kmers(file=temp_file_path, sortkeys=True)

    with open(temp_file_path, "r") as f:
        lines = f.readlines()

    # Expected output sorted by canonical kmers
    expected_lines = [
        f"AAAA\t2\n",  # 'AAAA'/'TTTT'
        f"AATT\t1\n",
        f"CCCC\t2\n",
    ]

    assert lines == expected_lines, f"Expected {expected_lines}, but got {lines}"

    # Cleanup
    remove(temp_file_path)


def test_dump_kmers_sortkeys(kmer_count_table):
    """Test the dump_kmers function with sortkeys=True.

    This test verifies if the function sorts by canonical k-mers when `sortkeys` is set to True.
    """
    result = kmer_count_table.dump_kmers(file=None, sortkeys=True)

    # Expected output sorted by canonical kmer
    expected = [
        ("AAAA", 2),  # 'AAAA'/'TTTT'
        ("AATT", 1),
        ("CCCC", 2),
    ]

    assert result == expected, f"Expected {expected}, but got {result}"


def test_dump_kmers_invalid_file_path(kmer_count_table):
    """Test that passing an invalid file path raises an error."""
    with pytest.raises(OSError):
        kmer_count_table.dump_kmers(file="", sortkeys=True)


def test_dump_kmers_empty_table(empty_kmer_count_table):
    """Test the dump_kmers function on an empty KmerCountTable.

    This test checks that the function handles an empty table correctly.
    """
    # Test that calling dump_kmers without file returns an empty list
    result = empty_kmer_count_table.dump_kmers(file=None, sortkeys=False)
    assert result == [], "Expected an empty list from an empty KmerCountTable"

    # Test that calling dump_kmers with a file writes nothing to the file
    with tempfile.NamedTemporaryFile(delete=False) as temp_file:
        temp_file_path = temp_file.name

    empty_kmer_count_table.dump_kmers(file=temp_file_path, sortkeys=False)

    with open(temp_file_path, "r") as f:
        lines = f.readlines()

    assert lines == [], "Expected an empty file for an empty KmerCountTable"

    # Cleanup
    remove(temp_file_path)
