import oxli
import pytest
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


def test_iter_dunder_method():
    """KmerCountTable should be iterable, yield hash:count pairs"""
    pass


def test_next_dunder_method():
    """Select next key in generator"""
    pass


def test_getitem_dunder_method():
    """Query an object to using the indexing syntax (obj[key])"""
    # Same behaviour as .get()
    pass


def test_setitem_dunder_method():
    """Set values using the indexing syntax (obj[key] = value)"""
    pass
