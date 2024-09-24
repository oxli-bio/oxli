import pytest

import oxli


# Helper function, create tables.
def create_sample_kmer_table(ksize, kmers):
    table = oxli.KmerCountTable(ksize)
    for kmer in kmers:
        table.count(kmer)
    return table


def test_basic():
    "string containing only forward canonical kmers."
    seq = "ATAAACC"             # all forward k-mers
    cg = oxli.KmerCountTable(ksize=4)

    x = cg.kmers_and_hashes(seq, False)
    assert x == [('ATAA', 179996601836427478),
                 ('TAAA', 15286642655859448092),
                 ('AAAC', 9097280691811734508),
                 ('AACC', 6779379503393060785)]


def test_basic_rc():
    "string containing only reverse canonical kmers."
    seq = "GGTTTAT"
    cg = oxli.KmerCountTable(ksize=4)

    x = cg.kmers_and_hashes(seq, False)
    print(x)
    assert x == [('AACC', 6779379503393060785),
                 ('AAAC', 9097280691811734508),
                 ('TAAA', 15286642655859448092),
                 ('ATAA', 179996601836427478)]


def test_basic_mixed():
    "string containing forward and reverse canonical kmers."
    seq = "ACGTTG"
    cg = oxli.KmerCountTable(ksize=4)

    x = cg.kmers_and_hashes(seq, False)
    print(x)
    assert x == [('ACGT', 2597925387403686983),
                 ('AACG', 7952982457453691616),
                 ('CAAC', 7315150081962684964)]


def test_basic_lower():
    "Test that sequences are turned into uppercase appropriately."
    seq = "acgttg"
    cg = oxli.KmerCountTable(ksize=4)

    x = cg.kmers_and_hashes(seq, False)
    print(x)
    assert x == [('ACGT', 2597925387403686983),
                 ('AACG', 7952982457453691616),
                 ('CAAC', 7315150081962684964)]


def test_bad_kmers_raise_error():
    "Test that sequences are turned into uppercase appropriately."
    seq = "acxttg"
    cg = oxli.KmerCountTable(ksize=4)

    # CTB: would be nice to turn this into a better error message,
    # ideally one containing the bad k-mer and/or location.
    try:
        x = cg.kmers_and_hashes(seq, False)
        assert False, "this should fail"
    except:
        pass


def test_bad_kmers_allowed():
    "Test that sequences are turned into uppercase appropriately."
    seq = "aattxttgg"
    cg = oxli.KmerCountTable(ksize=4)

    x = cg.kmers_and_hashes(seq, True)
    print(x)
    assert x == [('AATT', 382727017318141683),
                 ('', 0), ('', 0), ('', 0), ('', 0),
                 ('CCAA', 1798905482136869687)]
