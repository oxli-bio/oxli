import pytest

import oxli


# Helper function, create tables.
def create_sample_kmer_table(ksize, kmers):
    table = oxli.KmerCountTable(ksize)
    for kmer in kmers:
        table.count(kmer)
    return table


# Adding Kmers
def test_count():
    # yo dawg it works
    cg = oxli.KmerCountTable(4)
    kmer = "ATCG"

    assert cg.get(kmer) == 0
    assert cg.count(kmer) == 1
    assert cg.get(kmer) == 1


def test_count_hash():
    kmer = "TAAACCCTAACCCTAACCCTAACCCTAACCC"
    cg = oxli.KmerCountTable(ksize=31)
    hashkey = cg.hash_kmer(kmer)

    assert cg.get_hash(hashkey) == 0
    assert cg.count_hash(hashkey) == 1
    assert cg.get_hash(hashkey) == 1


def test_hash_rc():
    table = create_sample_kmer_table(3, ["AAA", "TTT", "AAC"])
    hash_aaa = table.hash_kmer("AAA")  # 10679328328772601858
    hash_ttt = table.hash_kmer("TTT")  # 10679328328772601858

    assert hash_aaa == hash_ttt, "Hash should be same for reverse complement."


def test_wrong_ksize():
    # but only with the right ksize
    cg = oxli.KmerCountTable(3)
    kmer = "ATCG"

    with pytest.raises(ValueError):
        cg.count(kmer)

    with pytest.raises(ValueError):
        cg.get(kmer)


def test_consume():
    # test basic consume
    cg = oxli.KmerCountTable(4)
    kmer = "ATCG"

    assert cg.consume(kmer) == 1
    assert cg.get("ATCG") == 1


def test_consume_2():
    # test reverse complement
    cg = oxli.KmerCountTable(4)
    seq = "ATCGG"

    assert cg.consume(seq) == 2
    assert cg.get("ATCG") == 1
    assert cg.get("TCGG") == 1
    assert cg.get("CCGA") == 1  # reverse complement!


def test_consume_bad_DNA():
    # test an invalid base in last position
    cg = oxli.KmerCountTable(4)
    seq = "ATCGGX"
    with pytest.raises(ValueError, match="bad k-mer encountered at position 2"):
        cg.consume(seq, skip_bad_kmers=False)


def test_consume_bad_DNA_2():
    # test an invalid base in first position
    cg = oxli.KmerCountTable(4)
    seq = "XATCGG"
    with pytest.raises(ValueError, match="bad k-mer encountered at position 0"):
        cg.consume(seq, skip_bad_kmers=False)


def test_consume_bad_DNA_ignore():
    # we can ignore bad DNA
    cg = oxli.KmerCountTable(4)
    seq = "XATCGG"
    print(cg.consume(seq, skip_bad_kmers=True))
    assert cg.get("ATCG") == 1
    assert cg.get("TCGG") == 1
    assert cg.get("CCGA") == 1  # rc


def test_consume_bad_DNA_ignore_is_default():
    # ignoring bad DNA is default
    cg = oxli.KmerCountTable(4)
    seq = "XATCGG"
    print(cg.consume(seq))
    assert cg.get("ATCG") == 1
    assert cg.get("TCGG") == 1
    assert cg.get("CCGA") == 1  # rc


# Getting counts
def test_count_vs_counthash():
    # test a bug reported by adam taranto: count and get should work together!
    kmer = "TAAACCCTAACCCTAACCCTAACCCTAACCC"
    cg = oxli.KmerCountTable(ksize=31)
    hashkey = cg.hash_kmer(kmer)

    assert cg.get(kmer) == 0
    assert cg.count(kmer) == 1
    assert cg.count(kmer) == 2
    assert cg.get(kmer) == 2
    assert cg.count_hash(hashkey) == 3

    x = cg.get(kmer)
    assert x == 3, x


def test_get_hash():
    """Retrieve counts using hash key."""
    table = create_sample_kmer_table(3, ["AAA", "TTT", "AAC"])
    # Find hash of kmer 'AAA'
    hash_aaa = table.hash_kmer("AAA")  # 10679328328772601858
    # Lookup counts for hash of 'AAA' and rc 'TTT'
    count_aaa = table.get_hash(hash_aaa)
    assert count_aaa == 2, "Hash count for 'AAA' should be 2"

    # Test single kmer
    hash_aac = table.hash_kmer("AAC")  # 6579496673972597301
    count_aac = table.get_hash(hash_aac)
    assert count_aac == 1, "Hash count for 'AAC' should be 1"

    # Test for kmer that is not in table
    hash_aag = table.hash_kmer("AAG")  # 12774992397053849803
    count_aag = table.get_hash(hash_aag)
    assert count_aag == 0, "Missing kmer count for 'AAG' should be 0"


def test_get_hash_array():
    """
    Get vector of counts corresponding to vector of hash keys.
    """
    table = create_sample_kmer_table(3, ["AAA", "TTT", "AAC"])
    hash_aaa = table.hash_kmer("AAA")
    hash_aac = table.hash_kmer("AAC")
    hash_ggg = table.hash_kmer("GGG")  # key not in table

    hash_keys = [hash_aaa, hash_aac, hash_ggg]
    hash_keys_rev = [hash_ggg, hash_aac, hash_aaa]

    counts = table.get_hash_array(hash_keys)
    rev_counts = table.get_hash_array(hash_keys_rev)

    assert counts == [2, 1, 0], (
        "Hash array counts should match the counts of 'AAA' and 'AAC' and return zero for 'GGG'."
    )
    assert rev_counts == [0, 1, 2], "Count should be in same order as input list"


# def test_get_array():
#    """
#    Get vector of counts corresponding to vector of kmers.
#    """
#    pass
