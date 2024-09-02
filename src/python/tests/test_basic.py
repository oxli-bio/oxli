import pytest
import oxli

def test_simple():
    # yo dawg it works
    cg = oxli.KmerCountTable(4)
    kmer = "ATCG"

    assert cg.get(kmer) == 0
    assert cg.count(kmer) == 1
    assert cg.get(kmer) == 1


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
    assert cg.get("CCGA") == 1 # reverse complement!


def test_consume_bad_DNA():
    # test an invalid base in last position
    cg = oxli.KmerCountTable(4)
    seq = "ATCGGX"
    with pytest.raises(ValueError,
                       match="bad k-mer encountered at position 2"):
        cg.consume(seq, allow_bad_kmers=False)


def test_consume_bad_DNA_2():
    # test an invalid base in first position
    cg = oxli.KmerCountTable(4)
    seq = "XATCGG"
    with pytest.raises(ValueError,
                       match="bad k-mer encountered at position 0"):
        cg.consume(seq, allow_bad_kmers=False)


def test_consume_bad_DNA_ignore():
    # we can ignore bad DNA
    cg = oxli.KmerCountTable(4)
    seq = "XATCGG"
    print(cg.consume(seq, allow_bad_kmers=True))
    assert cg.get("ATCG") == 1
    assert cg.get("TCGG") == 1
    assert cg.get("CCGA") == 1 # rc


def test_consume_bad_DNA_ignore_is_default():
    # ignoring bad DNA is default
    cg = oxli.KmerCountTable(4)
    seq = "XATCGG"
    print(cg.consume(seq))
    assert cg.get("ATCG") == 1
    assert cg.get("TCGG") == 1
    assert cg.get("CCGA") == 1 # rc
