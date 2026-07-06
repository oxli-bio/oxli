import collections
import pytest

import oxli


# Helper function, create tables.
def create_sample_kmer_table(ksize, kmers):
    table = oxli.KmerCountTable(ksize)
    for kmer in kmers:
        table.count(kmer)
    return table


def test_consume_1_chunk():
    # test basic consume
    cg = oxli.KmerCountTable(4)
    kmer = "ATCG"

    cg.parallel_consume(kmer, 4)
    assert cg.get("ATCG") == 1


def test_consume_2():
    seq = "TAAACCCTAACCCTAACCCTAACCCTAACCC"
    parallel_chunk_size = len(seq) // 2

    cg1 = oxli.KmerCountTable(ksize=4)

    kh = cg1.kmers_and_hashes(seq, False)

    kmer_counter = collections.Counter()
    for kmer, hashval in kh:
        kmer_counter[kmer] += 1

    cg1.parallel_consume(seq, parallel_chunk_size)

    mismatch = False
    for kmer, count in kmer_counter.most_common():
        if cg1.get(kmer) != count:
            mismatch = True

    assert not mismatch


def test_parallel_consume_matches_consume():
    """parallel_consume produces identical counts to consume."""
    seq = "TAAACCCTAACCCTAACCCTAACCCTAACCC"

    cg1 = oxli.KmerCountTable(ksize=4)
    cg2 = oxli.KmerCountTable(ksize=4)

    n1 = cg1.consume(seq)
    n2 = cg2.parallel_consume(seq, 10)

    assert n1 == n2
    # All k-mer counts should match
    for hashval in cg1.hashes:
        assert cg1.get_hash(hashval) == cg2.get_hash(hashval)


def test_parallel_consume_consumed_attr():
    """consumed attribute tracks total bases processed."""
    seq = "TAAACCCTAACCCTAACCCTAACCCTAACCC"

    cg1 = oxli.KmerCountTable(ksize=4)
    cg2 = oxli.KmerCountTable(ksize=4)

    cg1.consume(seq)
    cg2.parallel_consume(seq, 10)

    assert cg1.consumed == cg2.consumed
    assert cg1.consumed == len(seq)


def test_parallel_consume2_matches_consume():
    """parallel_consume2 produces identical counts to consume."""
    seq = "TAAACCCTAACCCTAACCCTAACCCTAACCC"

    cg1 = oxli.KmerCountTable(ksize=4)
    cg2 = oxli.KmerCountTable(ksize=4)

    n1 = cg1.consume(seq)
    n2 = cg2.parallel_consume2(seq, 10)

    assert n1 == n2
    # All k-mer counts should match
    for hashval in cg1.hashes:
        assert cg1.get_hash(hashval) == cg2.get_hash(hashval)


def test_parallel_consume2_consumed_attr():
    """parallel_consume2 consumed attribute tracks total bases."""
    seq = "TAAACCCTAACCCTAACCCTAACCCTAACCC"

    cg1 = oxli.KmerCountTable(ksize=4)
    cg2 = oxli.KmerCountTable(ksize=4)

    cg1.consume(seq)
    cg2.parallel_consume2(seq, 10)

    assert cg1.consumed == cg2.consumed
    assert cg2.consumed == len(seq)


def test_parallel_consume_short_seq():
    """Sequence shorter than chunk_size works correctly."""
    seq = "ATCG"
    cg = oxli.KmerCountTable(4)
    n = cg.parallel_consume(seq, 50000)
    assert n == 1
    assert cg.get("ATCG") == 1


def test_parallel_consume_with_store_kmers():
    """parallel_consume preserves hash-to-kmer mappings."""
    seq = "TAAACCCTAACCCTAACCCTAACCCTAACCC"

    cg1 = oxli.KmerCountTable(ksize=4, store_kmers=True)
    cg2 = oxli.KmerCountTable(ksize=4, store_kmers=True)

    cg1.consume(seq)
    cg2.parallel_consume(seq, 10)

    # Both should have the same k-mer set
    for hashval in cg1.hashes:
        kmer1 = cg1.unhash(hashval)
        kmer2 = cg2.unhash(hashval)
        assert kmer1 == kmer2
