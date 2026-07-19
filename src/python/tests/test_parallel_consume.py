import collections
import pytest
import random
import string

import oxli


# Helper function, create tables.
def create_sample_kmer_table(ksize, kmers):
    table = oxli.KmerCountTable(ksize)
    for kmer in kmers:
        table.count(kmer)
    return table


def random_dna(length, seed=42):
    """Generate a random DNA sequence of the given length."""
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


# ── basic correctness ─────────────────────────────────────────────────────────


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


# ── chunk-size edge cases ─────────────────────────────────────────────────────


def test_parallel_consume_chunk_equals_seq_len():
    """chunk_size == len(seq) processes the whole sequence in one chunk."""
    seq = "TAAACCCTAACCCTAACCCTAACCCTAACCC"

    cg1 = oxli.KmerCountTable(ksize=4)
    cg2 = oxli.KmerCountTable(ksize=4)

    cg1.consume(seq)
    cg2.parallel_consume(seq, len(seq))

    assert cg1.consumed == cg2.consumed
    for hashval in cg1.hashes:
        assert cg1.get_hash(hashval) == cg2.get_hash(hashval)


def test_parallel_consume_chunk_larger_than_seq():
    """chunk_size > len(seq) behaves like a single consume call."""
    seq = "TAAACCCTAACCCTAACCCTAACCCTAACCC"

    cg1 = oxli.KmerCountTable(ksize=4)
    cg2 = oxli.KmerCountTable(ksize=4)

    n1 = cg1.consume(seq)
    n2 = cg2.parallel_consume(seq, len(seq) * 10)

    assert n1 == n2
    for hashval in cg1.hashes:
        assert cg1.get_hash(hashval) == cg2.get_hash(hashval)


def test_parallel_consume_chunk_equals_ksize():
    """chunk_size equal to ksize (minimum valid chunk) works correctly."""
    seq = "TAAACCCTAACCCTAACCCTAACCCTAACCC"
    ksize = 4

    cg1 = oxli.KmerCountTable(ksize=ksize)
    cg2 = oxli.KmerCountTable(ksize=ksize)

    cg1.consume(seq)
    cg2.parallel_consume(seq, ksize)

    assert cg1.consumed == cg2.consumed
    for hashval in cg1.hashes:
        assert cg1.get_hash(hashval) == cg2.get_hash(hashval)


def test_parallel_consume_many_small_chunks():
    """Many small chunks still produce the correct result."""
    seq = "TAAACCCTAACCCTAACCCTAACCCTAACCC"
    ksize = 4

    cg1 = oxli.KmerCountTable(ksize=ksize)
    cg2 = oxli.KmerCountTable(ksize=ksize)

    n1 = cg1.consume(seq)
    # chunk_size=5 → many chunks for a 31-base sequence
    n2 = cg2.parallel_consume(seq, 5)

    assert n1 == n2
    for hashval in cg1.hashes:
        assert cg1.get_hash(hashval) == cg2.get_hash(hashval)


# ── k-mer size variations ─────────────────────────────────────────────────────


@pytest.mark.parametrize("ksize", [3, 7, 15, 21, 31])
def test_parallel_consume_various_ksizes(ksize):
    """parallel_consume matches consume for various k-mer sizes."""
    seq = random_dna(500, seed=ksize)

    cg1 = oxli.KmerCountTable(ksize=ksize)
    cg2 = oxli.KmerCountTable(ksize=ksize)

    n1 = cg1.consume(seq)
    n2 = cg2.parallel_consume(seq, 100)

    assert n1 == n2
    for hashval in cg1.hashes:
        assert cg1.get_hash(hashval) == cg2.get_hash(hashval), (
            f"Count mismatch for ksize={ksize}"
        )


# ── longer / random sequences ─────────────────────────────────────────────────


@pytest.mark.parametrize(
    "seq_len,chunk_size",
    [
        (1000, 100),
        (1000, 317),  # chunk_size not a divisor of seq_len
        (5000, 500),
        (10_000, 1000),
    ],
)
def test_parallel_consume_random_seq(seq_len, chunk_size):
    """parallel_consume matches consume for random sequences of various sizes."""
    seq = random_dna(seq_len)

    cg1 = oxli.KmerCountTable(ksize=21)
    cg2 = oxli.KmerCountTable(ksize=21)

    n1 = cg1.consume(seq)
    n2 = cg2.parallel_consume(seq, chunk_size)

    assert n1 == n2
    assert cg1.consumed == cg2.consumed == len(seq)
    for hashval in cg1.hashes:
        assert cg1.get_hash(hashval) == cg2.get_hash(hashval)


# ── bad k-mer handling ────────────────────────────────────────────────────────


def test_parallel_consume_skip_bad_kmers():
    """Bad k-mers are skipped when skip_bad_kmers=True (default)."""
    seq = "ATCGNATCGATCG"  # N in middle

    cg1 = oxli.KmerCountTable(ksize=4)
    cg2 = oxli.KmerCountTable(ksize=4)

    n1 = cg1.consume(seq, skip_bad_kmers=True)
    n2 = cg2.parallel_consume(seq, 5, skip_bad_kmers=True)

    # Both should have the same k-mers (bad ones skipped)
    assert cg1.consumed == cg2.consumed
    for hashval in cg1.hashes:
        assert cg1.get_hash(hashval) == cg2.get_hash(hashval)


# ── sum_counts sanity ─────────────────────────────────────────────────────────


def test_parallel_consume_sum_counts():
    """sum_counts equals the number returned by parallel_consume for clean DNA."""
    seq = random_dna(200)
    ksize = 7

    cg = oxli.KmerCountTable(ksize=ksize)
    n = cg.parallel_consume(seq, 50)

    # sum_counts may exceed n when k-mers repeat; n is the number of k-mer
    # positions processed, sum_counts is the total of all counts.
    assert cg.sum_counts == n


# ── default chunk_size ────────────────────────────────────────────────────────


def test_parallel_consume_default_chunk_size():
    """parallel_consume uses a sensible default chunk_size."""
    seq = random_dna(500)

    cg1 = oxli.KmerCountTable(ksize=21)
    cg2 = oxli.KmerCountTable(ksize=21)

    cg1.consume(seq)
    cg2.parallel_consume(seq)  # default chunk_size=50_000

    for hashval in cg1.hashes:
        assert cg1.get_hash(hashval) == cg2.get_hash(hashval)
