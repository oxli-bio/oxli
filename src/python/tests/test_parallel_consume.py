import collections
import pytest

import oxli


# Helper function, create tables.
def create_sample_kmer_table(ksize, kmers):
    table = oxli.KmerCountTable(ksize)
    for kmer in kmers:
        table.count(kmer)
    return table

# @CTB check 'consumed'

def test_consume_1_chunk():
    # test basic consume
    cg = oxli.KmerCountTable(4)
    kmer = "ATCG"

    cg.parallel_consume(kmer, 4)
    assert cg.get("ATCG") == 1


def test_consume_2():
    seq = "TAAACCCTAACCCTAACCCTAACCCTAACCC"
    #parallel_chunk_size = 4
    parallel_chunk_size = len(seq) // 2

    cg1 = oxli.KmerCountTable(ksize=4)
    cg2 = oxli.KmerCountTable(ksize=4)

    kh = cg1.kmers_and_hashes(seq, False)

    kmer_counter = collections.Counter()
    hash_counter = collections.Counter()
    for kmer, hashval in kh:
        kmer_counter[kmer] += 1
        hash_counter[kmer] += 1

    cg1.parallel_consume(seq, parallel_chunk_size)

    mismatch = False
    for kmer, count in kmer_counter.most_common():
        print(cg1.get(kmer) == count, kmer, count, cg1.get(kmer))
        if cg1.get(kmer) != count:
            mismatch = True

    assert not mismatch
