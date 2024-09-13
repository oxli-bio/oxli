import oxli
import pytest
from test_basic import create_sample_kmer_table

# Test attributes

def test_hashes_attribute():
    table = create_sample_kmer_table(3, ["AAA", "TTT", "AAC"])
    hashes = table.hashes
    hash_aaa = table.hash_kmer("AAA")  # 10679328328772601858
    hash_ttt = table.hash_kmer("TTT")  # 10679328328772601858
    hash_aac = table.hash_kmer("AAC")  # 6579496673972597301

    expected_hashes = set(
        [hash_aaa, hash_ttt, hash_aac]
    )  # {10679328328772601858, 6579496673972597301}
    assert (
        set(hashes) == expected_hashes
    ), ".hashes attribute should match the expected set of hash keys"


def test_version_attr():
    '''Check version attribute matches current version.'''
    pass

def test_total_consumed_seq_len_attr():
    '''Should log total seq len consumed.'''
    # Individual kmers
    # Long seqs with multiple kmers
    # Exclude invalid kmers?
    pass