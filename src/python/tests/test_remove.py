import pytest
import oxli

@pytest.fixture
def setup_kmer_table():
    """Fixture to set up a KmerCountTable with ksize=4 and some initial k-mers"""
    kct = oxli.KmerCountTable(ksize=4)
    kct.count("AAAA")  # Hash of canonical form will be used (AAAA)
    kct.count("CCCC")  # CCCC
    kct.count("ATAT")  # ATAT
    kct.count("GGGG")  # Should map to CCCC
    kct.count("TTTT")  # Should map to AAAA
    kct.count("CCCC")  # Increment count for CCCC/GGGG
    # AAAA/TTTT = 2
    # ATAT = 1
    # CCCC/GGGG = 3
    return kct

def test_drop(setup_kmer_table):
    """
    Test the drop method to remove a k-mer by its string representation.
    Edge case: Dropping a k-mer that doesn't exist.
    """
    kct = setup_kmer_table

    # Drop "GGGG" which exists, and check it's removed
    kct.drop("GGGG")
    assert kct.get("GGGG") == 0, "Expected 'GGGG' to be removed."

    # Drop "AAAA", should remove both "AAAA" and "TTTT" (same canonical form)
    kct.drop("AAAA")
    assert kct.get("AAAA") == 0, "Expected 'AAAA' (and 'TTTT') to be removed."

    # Edge case: Drop a k-mer that doesn't exist, e.g., "GGGA"
    kct.drop("GGGA")  # "GGGA" not present in the table
    assert kct.get("GGGA") == 0 # "GGGA" not present in the table
    
    # Raise error if kmer longer than table k len.
    with pytest.raises(ValueError):
        kct.drop("GGGAA")  # "GGGAA" longer than table k

def test_drop_hash(setup_kmer_table):
    """
    Test the drop_hash method to remove a k-mer by its hash.
    Edge case: Dropping a hash that doesn't exist.
    """
    kct = setup_kmer_table

    # Drop by the hash for "CCCC", and check it's removed
    hashval = kct.hash_kmer("CCCC")
    kct.drop_hash(hashval)
    assert kct.get_hash(hashval) == 0, "Expected 'CCCC' and 'GGGG' to be removed."
    assert kct.get("CCCC") == 0, "Expected 'CCCC' to be removed."
    assert kct.get("GGGG") == 0, "Expected 'GGGG' to be removed."

    # Edge case: Drop a hash that doesn't exist
    non_existent_hash = 999999999
    kct.drop_hash(non_existent_hash)  # Should not raise an error
    assert kct.get_hash(non_existent_hash) == 0, "Expected non-existent hash removal to succeed."

def test_mincut(setup_kmer_table):
    """
    Test the mincut method to remove all k-mers with counts less than a given threshold.
    Edge cases: Threshold is higher than all counts, no k-mers to remove.
    """
    kct = setup_kmer_table

    # Set a threshold that only removes k-mers with counts < 2
    removed = kct.mincut(3)
    assert removed == 2, "Expected 2 k-mers to be removed ('ATAT' and 'AAAA/TTTT')."
    assert kct.get("GGGG") == 3, "Expected 'GGGG/CCCC' to remain."
    
    # Edge case: Threshold is higher than all k-mer counts (remove everything)
    removed = kct.mincut(10)
    assert removed == 1, "Expected all remaining k-mers to be removed ('GGGG/CCCC')."
    assert len(kct.hashes) == 0, "Expected no k-mers left after removing all."

def test_maxcut(setup_kmer_table):
    """
    Test the maxcut method to remove all k-mers with counts greater than a given threshold.
    Edge case: Threshold is lower than all counts, no k-mers to remove.
    """
    kct = setup_kmer_table

    # Set a threshold that only removes k-mers with counts > 1 (GGGG)
    removed = kct.maxcut(2)
    assert removed == 1, "Expected 'CCCC/GGGG' to be removed."
    assert kct.get("GGGG") == 0, "Expected 'CCCC/GGGG' to be removed."
    assert kct.get("AAAA") == 2, "Should not remove kmers with exact maxcut value, only greater."

    # Edge case: Threshold is higher than all k-mer counts (remove none)
    removed = kct.maxcut(10)
    assert removed == 0, "Expected no k-mers to be removed since all counts are < 10."
    assert len(kct.hashes) == 2, "Expected 2 records with counts < 10 to remain in the table."
    
    # Edge case: Threshold is lower than all k-mer counts (remove all)
    removed = kct.maxcut(0)
    assert removed == 2, "Expected no k-mers to be removed since all counts are > 0."
    assert len(kct.hashes) == 0, "Expected 0 records with counts < 1 to remain in the table."