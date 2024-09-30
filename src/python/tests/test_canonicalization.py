import pytest

import oxli

def test_canon_kmer():
    """Test the canon() function to ensure it returns the lexicographically smaller k-mer."""
    kmer_table = oxli.KmerCountTable(ksize=4, store_kmers=True)
    
    # Test k-mer with reverse complement
    assert kmer_table.canon("AAAA") == "AAAA", "Expected canonical form to be 'AAAA'"
    assert kmer_table.canon("TTTT") == "AAAA", "Expected canonical form to be 'AAAA'"
    assert kmer_table.canon("ATCG") == "ATCG", "Expected canonical form to be 'CGAT'"
    assert kmer_table.canon("CGAT") == "ATCG", "Expected canonical form to be 'CGAT'"

def test_count_with_canonical_kmer():
    """Test the count() function to ensure it stores the canonical k-mer."""
    kmer_table = oxli.KmerCountTable(ksize=4, store_kmers=True)
    kmer = "TTTT"
    # Count a k-mer and its reverse complement
    kmer_table.count(kmer)
    kmer_table.count(kmer)

    # Check that the canonical k-mer is stored
    hashval = kmer_table.hash_kmer(kmer)
    assert kmer_table.unhash(hashval) == "AAAA", "Expected canonical k-mer 'AAAA'"

    # Check that the count for the canonical k-mer is correct (should be 2)
    assert kmer_table.get_hash(hashval) == 2, "Expected count of 2 for k-mer 'AAAA'"
    

def test_canon_invalid_kmer_size():
    """
    Test that canon() raises a ValueError when the k-mer length does not match the expected ksize.
    """
    kmer_table = oxli.KmerCountTable(ksize=4, store_kmers=True)  # Create a KmerCountTable with ksize=4
    
    # K-mer too short
    with pytest.raises(ValueError, match="kmer size does not match count table ksize"):
        kmer_table.canon("AAA")  # 3-mer for a 4-mer table should raise an error

    # K-mer too long
    with pytest.raises(ValueError, match="kmer size does not match count table ksize"):
        kmer_table.canon("AAAAA")  # 5-mer for a 4-mer table should raise an error


def test_canon_invalid_dna_characters():
    """
    Test that canon() raises a ValueError when the k-mer contains non-DNA characters.
    """
    kmer_table = oxli.KmerCountTable(ksize=4, store_kmers=True)  # Create a KmerCountTable with ksize=4
    
    # Test lowercase conversion
    canon_g= kmer_table.canon("gggg")
    assert canon_g == "CCCC", "Lowercase gggg should be converted to CCCC"
    
    # K-mer with non-DNA character 'X'
    with pytest.raises(ValueError, match="kmer contains invalid characters"):
        kmer_table.canon("ATXG")  # Invalid character 'X' should raise an error

    # K-mer with lowercase and invalid character 'B'
    with pytest.raises(ValueError, match="kmer contains invalid characters"):
        kmer_table.canon("aTbG")  # Lowercase is fine, but 'b' is not a valid DNA character
