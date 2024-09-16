import pytest
import tempfile
from os import remove
from oxli import KmerCountTable

@pytest.fixture
def kmer_count_table():
    """Fixture to set up a KmerCountTable instance with sample data."""
    kct = KmerCountTable(ksize=4)
    kct.count("AAAA")  # 17832910516274425539
    kct.count("TTTT")  # 17832910516274425539
    kct.count("AATT")  # 382727017318141683
    kct.count("GGGG")  # 73459868045630124
    kct.count("GGGG")  # 73459868045630124
    return kct

@pytest.fixture
def empty_kmer_count_table():
    """Fixture to set up an empty KmerCountTable instance."""
    return KmerCountTable(ksize=4)

def test_dump_hashes_return_vector(kmer_count_table):
    """Test the dump_hashes function when not writing to a file.

    This test checks if the function returns the correct list of (hash, count) tuples.
    """
    result = kmer_count_table.dump_hashes(file=None, sortkeys=False)
    
    # Expected output sorted by count then hash (default behavior)
    expected = [
        (382727017318141683, 1),   # 'AATT'
        (73459868045630124, 2),    # 'GGGG'
        (17832910516274425539, 2)  # 'AAAA'/'TTTT'
    ]
    
    assert result == expected, f"Expected {expected}, but got {result}"

def test_dump_hashes_write_to_file(kmer_count_table):
    """Test the dump_hashes function when writing to a file.

    This test checks if the function correctly writes the hash:count pairs to a file.
    """
    with tempfile.NamedTemporaryFile(delete=False) as temp_file:
        temp_file_path = temp_file.name
    
    kmer_count_table.dump_hashes(file=temp_file_path, sortkeys=False)
    
    with open(temp_file_path, 'r') as f:
        lines = f.readlines()
    
    # Expected output sorted by count then hash (default behavior)
    expected_lines = [
        f"{382727017318141683}\t1\n",   # 'AATT'
        f"{73459868045630124}\t2\n",    # 'GGGG'
        f"{17832910516274425539}\t2\n"  # 'AAAA'/'TTTT'
    ]
    
    assert lines == expected_lines, f"Expected {expected_lines}, but got {lines}"
    
    # Cleanup
    remove(temp_file_path)

def test_dump_hashes_sortkeys(kmer_count_table):
    """Test the dump_hashes function with sortkeys=True.

    This test verifies if the function sorts by hash keys when `sortkeys` is set to True.
    """
    result = kmer_count_table.dump_hashes(file=None, sortkeys=True)
    
    # Expected output sorted by hash key
    expected = [
        (73459868045630124, 2),    # 'GGGG'
        (382727017318141683, 1),   # 'AATT'
        (17832910516274425539, 2)  # 'AAAA'/'TTTT'        
    ]
    
    assert result == expected, f"Expected {expected}, but got {result}"

def test_dump_hash_empty_table(empty_kmer_count_table):
    """Test the dump_hashes function on an empty KmerCountTable.

    This test checks that the function handles an empty table correctly.
    """
    # Test that calling dump_hashes without file returns an empty list
    result = empty_kmer_count_table.dump_hashes(file=None, sortkeys=False)
    assert result == [], "Expected an empty list from an empty KmerCountTable"
    
    # Test that calling dump_hashes with a file writes nothing to the file
    with tempfile.NamedTemporaryFile(delete=False) as temp_file:
        temp_file_path = temp_file.name
    
    empty_kmer_count_table.dump_hashes(file=temp_file_path, sortkeys=False)
    
    with open(temp_file_path, 'r') as f:
        lines = f.readlines()
    
    assert lines == [], "Expected an empty file for an empty KmerCountTable"
    
    # Cleanup
    remove(temp_file_path)