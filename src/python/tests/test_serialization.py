import gzip
import json
import tempfile
import pytest

from oxli import KmerCountTable
from os import remove
from test_attr import get_version_from_cargo_toml

CURRENT_VERSION = get_version_from_cargo_toml()

@pytest.fixture
def sample_kmer_table():
    """Fixture that provides a sample KmerCountTable object."""
    table = KmerCountTable(ksize=4)
    table.count("AAAA")
    table.count("TTTT")
    return table

@pytest.fixture
def temp_file():
    """Fixture that provides a temporary file path for testing."""
    with tempfile.NamedTemporaryFile(delete=False, suffix='.json.gz') as temp:
        yield temp.name
    # Remove the file after the test is done
    remove(temp.name)

def test_serialize_json(sample_kmer_table):
    """
    Test case for the `serialize_json` function.
    
    This test verifies that the `serialize_json` function correctly serializes a
    KmerCountTable object into a JSON string.
    """
    # Serialize the KmerCountTable object to JSON
    json_data = sample_kmer_table.serialize_json()

    # Convert back to dict to verify correctness
    json_dict = json.loads(json_data)

    # Check that essential attributes exist
    assert "counts" in json_dict, "Counts should be serialized."
    assert json_dict["ksize"] == 4, "Ksize should be correctly serialized."
    assert sample_kmer_table.version == json_dict["version"], "Version should be serialized."

def test_save_load_roundtrip(sample_kmer_table, temp_file):
    """
    Test the save and load functionality.

    This test saves a KmerCountTable object to a file, then loads it back and
    verifies that the data in the loaded object matches the original.
    """
    # Save the sample KmerCountTable to a Gzip file
    sample_kmer_table.save(temp_file)

    # Load the KmerCountTable from the file
    loaded_table = KmerCountTable.load(temp_file)

    # Verify that the loaded data matches the original
    assert loaded_table.get("AAAA") == sample_kmer_table.get("AAAA"), "Counts should be preserved after loading."
    assert loaded_table.get("TTTT") == sample_kmer_table.get("TTTT"), "Counts for reverse complement should be preserved."
    assert list(loaded_table) == list(sample_kmer_table), "All records in same order."
    
def test_version_warning_on_load_stderr(sample_kmer_table, temp_file, capfd):
    """
    Test that a warning is issued if the loaded object's version is different from the current Oxli version.

    Uses pytest's capsys fixture to capture stderr output.
    """
    # Save the table to a file
    sample_kmer_table.save(temp_file)

    # Mock the current version to simulate a version mismatch
    mock_json = sample_kmer_table.serialize_json().replace(CURRENT_VERSION, "0.0.1")
    with gzip.open(temp_file, 'wt') as f:
        json.dump(json.loads(mock_json), f)

    # Capture stderr output
    loaded_table = KmerCountTable.load(temp_file)
    captured = capfd.readouterr()

    # Check stderr for the version mismatch warning
    assert "Version mismatch" in captured.err
    assert f"loaded version is 0.0.1, but current version is {CURRENT_VERSION}" in captured.err
