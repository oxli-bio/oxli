import gzip
import json

import pytest

from oxli import KmerCountTable
from test_attr import get_version_from_cargo_toml

CURRENT_VERSION = get_version_from_cargo_toml()


@pytest.fixture
def sample_kmer_table():
    """Fixture that provides a sample KmerCountTable object."""
    table = KmerCountTable(ksize=4)
    table.count("AAAA")
    table.count("TTTT")
    return table


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
    assert (
        sample_kmer_table.version == json_dict["version"]
    ), "Version should be serialized."


def test_save_load_roundtrip(sample_kmer_table, tmp_path):
    """
    Test the save and load functionality.

    This test saves a KmerCountTable object to a file, then loads it back and
    verifies that the data in the loaded object matches the original.
    """
    temp_file = str(tmp_path / "save.json")

    # Save the sample KmerCountTable to a Gzip file
    sample_kmer_table.save(temp_file)

    # Load the KmerCountTable from the file
    loaded_table = KmerCountTable.load(temp_file)

    # Verify that the loaded data matches the original
    assert loaded_table.get("AAAA") == sample_kmer_table.get(
        "AAAA"
    ), "Counts should be preserved after loading."
    assert loaded_table.get("TTTT") == sample_kmer_table.get(
        "TTTT"
    ), "Counts for reverse complement should be preserved."
    assert list(loaded_table) == list(sample_kmer_table), "All records in same order."


def test_version_warning_on_load_stderr(sample_kmer_table, tmp_path, capfd):
    """
    Test that a warning is issued if the loaded object's version is different from the current Oxli version.

    Uses pytest's capsys fixture to capture stderr output.
    """
    temp_file = str(tmp_path / "save.json")

    # Save the table to a file
    sample_kmer_table.save(temp_file)

    # Mock the current version to simulate a version mismatch
    mock_json = sample_kmer_table.serialize_json().replace(CURRENT_VERSION, "0.0.1")
    with gzip.open(temp_file, "wt") as f:
        json.dump(json.loads(mock_json), f)

    # Capture stderr output
    loaded_table = KmerCountTable.load(temp_file)
    captured = capfd.readouterr()

    # Check stderr for the version mismatch warning
    assert "Version mismatch" in captured.err
    assert (
        f"loaded version is 0.0.1, but current version is {CURRENT_VERSION}"
        in captured.err
    )


def test_load_bad_json(tmp_path, capfd):
    """
    Test that failure happens appropriately when trying to load a bad
    JSON file.
    """
    temp_file = str(tmp_path / "bad.json")

    with open(temp_file, "wt") as fp:
        fp.write("hello, world")

    with pytest.raises(RuntimeError, match="Deserialization error:"):
        tb = KmerCountTable.load(temp_file)


def test_save_bad_path(sample_kmer_table, tmp_path, capfd):
    """
    Test that failure happens appropriately when trying to save to a bad
    location.
    """
    temp_file = str(tmp_path / "noexist" / "save.json")

    with pytest.raises(OSError, match="No such file or directory"):
        sample_kmer_table.save(temp_file)
