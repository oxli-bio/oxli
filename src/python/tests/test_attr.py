from pathlib import Path

import pytest
import toml

import oxli
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


def get_version_from_cargo_toml():
    # Path to Cargo.toml relative to the location of the test file
    cargo_toml_path = Path(__file__).resolve().parents[3] / "Cargo.toml"

    if not cargo_toml_path.exists():
        raise FileNotFoundError(f"{cargo_toml_path} not found")

    with cargo_toml_path.open("r") as f:
        cargo_toml = toml.load(f)

    return cargo_toml["package"]["version"]


def test_kmer_count_table_version():
    # Create an instance of KmerCountTable with a k-mer size
    kmer_table = oxli.KmerCountTable(ksize=31)

    # Get the expected version from Cargo.toml
    expected_version = get_version_from_cargo_toml()

    # Check if the version attribute matches the expected version
    assert (
        kmer_table.version == expected_version
    ), f"Expected version {expected_version}, but got {kmer_table.version}"


# Test consumed bases tracker
def test_initial_consumed():
    kmer_table = oxli.KmerCountTable(ksize=31)
    assert kmer_table.consumed == 0, "Initial consumed should be 0"


def test_consumed_after_count():
    kmer_table = oxli.KmerCountTable(ksize=16)
    kmer_table.count("ACGTACGTACGTACGT")  # Length is 16
    assert (
        kmer_table.consumed == 16
    ), "consumed should be updated to 16 after counting k-mer"


def test_consumed_after_consume():
    kmer_table = oxli.KmerCountTable(ksize=16)
    kmer_table.consume("ACGTACGXACGTACGT", allow_bad_kmers=True)  # Length is 16
    assert (
        kmer_table.consumed == 16
    ), "consumed should be updated to 16 after consuming sequence"


def test_consumed_after_multiple_operations():
    kmer_table = oxli.KmerCountTable(ksize=16)
    kmer_table.count("ACGTACGTACGTACGT")  # Length is 16
    kmer_table.consume("GCTAGCTAGCTA")  # Length is 12, but no kmers added as > 16
    assert (
        kmer_table.consumed == 28
    ), "consumed should be updated to 28 after multiple operations"


# Test total counts attribute
def test_sum_counts_initial():
    kmer_table = oxli.KmerCountTable(ksize=16)
    assert kmer_table.sum_counts == 0, "Initial sum_counts should be 0"


def test_sum_counts_after_count():
    kmer_table = oxli.KmerCountTable(ksize=16)
    kmer_table.count("ACGTACGTACGTACGT")  # Counts as 1
    assert (
        kmer_table.sum_counts == 1
    ), "sum_counts should be updated to 1 after counting k-mer"


def test_sum_counts_after_consume():
    kmer_table = oxli.KmerCountTable(ksize=16)
    kmer_table.consume("ACGTACGTACGTACGTA")  # Counts as 2 k-mers
    assert (
        kmer_table.sum_counts == 2
    ), "sum_counts should be updated after consuming sequence"


def test_sum_counts_after_multiple_operations():
    kmer_table = oxli.KmerCountTable(ksize=16)
    kmer_table.count("ACGTACGTACGTACGT")  # Counts as 1
    kmer_table.consume("ACGTACGTACGTACGTA")  # Counts as 2 k-mers
    assert (
        kmer_table.sum_counts == 3
    ), "sum_counts should be updated after multiple operations"
