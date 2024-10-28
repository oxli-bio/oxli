import pytest

import oxli


def test_basic():
    "string containing only forward canonical kmers."
    seq = "ATAAACC"  # all forward k-mers
    cg = oxli.KmerCountTable(ksize=4)

    x = cg.kmers_and_hashes(seq, False)
    assert x == [
        ("ATAA", 179996601836427478),
        ("TAAA", 15286642655859448092),
        ("AAAC", 9097280691811734508),
        ("AACC", 6779379503393060785),
    ]


def test_basic_rc():
    "string containing only reverse canonical kmers."
    seq = "GGTTTAT"
    cg = oxli.KmerCountTable(ksize=4)

    x = cg.kmers_and_hashes(seq, False)
    print(x)
    assert x == [
        ("AACC", 6779379503393060785),
        ("AAAC", 9097280691811734508),
        ("TAAA", 15286642655859448092),
        ("ATAA", 179996601836427478),
    ]


def test_basic_mixed():
    "string containing forward and reverse canonical kmers."
    seq = "ACGTTG"
    cg = oxli.KmerCountTable(ksize=4)

    x = cg.kmers_and_hashes(seq, False)
    print(x)
    assert x == [
        ("ACGT", 2597925387403686983),
        ("AACG", 7952982457453691616),
        ("CAAC", 7315150081962684964),
    ]

    for kmer, hashval_rs in x:
        print(kmer, hashval_rs, cg.hash_kmer(kmer))
        assert cg.hash_kmer(kmer) == hashval_rs


def test_basic_lower():
    "Test that sequences are turned into uppercase appropriately."
    seq = "acgttg"
    cg = oxli.KmerCountTable(ksize=4)

    x = cg.kmers_and_hashes(seq, False)
    print(x)
    assert x == [
        ("ACGT", 2597925387403686983),
        ("AACG", 7952982457453691616),
        ("CAAC", 7315150081962684964),
    ]


def test_bad_kmers_raise_warning(capfd):
    "Test that bad k-mers print warning with info"
    seq = "acxttg"
    cg = oxli.KmerCountTable(ksize=4)

    # Capture stderr output
    x = cg.kmers_and_hashes(seq, False)
    captured = capfd.readouterr()

    # Check for warning in stderr
    assert f"bad k-mer at position 1: ACXT" in captured.err


def test_bad_kmers_raise_warning_2(capfd):
    "Test bad k-mers raise the right error even when not at beginning :)"
    seq = "aattxttgg"
    cg = oxli.KmerCountTable(ksize=4)

    # Capture stderr output
    x = cg.kmers_and_hashes(seq, False)
    captured = capfd.readouterr()

    # Check for warning in stderr
    assert f"bad k-mer at position 2: ATTX" in captured.err


def test_report_bad_kmers():
    "Test that bad k-mers are reported as (" ",0) when skip_bad_kmers is False"
    seq = "aattxttgg"
    cg = oxli.KmerCountTable(ksize=4)

    x = cg.kmers_and_hashes(seq, False)
    print(x)
    assert x == [
        ("AATT", 382727017318141683),
        ("", 0),
        ("", 0),
        ("", 0),
        ("", 0),
        ("CCAA", 1798905482136869687),
    ]


def test_skip_bad_kmers():
    "Test that bad k-mers are ommited when skip_bad_kmers is True"
    seq = "aattxttgg"
    cg = oxli.KmerCountTable(ksize=4)

    x = cg.kmers_and_hashes(seq, True)
    print(x)
    assert x == [
        ("AATT", 382727017318141683),
        ("CCAA", 1798905482136869687),
    ]


# Tests for hash:kmer storage and retreival


def test_count_saves_kmer():
    """Test that count() stores k-mers and their corresponding hashes when store_kmers=True."""
    kmer = "AAAA"
    cg = oxli.KmerCountTable(ksize=4, store_kmers=True)

    # Call count() on a k-mer
    count = cg.count(kmer)

    # Check that the k-mer was counted
    assert count == 1, f"Expected count to be 1 after first insertion, but got {count}"

    # Hash value of the k-mer should now exist in the hash_to_kmer map
    hashval = cg.hash_kmer(kmer)

    # Check that the k-mer is stored correctly in the hash_to_kmer map
    stored_kmer = cg.unhash(hashval)
    assert (
        stored_kmer == kmer
    ), f"Expected stored k-mer to be {kmer}, but got {stored_kmer}"


def test_count_saves_canonical_kmer():
    """Test that count() stores correct canonical form of k-mers and their corresponding hashes when store_kmers=True."""
    cg = oxli.KmerCountTable(ksize=4, store_kmers=True)
    kmer = "TTTT"
    canon_kmer = "AAAA"

    # Call count() on a k-mer
    cg.count(kmer)

    # Hash value of the k-mer should now exist in the hash_to_kmer map
    hashval = cg.hash_kmer(kmer)

    # Check that the k-mer is stored correctly in the hash_to_kmer map
    stored_kmer = cg.unhash(hashval)

    assert (
        stored_kmer == canon_kmer
    ), f"Expected stored k-mer to be {canon_kmer}, but got {stored_kmer}"


def test_consume_saves_kmers():
    """Test that consume() processes a sequence and stores k-mers and their hashes."""
    seq = "ACGTTG"
    cg = oxli.KmerCountTable(ksize=4, store_kmers=True)

    # Consume the sequence, expecting 3 k-mers ("ACGT", "AACG", "CAAC")
    n_kmers = cg.consume(seq)

    # Check that 3 k-mers were processed
    assert n_kmers == 3, f"Expected to consume 3 k-mers, but got {n_kmers}"

    # Check that all k-mers are stored in the hash_to_kmer map
    for kmer in ["ACGT", "AACG", "CAAC"]:
        hashval = cg.hash_kmer(kmer)
        stored_kmer = cg.unhash(hashval)
        assert (
            stored_kmer == kmer
        ), f"Expected stored k-mer to be {kmer}, but got {stored_kmer}"


def test_count_increments_kmer():
    """Test that count() increments the count of a k-mer when called multiple times."""
    kmer = "AAAA"
    rev_kmer = "TTTT"
    cg = oxli.KmerCountTable(ksize=4, store_kmers=True)

    # Call count() twice on the same k-mer
    count1 = cg.count(kmer)
    count2 = cg.count(rev_kmer)

    # Check that the count has incremented
    assert (
        count1 == 1
    ), f"Expected count to be 1 after first insertion, but got {count1}"
    assert (
        count2 == 2
    ), f"Expected count to be 2 after second insertion, but got {count2}"

    # Ensure the k-mer is still stored correctly in hash_to_kmer
    hashval = cg.hash_kmer(kmer)
    stored_kmer = cg.unhash(hashval)
    assert (
        stored_kmer == kmer
    ), f"Expected stored k-mer to be {kmer}, but got {stored_kmer}"


def test_consume_increments_kmers():
    """Test that consume() increments k-mer counts when the same k-mers are encountered."""
    sequence = "AAAAACCCC"  # Contains overlapping "AAAA" twice
    cg = oxli.KmerCountTable(ksize=4, store_kmers=True)

    # Consume the sequence, expecting 6 k-mers (AAAA, AAAA, AAAC, AACC, ACCC, CCCC)
    n_kmers = cg.consume(sequence)

    # Check that 6 k-mers were processed
    assert n_kmers == 6, f"Expected to consume 6 k-mers, but got {n_kmers}"

    # Check that the count for "AAAA" is now 2
    assert cg.get("AAAA") == 2, "Expected count for 'AAAA' to be 2"


def test_unhash_invalid_kmer():
    """Test that unhash() raises an error when given an invalid hash."""
    cg = oxli.KmerCountTable(ksize=4, store_kmers=True)
    cg.count("AAAA")

    invalid_hash = 1234567890  # A hash that doesn't exist

    # Expecting an exception when trying to unhash an invalid value
    with pytest.raises(
        KeyError, match=f"Warning: Hash {invalid_hash} not found in table."
    ):
        cg.unhash(invalid_hash)


def test_unhash_no_kmer_table():
    """Test that unhash() raises an error when used on a count table without kmer tracking."""
    cg = oxli.KmerCountTable(ksize=3, store_kmers=False)
    kmer = "AAA"
    cg.count(kmer)

    real_hash = cg.hash_kmer(kmer)

    # Expecting an exception when trying to unhash an invalid value
    with pytest.raises(ValueError, match="K-mer storage is not enabled."):
        cg.unhash(real_hash)


def test_consume_invalid_kmers(capfd):
    """Test that consume() processes a sequence and stores k-mers and their hashes."""
    seq = "XAAAAAXGGGG"
    cg = oxli.KmerCountTable(ksize=3, store_kmers=True)

    # Consume the sequence, expecting 5 k-mers ("AAA", "AAA", "AAA", "GGG", "GGG")
    n_kmers = cg.consume(seq)  # [(10679328328772601858, 3), (12126843654075378313, 2)]
    # Capture stderr warnings for bad kmers
    captured = capfd.readouterr()

    # Check for warnings in stderr
    assert "bad k-mer at position 1: XAA" in captured.err
    assert "bad k-mer at position 5: AAX" in captured.err
    assert "bad k-mer at position 6: AXG" in captured.err
    assert "bad k-mer at position 7: XGG" in captured.err

    # Check that 5 k-mers were processed
    assert n_kmers == 5, f"Expected to consume 2 k-mers, but got {n_kmers}"

    # Check 2 distinct kmers
    assert len(cg) == 2, "Expected exactly 2 distinct kmers"

    # Check that all k-mers are stored in the hash_to_kmer map
    for kmer in ["AAA", "CCC"]:
        hashval = cg.hash_kmer(kmer)
        stored_kmer = cg.unhash(hashval)
        assert (
            stored_kmer == kmer
        ), f"Expected stored k-mer to be {kmer}, but got {stored_kmer}"
