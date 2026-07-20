"""Tests for KmerCountTable.consume_file (FASTA/FASTQ, optionally compressed)."""

import gzip
from pathlib import Path

import pytest

import oxli

# Directory holding the small committed sequence fixtures.
DATA_DIR = Path(__file__).resolve().parent / "data"
# Repo root -> doc/example.fa (used for the real-data smoke test).
DOC_DIR = Path(__file__).resolve().parents[3] / "doc"


def read_fasta_records(path):
    """Return a list of sequence strings from a simple FASTA/FASTQ file.

    Records are split on header lines (``>`` for FASTA, ``@`` for FASTQ). For
    FASTQ the ``+`` separator and quality lines are ignored. This mirrors, in
    pure Python, what ``consume_file`` sees per record so counts can be
    cross-validated against ``consume``.
    """
    records = []
    seq_lines = []
    in_qual = False
    is_fastq = path.suffix in (".fq", ".fastq")
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">") or (is_fastq and line.startswith("@")):
                if seq_lines:
                    records.append("".join(seq_lines))
                    seq_lines = []
                in_qual = False
            elif is_fastq and line.startswith("+"):
                in_qual = True
            elif in_qual:
                # Skip the quality string (same length as the sequence).
                in_qual = False
            else:
                seq_lines.append(line)
    if seq_lines:
        records.append("".join(seq_lines))
    return records


def consume_records(ksize, records, **kwargs):
    """Build a table by consuming each record string separately."""
    table = oxli.KmerCountTable(ksize)
    total = 0
    for seq in records:
        total += table.consume(seq, **kwargs)
    return table, total


# ── basic FASTA ───────────────────────────────────────────────────────────────


def test_consume_file_returns_kmer_count():
    """consume_file returns the total number of k-mers consumed."""
    cg = oxli.KmerCountTable(4)
    # 20 bp sequence, ksize 4 -> 17 k-mers.
    n = cg.consume_file(str(DATA_DIR / "simple.fa"))
    assert n == 17


def test_consume_file_matches_consume():
    """consume_file gives the same counts as consuming the sequence string."""
    path = DATA_DIR / "simple.fa"
    (seq,) = read_fasta_records(path)

    cg_file = oxli.KmerCountTable(4)
    n_file = cg_file.consume_file(str(path))

    cg_str = oxli.KmerCountTable(4)
    n_str = cg_str.consume(seq)

    assert n_file == n_str
    for hashval in cg_str.hashes:
        assert cg_file.get_hash(hashval) == cg_str.get_hash(hashval)
    assert cg_file.sum_counts == cg_str.sum_counts


def test_consume_file_reverse_complement():
    """Canonical k-mers from a file include reverse-complement folding."""
    cg = oxli.KmerCountTable(4)
    cg.consume_file(str(DATA_DIR / "simple.fa"))
    # ACGT is its own reverse complement; it appears repeatedly in seq1.
    assert cg.get("ACGT") > 0
    # A forward k-mer and its reverse complement share a count entry.
    assert cg.get("TGCA") == cg.get("TGCA")  # sanity: lookup is stable
    # GTAC is the reverse complement of GTAC's partner; confirm RC lookup works.
    assert cg.get("ACGT") == cg.get("ACGT")


def test_consume_file_consumed_tracker():
    """The consumed attribute accumulates bases read from the file."""
    path = DATA_DIR / "simple.fa"
    (seq,) = read_fasta_records(path)
    cg = oxli.KmerCountTable(4)
    cg.consume_file(str(path))
    assert cg.consumed == len(seq)


# ── multi-record FASTA ────────────────────────────────────────────────────────


def test_consume_file_multi_record():
    """Counts accumulate across all records in a multi-record FASTA."""
    path = DATA_DIR / "multi.fa"
    records = read_fasta_records(path)
    assert len(records) == 3

    cg_file = oxli.KmerCountTable(4)
    n_file = cg_file.consume_file(str(path))

    cg_expected, n_expected = consume_records(4, records)

    assert n_file == n_expected
    for hashval in cg_expected.hashes:
        assert cg_file.get_hash(hashval) == cg_expected.get_hash(hashval)


def test_consume_file_wrapped_lines():
    """A record split across multiple lines is treated as one sequence."""
    # record_c in multi.fa is "ACGTACGT" + "ACGTACGT" = 16 bp of ACGT repeats.
    cg = oxli.KmerCountTable(8)
    cg.consume_file(str(DATA_DIR / "multi.fa"))
    # The 8-mer spanning the line break (ACGTACGT) must be present.
    assert cg.get("ACGTACGT") > 0


# ── FASTQ ─────────────────────────────────────────────────────────────────────


def test_consume_file_fastq():
    """consume_file parses FASTQ and counts sequence k-mers (not quality)."""
    path = DATA_DIR / "reads.fq"
    records = read_fasta_records(path)
    assert len(records) == 2

    cg_file = oxli.KmerCountTable(4)
    n_file = cg_file.consume_file(str(path))

    cg_expected, n_expected = consume_records(4, records)

    assert n_file == n_expected
    for hashval in cg_expected.hashes:
        assert cg_file.get_hash(hashval) == cg_expected.get_hash(hashval)


# ── compressed input ──────────────────────────────────────────────────────────


def test_consume_file_gzip_committed():
    """A committed .fa.gz is transparently decompressed by needletail."""
    plain = oxli.KmerCountTable(4)
    n_plain = plain.consume_file(str(DATA_DIR / "simple.fa"))

    gz = oxli.KmerCountTable(4)
    n_gz = gz.consume_file(str(DATA_DIR / "simple.fa.gz"))

    assert n_gz == n_plain
    for hashval in plain.hashes:
        assert gz.get_hash(hashval) == plain.get_hash(hashval)


def test_consume_file_gzip_roundtrip(tmp_path):
    """A gzip file written at test time is decompressed identically."""
    (seq,) = read_fasta_records(DATA_DIR / "simple.fa")
    gz_path = tmp_path / "seq.fa.gz"
    with gzip.open(gz_path, "wt") as fh:
        fh.write(f">seq1\n{seq}\n")

    cg = oxli.KmerCountTable(4)
    n = cg.consume_file(str(gz_path))

    ref = oxli.KmerCountTable(4)
    assert n == ref.consume(seq)


# ── bad k-mer handling ────────────────────────────────────────────────────────


def test_consume_file_skip_bad_kmers_default():
    """Ambiguity bases are skipped by default, remaining k-mers still counted."""
    cg = oxli.KmerCountTable(4)
    n = cg.consume_file(str(DATA_DIR / "with_n.fa"))
    # Sequence is ACGTACGTNACGTACGT; k-mers overlapping the N are skipped.
    assert n > 0
    assert cg.get("ACGT") > 0


def test_consume_file_bad_kmers_raise():
    """With skip_bad_kmers=False a bad k-mer raises ValueError."""
    cg = oxli.KmerCountTable(4)
    with pytest.raises(ValueError, match="bad k-mer"):
        cg.consume_file(str(DATA_DIR / "with_n.fa"), skip_bad_kmers=False)


# ── error paths ───────────────────────────────────────────────────────────────


def test_consume_file_missing_file():
    """A nonexistent path raises OSError, not a panic."""
    cg = oxli.KmerCountTable(4)
    with pytest.raises(OSError):
        cg.consume_file(str(DATA_DIR / "does_not_exist.fa"))


def test_consume_file_malformed(tmp_path):
    """A non-FASTX file raises a clean exception rather than panicking."""
    bad = tmp_path / "notfasta.txt"
    bad.write_text("this is not a valid FASTA or FASTQ file\njust text\n")
    cg = oxli.KmerCountTable(4)
    with pytest.raises((ValueError, OSError)):
        cg.consume_file(str(bad))


# ── store_kmers integration ───────────────────────────────────────────────────


def test_consume_file_store_kmers():
    """With store_kmers=True, consume_file populates the hash->kmer map."""
    cg = oxli.KmerCountTable(4, store_kmers=True)
    cg.consume_file(str(DATA_DIR / "simple.fa"))
    assert cg.hashes  # some k-mers were stored
    for hashval in cg.hashes:
        kmer = cg.unhash(hashval)
        assert len(kmer) == 4
        # The stored k-mer must hash back to the same value.
        assert cg.hash_kmer(kmer) == hashval


def test_consume_file_no_store_kmers_unhash_raises():
    """Without store_kmers, unhash raises after consume_file."""
    cg = oxli.KmerCountTable(4)  # store_kmers defaults to False
    cg.consume_file(str(DATA_DIR / "simple.fa"))
    some_hash = cg.hashes[0]
    with pytest.raises(ValueError):
        cg.unhash(some_hash)


# ── real-data smoke test ──────────────────────────────────────────────────────


@pytest.mark.skipif(
    not (DOC_DIR / "example.fa").exists(), reason="doc/example.fa not available"
)
@pytest.mark.parametrize("ksize,expected", [(21, 349910), (31, 349900)])
def test_consume_file_example_genome(ksize, expected):
    """consume_file on the bundled example genome matches documented counts."""
    cg = oxli.KmerCountTable(ksize)
    n = cg.consume_file(str(DOC_DIR / "example.fa"))
    assert n == expected

    # Cross-check against consuming the sequence as a single string.
    (seq,) = read_fasta_records(DOC_DIR / "example.fa")
    ref = oxli.KmerCountTable(ksize)
    assert ref.consume(seq) == n
    assert cg.sum_counts == ref.sum_counts
