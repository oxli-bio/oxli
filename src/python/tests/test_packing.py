"""Tests for the 2-bit packed canonical-k-mer storage (``store_kmers=True``).

The packed representation is an internal detail, but it must decode back to
exactly the canonical strings the previous ``String`` storage produced. These
tests exercise a range of k-mer sizes — including values that do and do not fall
on a 4-base byte boundary — plus reverse-complement canonicalisation.
"""

import pytest

from oxli import KmerCountTable


def _canonical(kmer: str) -> str:
    """Reference canonical form: lexicographically smaller of kmer and its RC."""
    complement = str.maketrans('ACGT', 'TGCA')
    rc = kmer.translate(complement)[::-1]
    return min(kmer, rc)


# k values chosen to straddle the 4-bases-per-byte packing boundary
# (e.g. 4, 8, 16, 32 are exact; 3, 5, 21, 31, 33 have a partial final byte).
@pytest.mark.parametrize('ksize', [3, 4, 5, 8, 15, 16, 21, 31, 32, 33, 64])
def test_unhash_roundtrips_across_ksizes(ksize):
    """Every stored hash unhashes to a canonical k-mer that re-hashes to itself."""
    # A deterministic, non-repetitive sequence a few k-mers long.
    bases = 'ACGTACGTGGCCAATTACGTACGTTTGGCCAAACGTACGTACGTGGCCAATT' * 3
    seq = bases[: ksize + 25]

    table = KmerCountTable(ksize=ksize, store_kmers=True)
    table.consume(seq)

    assert table.hashes, 'expected at least one k-mer stored'
    for h in table.hashes:
        kmer = table.unhash(h)
        assert len(kmer) == ksize, 'decoded k-mer has the wrong length'
        assert set(kmer) <= set('ACGT'), 'decoded k-mer has non-DNA characters'
        # unhash returns the canonical form...
        assert kmer == _canonical(kmer)
        # ...and it re-hashes to the same value (packing is lossless).
        assert table.hash_kmer(kmer) == h


def test_unhash_returns_canonical_reverse_complement():
    """A k-mer and its reverse complement share one canonical stored form."""
    table = KmerCountTable(ksize=4, store_kmers=True)
    # 'GGGG' canonicalises to 'CCCC' (its RC, which is lexicographically smaller).
    table.count('GGGG')
    h = table.hash_kmer('GGGG')
    assert table.hash_kmer('CCCC') == h  # same hash for RC pair
    assert table.unhash(h) == 'CCCC'


def test_dump_kmers_matches_unhash():
    """dump_kmers decodes the packed map consistently with unhash."""
    table = KmerCountTable(ksize=5, store_kmers=True)
    table.consume('ACGTACGTAACCGGTT')

    dumped = dict(table.dump_kmers())
    for h in table.hashes:
        assert table.unhash(h) in dumped
