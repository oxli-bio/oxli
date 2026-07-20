"""Type stubs for the oxli k-mer counting extension module.

The runtime implementation lives in the compiled Rust module
(``oxli.oxli``); this stub documents the public :class:`KmerCountTable`
API with type hints and numpydoc-style docstrings for editors and type
checkers.
"""

from typing import Iterator

__all__ = ["KmerCountTable"]

class KmerCountTable:
    """A hash table mapping canonical k-mer hashes to their counts.

    Sequences are decomposed into overlapping k-mers of length ``ksize``.
    Each k-mer is reduced to its canonical form (the lexicographically
    smaller of the k-mer and its reverse complement) and hashed with a
    seeded MurmurHash3 (sourmash-compatible). Counts are stored per hash.

    Parameters
    ----------
    ksize : int
        Length of the k-mers to count. Must be positive.
    store_kmers : bool, optional
        If ``True``, also retain a mapping of hash to canonical k-mer
        string so that hashes can be reversed with :meth:`unhash` and
        listed with :meth:`dump_kmers`. Increases memory use. Default is
        ``False``.

    Attributes
    ----------
    ksize : int
        The k-mer size the table was created with.
    version : str
        The oxli version string the table was created with.
    consumed : int
        Total number of sequence bases processed by ``count``/``consume``
        family methods.
    hashes : list of int
        All hash keys currently stored in the table.
    min : int
        Smallest count in the table.
    max : int
        Largest count in the table.
    sum_counts : int
        Sum of all counts in the table.

    Examples
    --------
    >>> import oxli
    >>> kct = oxli.KmerCountTable(ksize=4)
    >>> kct.consume("ACGTACGT")
    5
    >>> kct.get("ACGT")
    2
    """

    ksize: int

    def __init__(self, ksize: int, store_kmers: bool = ...) -> None: ...

    # ── properties ────────────────────────────────────────────────────────────
    @property
    def version(self) -> str:
        """str: oxli version the table was created with."""
        ...

    @property
    def consumed(self) -> int:
        """int: Total bases processed by count/consume methods."""
        ...

    @property
    def hashes(self) -> list[int]:
        """list of int: All hash keys stored in the table."""
        ...

    @property
    def min(self) -> int:
        """int: Minimum count value in the table."""
        ...

    @property
    def max(self) -> int:
        """int: Maximum count value in the table."""
        ...

    @property
    def sum_counts(self) -> int:
        """int: Sum of all counts in the table."""
        ...

    # ── hashing ───────────────────────────────────────────────────────────────
    def hash_kmer(self, kmer: str) -> int:
        """Compute the canonical hash of a single k-mer.

        Parameters
        ----------
        kmer : str
            A k-mer of length ``ksize``.

        Returns
        -------
        int
            The MurmurHash3 hash of the canonical k-mer.

        Raises
        ------
        ValueError
            If ``len(kmer)`` does not equal ``ksize``.
        """
        ...

    def unhash(self, hash: int) -> str:
        """Return the canonical k-mer stored for a hash.

        Only available when the table was created with ``store_kmers=True``.

        Parameters
        ----------
        hash : int
            A hash value present in the table.

        Returns
        -------
        str
            The canonical k-mer string mapped to ``hash``.

        Raises
        ------
        ValueError
            If the table was not created with ``store_kmers=True``.
        KeyError
            If ``hash`` is not present in the table.
        """
        ...

    def canon(self, kmer: str) -> str:
        """Return the canonical form of a k-mer.

        The canonical form is the lexicographically smaller of the k-mer
        and its reverse complement.

        Parameters
        ----------
        kmer : str
            A k-mer of length ``ksize`` containing only ``A``/``C``/``G``/``T``
            (case-insensitive).

        Returns
        -------
        str
            The canonical (upper-case) k-mer.

        Raises
        ------
        ValueError
            If the length does not match ``ksize`` or a non-DNA base is present.
        """
        ...

    # ── counting ──────────────────────────────────────────────────────────────
    def count_hash(self, hashval: int) -> int:
        """Increment the count for a hash value by one.

        Parameters
        ----------
        hashval : int
            The hash whose count should be incremented.

        Returns
        -------
        int
            The new count for ``hashval``.
        """
        ...

    def count(self, kmer: str) -> int:
        """Increment the count of a single k-mer by one.

        Parameters
        ----------
        kmer : str
            A k-mer of length ``ksize``.

        Returns
        -------
        int
            The new count for the k-mer.

        Raises
        ------
        ValueError
            If ``len(kmer)`` does not equal ``ksize``.
        """
        ...

    def consume(self, seq: str, skip_bad_kmers: bool = ...) -> int:
        """Count all k-mers in a DNA string.

        Parameters
        ----------
        seq : str
            A DNA sequence.
        skip_bad_kmers : bool, optional
            If ``True`` (default), k-mers containing non-DNA characters are
            silently skipped. If ``False``, the first bad k-mer raises an
            error.

        Returns
        -------
        int
            The number of k-mers consumed.

        Raises
        ------
        ValueError
            If ``skip_bad_kmers`` is ``False`` and a bad k-mer is encountered.
        """
        ...

    def consume_bytes(self, seq: bytes, skip_bad_kmers: bool = ...) -> int:
        """Count all k-mers in a sequence given as raw bytes.

        Lower-level counterpart to :meth:`consume` used internally by
        :meth:`consume_file`.

        Parameters
        ----------
        seq : bytes
            A DNA sequence as ASCII bytes.
        skip_bad_kmers : bool, optional
            See :meth:`consume`. Default is ``True``.

        Returns
        -------
        int
            The number of k-mers consumed.

        Raises
        ------
        ValueError
            If ``skip_bad_kmers`` is ``False`` and a bad k-mer is encountered,
            or if the bytes are not valid UTF-8 while ``store_kmers`` is enabled.
        """
        ...

    def consume_file(self, filename: str, skip_bad_kmers: bool = ...) -> int:
        """Count all k-mers from a FASTA/FASTQ file.

        The file is parsed with needletail, which auto-detects the format
        and transparently decompresses gzip/bzip2/xz inputs. Each record is
        normalized (upper-cased, line breaks stripped) before counting.

        Parameters
        ----------
        filename : str
            Path to a FASTA or FASTQ file, optionally compressed.
        skip_bad_kmers : bool, optional
            See :meth:`consume`. Default is ``True``.

        Returns
        -------
        int
            Total number of k-mers consumed across all records.

        Raises
        ------
        OSError
            If the file cannot be opened.
        ValueError
            If a record is malformed, or ``skip_bad_kmers`` is ``False`` and a
            bad k-mer is encountered.
        """
        ...

    def parallel_consume(
        self, seq: str, chunk_size: int = ..., skip_bad_kmers: bool = ...
    ) -> int:
        """Count k-mers in a DNA string in parallel using overlapping chunks.

        The sequence is split into chunks overlapping by ``ksize - 1`` bases so
        that no k-mer spanning a boundary is missed; chunks are processed
        concurrently and merged. The result is identical to :meth:`consume`.

        Parameters
        ----------
        seq : str
            A DNA sequence.
        chunk_size : int, optional
            Target number of k-mers per chunk (clamped to at least ``ksize``).
            Default is ``50000``.
        skip_bad_kmers : bool, optional
            See :meth:`consume`. Default is ``True``.

        Returns
        -------
        int
            The number of k-mers consumed.
        """
        ...

    def kmers_and_hashes(
        self, seq: str, skip_bad_kmers: bool
    ) -> list[tuple[str, int]]:
        """Return the canonical k-mer and hash for each window of a sequence.

        Parameters
        ----------
        seq : str
            A DNA sequence.
        skip_bad_kmers : bool
            If ``True``, bad k-mers yield ``("", 0)``; if ``False``, the first
            bad k-mer raises an error.

        Returns
        -------
        list of tuple of (str, int)
            One ``(canonical_kmer, hash)`` pair per k-mer window.

        Raises
        ------
        ValueError
            If ``skip_bad_kmers`` is ``False`` and a bad k-mer is encountered.
        """
        ...

    # ── retrieval ─────────────────────────────────────────────────────────────
    def get(self, kmer: str) -> int:
        """Return the stored count for a k-mer.

        Parameters
        ----------
        kmer : str
            A k-mer of length ``ksize``.

        Returns
        -------
        int
            The count for the canonical form of ``kmer`` (``0`` if absent).

        Raises
        ------
        ValueError
            If ``len(kmer)`` does not equal ``ksize``.
        """
        ...

    def get_hash(self, hashval: int) -> int:
        """Return the stored count for a hash value.

        Parameters
        ----------
        hashval : int
            A hash value.

        Returns
        -------
        int
            The count for ``hashval`` (``0`` if absent).
        """
        ...

    def get_hash_array(self, hash_keys: list[int]) -> list[int]:
        """Return counts for a list of hash values.

        Parameters
        ----------
        hash_keys : list of int
            Hash values to look up.

        Returns
        -------
        list of int
            Counts in the same order as ``hash_keys`` (``0`` for absent keys).
        """
        ...

    # ── mutation ──────────────────────────────────────────────────────────────
    def drop(self, kmer: str) -> None:
        """Remove a k-mer from the table.

        Parameters
        ----------
        kmer : str
            A k-mer of length ``ksize``.

        Raises
        ------
        ValueError
            If ``len(kmer)`` does not equal ``ksize``.
        """
        ...

    def drop_hash(self, hashval: int) -> None:
        """Remove a hash value from the table.

        Parameters
        ----------
        hashval : int
            The hash value to remove.
        """
        ...

    def mincut(self, min_count: int) -> int:
        """Drop all k-mers with a count below a threshold.

        Parameters
        ----------
        min_count : int
            K-mers with a count strictly less than ``min_count`` are removed.

        Returns
        -------
        int
            The number of k-mers removed.
        """
        ...

    def maxcut(self, max_count: int) -> int:
        """Drop all k-mers with a count above a threshold.

        Parameters
        ----------
        max_count : int
            K-mers with a count strictly greater than ``max_count`` are removed.

        Returns
        -------
        int
            The number of k-mers removed.
        """
        ...

    # ── serialization ─────────────────────────────────────────────────────────
    def serialize_json(self) -> str:
        """Serialize the table to a JSON string.

        Returns
        -------
        str
            A JSON representation of the table.
        """
        ...

    def save(self, filepath: str) -> None:
        """Save the table to a gzip-compressed JSON file.

        Parameters
        ----------
        filepath : str
            Destination path.

        Raises
        ------
        OSError
            If the file cannot be written.
        """
        ...

    @staticmethod
    def load(filepath: str) -> "KmerCountTable":
        """Load a table previously written with :meth:`save`.

        Parameters
        ----------
        filepath : str
            Path to a saved table (compression auto-detected).

        Returns
        -------
        KmerCountTable
            The loaded table.

        Raises
        ------
        OSError
            If the file cannot be read.
        ValueError
            If the file does not contain a valid serialized table.
        """
        ...

    def dump(
        self,
        file: str | None = ...,
        sortcounts: bool = ...,
        sortkeys: bool = ...,
    ) -> list[tuple[int, int]]:
        """Return or write ``(hash, count)`` pairs.

        Parameters
        ----------
        file : str or None, optional
            If given, write tab-separated pairs to this path; otherwise return
            them. Default is ``None``.
        sortcounts : bool, optional
            Sort by count (secondary sort by hash). Default is ``False``.
        sortkeys : bool, optional
            Sort by hash key. Default is ``False``.

        Returns
        -------
        list of tuple of (int, int)
            ``(hash, count)`` pairs.

        Raises
        ------
        ValueError
            If both ``sortcounts`` and ``sortkeys`` are ``True``.
        """
        ...

    def dump_kmers(
        self,
        file: str | None = ...,
        sortcounts: bool = ...,
        sortkeys: bool = ...,
    ) -> list[tuple[str, int]]:
        """Return or write ``(canonical_kmer, count)`` pairs.

        Requires the table to have been created with ``store_kmers=True``.

        Parameters
        ----------
        file : str or None, optional
            If given, write tab-separated pairs to this path; otherwise return
            them. Default is ``None``.
        sortcounts : bool, optional
            Sort by count (secondary sort by k-mer). Default is ``False``.
        sortkeys : bool, optional
            Sort by canonical k-mer. Default is ``False``.

        Returns
        -------
        list of tuple of (str, int)
            ``(canonical_kmer, count)`` pairs.

        Raises
        ------
        ValueError
            If the table was not created with ``store_kmers=True``, or if both
            ``sortcounts`` and ``sortkeys`` are ``True``.
        """
        ...

    def histo(self, zero: bool = ...) -> list[tuple[int, int]]:
        """Return a histogram of k-mer frequencies.

        Parameters
        ----------
        zero : bool, optional
            If ``True`` (default), include frequencies from ``0`` up to the
            maximum observed count, even where no k-mers have that frequency.

        Returns
        -------
        list of tuple of (int, int)
            ``(frequency, number_of_kmers)`` pairs.
        """
        ...

    # ── set operations ────────────────────────────────────────────────────────
    def hash_set(self) -> set[int]:
        """Return the set of hash keys in the table.

        Returns
        -------
        set of int
            All hash keys currently stored.
        """
        ...

    def union(self, other: "KmerCountTable") -> set[int]:
        """Return the union of hash keys with another table.

        Parameters
        ----------
        other : KmerCountTable
            The other table.

        Returns
        -------
        set of int
            Hashes present in either table.
        """
        ...

    def intersection(self, other: "KmerCountTable") -> set[int]:
        """Return the intersection of hash keys with another table.

        Parameters
        ----------
        other : KmerCountTable
            The other table.

        Returns
        -------
        set of int
            Hashes present in both tables.
        """
        ...

    def difference(self, other: "KmerCountTable") -> set[int]:
        """Return hashes present in this table but not another.

        Parameters
        ----------
        other : KmerCountTable
            The other table.

        Returns
        -------
        set of int
            Hashes in ``self`` but not ``other``.
        """
        ...

    def symmetric_difference(self, other: "KmerCountTable") -> set[int]:
        """Return hashes present in exactly one of two tables.

        Parameters
        ----------
        other : KmerCountTable
            The other table.

        Returns
        -------
        set of int
            Hashes in exactly one of ``self`` or ``other``.
        """
        ...

    def add(self, other: "KmerCountTable") -> tuple[int, int]:
        """Add the counts of another table into this one, in place.

        Parameters
        ----------
        other : KmerCountTable
            The table whose counts are added into ``self``.

        Returns
        -------
        tuple of (int, int)
            A summary of the merge (e.g. number of keys updated and added).

        Raises
        ------
        ValueError
            If the two tables have different ``ksize`` values.
        """
        ...

    # ── similarity metrics ────────────────────────────────────────────────────
    def jaccard(self, other: "KmerCountTable") -> float:
        """Return the Jaccard similarity of k-mer sets with another table.

        Parameters
        ----------
        other : KmerCountTable
            The other table.

        Returns
        -------
        float
            The Jaccard index in ``[0, 1]``.
        """
        ...

    def cosine(self, other: "KmerCountTable") -> float:
        """Return the cosine similarity of count vectors with another table.

        Parameters
        ----------
        other : KmerCountTable
            The other table.

        Returns
        -------
        float
            The cosine similarity in ``[0, 1]``.
        """
        ...

    # ── dunder methods ────────────────────────────────────────────────────────
    def __len__(self) -> int: ...
    def __getitem__(self, kmer: str) -> int: ...
    def __setitem__(self, kmer: str, count: int) -> None: ...
    def __iter__(self) -> Iterator[tuple[int, int]]: ...
    def __or__(self, other: "KmerCountTable") -> set[int]: ...
    def __and__(self, other: "KmerCountTable") -> set[int]: ...
    def __sub__(self, other: "KmerCountTable") -> set[int]: ...
    def __xor__(self, other: "KmerCountTable") -> set[int]: ...
