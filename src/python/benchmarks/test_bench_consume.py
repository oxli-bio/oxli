"""Benchmarks for sequence import / k-mer consumption."""

import oxli

KSIZES = [21, 31]


def test_consume_k21(benchmark, example_seq):
    """consume() the example genome at k=21 into a fresh table."""

    def run():
        table = oxli.KmerCountTable(21)
        table.consume(example_seq)

    benchmark(run)


def test_consume_k31(benchmark, example_seq):
    """consume() the example genome at k=31 into a fresh table."""

    def run():
        table = oxli.KmerCountTable(31)
        table.consume(example_seq)

    benchmark(run)


def test_consume_store_kmers(benchmark, example_seq):
    """consume() with store_kmers=True (drives KmersAndHashesIter path)."""

    def run():
        table = oxli.KmerCountTable(21, store_kmers=True)
        table.consume(example_seq)

    benchmark(run)


def test_consume_file_fasta(benchmark, example_fa_path):
    """consume_file() reading the plain FASTA genome via needletail."""

    def run():
        table = oxli.KmerCountTable(21)
        table.consume_file(example_fa_path)

    benchmark(run)


def test_consume_file_gzip(benchmark, example_gz_path):
    """consume_file() reading a gzip-compressed copy (decompression path)."""

    def run():
        table = oxli.KmerCountTable(21)
        table.consume_file(example_gz_path)

    benchmark(run)


def test_parallel_consume_default(benchmark, example_seq):
    """parallel_consume() with the default chunk size (50k)."""

    def run():
        table = oxli.KmerCountTable(21)
        table.parallel_consume(example_seq)

    benchmark(run)


def test_parallel_consume_chunk_10k(benchmark, example_seq):
    """parallel_consume() with a small chunk size (more, smaller chunks)."""

    def run():
        table = oxli.KmerCountTable(21)
        table.parallel_consume(example_seq, 10_000)

    benchmark(run)


def test_parallel_consume_chunk_200k(benchmark, example_seq):
    """parallel_consume() with a large chunk size (fewer, larger chunks)."""

    def run():
        table = oxli.KmerCountTable(21)
        table.parallel_consume(example_seq, 200_000)

    benchmark(run)
