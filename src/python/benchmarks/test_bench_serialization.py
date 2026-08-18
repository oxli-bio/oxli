"""Benchmarks for serialization, persistence, and dumping."""

import oxli


def test_serialize_json(benchmark, populated_table):
    """serialize_json() the whole table to a JSON string."""
    benchmark(lambda: populated_table.serialize_json())


def test_save(benchmark, populated_table, tmp_path):
    """save() the table to a gzip-compressed binary (bincode) file."""
    out = str(tmp_path / 'bench_save.oxli')

    def run():
        populated_table.save(out)

    benchmark(run)


def test_load(benchmark, populated_table, tmp_path):
    """load() a table back from a saved gzip-compressed binary file."""
    path = str(tmp_path / 'bench_load.oxli')
    populated_table.save(path)

    def run():
        oxli.KmerCountTable.load(path)

    benchmark(run)


def test_dump(benchmark, populated_table):
    """dump() (hash, count) pairs, sorted by count."""
    benchmark(lambda: populated_table.dump(sortcounts=True))


def test_dump_kmers(benchmark, stored_table):
    """dump_kmers() (canonical_kmer, count) pairs from a store_kmers table."""
    benchmark(lambda: stored_table.dump_kmers(sortcounts=True))
