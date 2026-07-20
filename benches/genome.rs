/// Benchmarks for oxli using the *Akkermansia muciniphila* ATCC BAA-835 genome
/// fragment bundled with the repository (`doc/example.fa`, ~350 kb).
///
/// The sequence is loaded from the repository at benchmark time, so the
/// benchmarks are fully deterministic and require no network access. This keeps
/// CodSpeed measurements stable and reproducible in CI.
///
/// # Strategy comparison
///
/// Two strategies are benchmarked:
///
/// * **`consume`** – single-threaded sequential k-mer counting.  Simple and
///   cache-friendly; best for small inputs or environments with a single CPU.
///
/// * **`parallel_consume`** – the sequence is split into overlapping chunks
///   (overlap = ksize − 1 so boundary k-mers are never missed).  Each chunk
///   is processed independently by a Rayon worker thread, and the resulting
///   per-chunk tables are merged serially into the main table.  Scales well
///   with available CPU cores; preferable for long sequences on multi-core
///   hardware.
///
/// Best practice: use `parallel_consume` when the sequence is large enough to
/// benefit from parallelism (typically >> `chunk_size`).  The default
/// `chunk_size` of 50 000 bases balances thread-dispatch overhead against
/// per-chunk work.  Smaller chunks increase parallelism but also increase the
/// overhead of table creation and serial merging; larger chunks reduce
/// parallelism but lower overhead.
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use oxli::KmerCountTable;
use std::fs::File;
use std::hint::black_box;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

// ── genome helpers ────────────────────────────────────────────────────────────

/// Path to the FASTA file bundled with the repository.
fn genome_path() -> PathBuf {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("doc");
    path.push("example.fa");
    path
}

/// Read every DNA sequence line from the bundled FASTA file and return them
/// concatenated as a single `String`.
fn load_genome_sequence() -> String {
    let path = genome_path();
    let file = File::open(&path)
        .unwrap_or_else(|e| panic!("[bench] Cannot open genome file {path:?}: {e}"));
    let reader = BufReader::new(file);

    let mut sequence = String::new();
    for line in reader.lines() {
        let line = line.expect("[bench] Error reading genome");
        if !line.starts_with('>') {
            sequence.push_str(line.trim());
        }
    }

    assert!(
        !sequence.is_empty(),
        "[bench] Genome sequence is empty after parsing."
    );

    sequence
}

// ── benchmark functions ───────────────────────────────────────────────────────

/// Benchmark `KmerCountTable::consume` (single-threaded) for a range of k sizes.
fn bench_consume(c: &mut Criterion) {
    let sequence = load_genome_sequence();

    let mut group = c.benchmark_group("consume");
    for ksize in [21u8, 31] {
        group.bench_with_input(BenchmarkId::new("akkermansia", ksize), &ksize, |b, &k| {
            b.iter(|| {
                let mut table = KmerCountTable::new(k, false);
                table
                    .consume(black_box(&sequence), true)
                    .expect("consume failed");
            });
        });
    }
    group.finish();
}

/// Benchmark `KmerCountTable::parallel_consume` with the default chunk size.
fn bench_parallel_consume(c: &mut Criterion) {
    let sequence = load_genome_sequence();

    let mut group = c.benchmark_group("parallel_consume");
    for ksize in [21u8, 31] {
        group.bench_with_input(BenchmarkId::new("akkermansia", ksize), &ksize, |b, &k| {
            b.iter(|| {
                let mut table = KmerCountTable::new(k, false);
                table
                    .parallel_consume(black_box(&sequence), 50_000, true)
                    .expect("parallel_consume failed");
            });
        });
    }
    group.finish();
}

/// Benchmark `parallel_consume` across different chunk sizes to find the sweet
/// spot for the bundled genome fragment at ksize = 21.
fn bench_parallel_consume_chunk_sizes(c: &mut Criterion) {
    let sequence = load_genome_sequence();

    let mut group = c.benchmark_group("parallel_consume_chunk_sizes");
    for chunk_size in [10_000usize, 50_000, 200_000] {
        group.bench_with_input(
            BenchmarkId::new("akkermansia_k21", chunk_size),
            &chunk_size,
            |b, &cs| {
                b.iter(|| {
                    let mut table = KmerCountTable::new(21, false);
                    table
                        .parallel_consume(black_box(&sequence), cs, true)
                        .expect("parallel_consume failed");
                });
            },
        );
    }
    group.finish();
}

// ── criterion entry points ────────────────────────────────────────────────────

criterion_group!(
    benches,
    bench_consume,
    bench_parallel_consume,
    bench_parallel_consume_chunk_sizes
);
criterion_main!(benches);