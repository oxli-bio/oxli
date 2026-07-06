/// Benchmarks for oxli using the Escherichia coli str. K-12 substr. MG1655 genome
/// (GenBank accession GCF_000005845.2).
///
/// The genome is fetched from the NCBI FTP site on first run and cached in the
/// system temp directory.  If the download fails (e.g. no network access), all
/// benchmarks in this file are silently skipped so that CI still passes.
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
///   with available CPU cores; preferable for long sequences (≥ hundreds of
///   kilobases) on multi-core hardware.
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
use std::process::Command;

// ── genome helpers ────────────────────────────────────────────────────────────

const GENOME_URL: &str = concat!(
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/",
    "GCF_000005845.2_ASM584v2/",
    "GCF_000005845.2_ASM584v2_genomic.fna.gz"
);

const GENOME_FILENAME: &str = "GCF_000005845.2_ASM584v2_genomic.fna.gz";

/// Returns the path where the cached genome file is stored.
fn genome_path() -> PathBuf {
    let mut path = std::env::temp_dir();
    path.push(GENOME_FILENAME);
    path
}

/// Attempt to download the genome with `curl` or `wget`.
/// Returns `true` on success.  Cleans up any partial file on failure.
fn try_download(dest: &PathBuf) -> bool {
    let dest_str = dest.to_string_lossy();

    // Try curl first (widely available on Linux/macOS/CI)
    if let Ok(status) = Command::new("curl")
        .args(["-fsSL", "-o", &dest_str, GENOME_URL])
        .status()
    {
        if status.success() {
            return true;
        }
    }

    // Fall back to wget
    if let Ok(status) = Command::new("wget")
        .args(["-q", "-O", &dest_str, GENOME_URL])
        .status()
    {
        if status.success() {
            return true;
        }
    }

    // Remove any partial file left behind by a failed download attempt
    let _ = std::fs::remove_file(dest);
    false
}

/// Read every DNA sequence line from a (possibly gzip-compressed) FASTA file
/// and return them concatenated as a single `String`.  Returns `None` if the
/// file cannot be opened or the download fails.
fn load_genome_sequence() -> Option<String> {
    let path = genome_path();

    if !path.exists() {
        eprintln!(
            "[bench] Genome not found at {:?}. Attempting download from NCBI…",
            path
        );
        if !try_download(&path) {
            eprintln!("[bench] Download failed. Genome benchmarks will be skipped.");
            return None;
        }
        eprintln!("[bench] Download complete.");
    }

    let file = File::open(&path)
        .map_err(|e| eprintln!("[bench] Cannot open genome file: {e}"))
        .ok()?;
    let reader = BufReader::new(file);
    let (decompressed, _fmt) = niffler::get_reader(Box::new(reader))
        .map_err(|e| eprintln!("[bench] Cannot decompress genome file: {e}"))
        .ok()?;
    let buf = BufReader::new(decompressed);

    let mut sequence = String::new();
    for line in buf.lines() {
        let line = line
            .map_err(|e| eprintln!("[bench] Error reading genome: {e}"))
            .ok()?;
        if !line.starts_with('>') {
            sequence.push_str(line.trim());
        }
    }

    if sequence.is_empty() {
        eprintln!("[bench] Genome sequence is empty after parsing.");
        return None;
    }

    Some(sequence)
}

// ── benchmark functions ───────────────────────────────────────────────────────

/// Benchmark `KmerCountTable::consume` (single-threaded) for a range of k sizes.
fn bench_consume(c: &mut Criterion) {
    let sequence = match load_genome_sequence() {
        Some(s) => s,
        None => return,
    };

    let mut group = c.benchmark_group("consume");
    for ksize in [21u8, 31] {
        group.bench_with_input(BenchmarkId::new("ecoli", ksize), &ksize, |b, &k| {
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
    let sequence = match load_genome_sequence() {
        Some(s) => s,
        None => return,
    };

    let mut group = c.benchmark_group("parallel_consume");
    for ksize in [21u8, 31] {
        group.bench_with_input(BenchmarkId::new("ecoli", ksize), &ksize, |b, &k| {
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
/// spot for the *E. coli* genome (~4.6 Mbp) at ksize = 21.
fn bench_parallel_consume_chunk_sizes(c: &mut Criterion) {
    let sequence = match load_genome_sequence() {
        Some(s) => s,
        None => return,
    };

    let mut group = c.benchmark_group("parallel_consume_chunk_sizes");
    for chunk_size in [10_000usize, 50_000, 200_000, 500_000] {
        group.bench_with_input(
            BenchmarkId::new("ecoli_k21", chunk_size),
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
