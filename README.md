<a href="https://opensource.org/licenses/BSD-3-Clause">
  <img src="https://img.shields.io/badge/License-BSD_3--Clause-blue.svg" align="left" height="20"/>
</a>
<a href="https://app.codspeed.io/oxli-bio/oxli?utm_source=badge">
  <img src="https://img.shields.io/endpoint?url=https://codspeed.io/badge.json" alt="CodSpeed" align="left" height="20"/>
</a> 
<br>

# oxli

oxli is a powerful Rust library with a simple Python interface for counting k-mers
in genomic sequencing data.

Use oxli to bring fast kmer counting and comparison operations to your Python projects!

This library is written on top of the
[sourmash](https://sourmash.readthedocs.io/)
[rust library](https://sourmash.readthedocs.io/), and the underlying
code for dealing with sequence data is well tested.

## Installation

### Quick setup

oxli is
[available on conda-forge for Linux, Mac OS X, and Windows](https://github.com/conda-forge/oxli-feedstock) for Python versions 3.10 through 3.14 (including the free-threaded 3.13t/3.14t builds):

```bash
conda install oxli
```

This will install the oxli library for Python.

### For developers

You can also try building oxli yourself and using it in [development mode](https://github.com/oxli-bio/oxli/wiki/For-Developers):

```bash
# Setup conda development env
mamba env create -f environment.yml -n oxli

# Install oxli in dev mode
pip install -e '.[test]'
```

## Getting Started

See the [the oxli Wiki](https://github.com/oxli-bio/oxli/wiki/Getting-Started) for documentation on the Python API.

### Basic Usage

Initialise a new `KmerCountTable`
```python
# Import oxli
from oxli import KmerCountTable

# Create new count table
kct = KmerCountTable(ksize=4) # Count 4-mers
```
Adding k-mer counts.

```python
# Add single k-mer with count()
kct.count("AAAA")
>>> 1

# Increment count
kct.count("AAAA")
>>> 2

# Forward and Reverse complement counted together
kct.count("TTTT")
>>> 3

# Add many k-mers from a longer sequence with consume
kct.consume("GGGGGGGGGG") # 7 x 4-mers of 'GGGG'
```

Lookup counts by k-mer.

```python
# Retrieve kmer counts
kct.get('GGGG') # Count for GGGG/CCCC
>>> 7
kct.get('AAAA') #Count for AAAA/TTTT
>>> 3
```

Extracting k-mers from files.

```python
# Screed for FASTA/FASTQ parsing
import screed

# Create new table
counts = KmerCountTable(ksize=21)

# Read fasta records and extract k-mers
for record in screed.open('doc/example.fa'):
    counts.consume(record.sequence)
>>> 349910
```

For convenience, `consume_file` reads a FASTA/FASTQ file directly using a fast
native (needletail) parser, so no external parsing library is required. It
transparently handles gzip/bzip2/xz-compressed files:

```python
counts = KmerCountTable(ksize=21)

# Count k-mers from every record in the file (plain or compressed)
counts.consume_file('doc/example.fa')
>>> 349910
```

### Saving and loading

Tables can be persisted to disk and reloaded:

```python
counts.save('counts.oxli')                 # gzip-compressed binary
reloaded = KmerCountTable.load('counts.oxli')
```

`save` writes a compact gzip-compressed binary format. `load` auto-detects the
format, so tables written by older oxli versions (gzip-JSON) still load. If you
need a text representation, `serialize_json()` returns the table as a JSON
string.


## Benchmarking

oxli has two complementary benchmark suites:

- **Rust / criterion** (`benches/genome.rs`) — micro-benchmarks the core
  `consume` / `parallel_consume` hot paths (and `kmers_and_hashes`, `cosine`,
  `add`) against the *Akkermansia muciniphila* genome fragment bundled at
  `doc/example.fa` (deterministic, no network access):

  ```bash
  make bench      # cargo bench
  ```

- **Python / [pytest-codspeed](https://github.com/CodSpeedHQ/pytest-codspeed)**
  (`src/python/benchmarks/`) — benchmarks the full `KmerCountTable` API as called
  from Python (sequence import, hashing, counting, retrieval, set operations,
  similarity metrics, mutation, serialization, histograms) over the committed
  `doc/example.fa`:

  ```bash
  pip install '.[test]'
  make bench-py   # pytest src/python/benchmarks --codspeed
  ```

  The benchmarks live outside the test path, so the normal `make test` run does
  not collect them. On pull requests the `CodSpeed` GitHub Actions workflow runs
  this suite and reports performance changes (requires the CodSpeed app and a
  `CODSPEED_TOKEN` secret to be configured on the repository).


## What's the history here?

First, oxli is channeling
[khmer](https://khmer.readthedocs.io/en/latest/), a package written by
@ctb and many others.  You shouldn't be too surprised to see useful
functionality from khmer making an appearance in oxli.

The khmer package was useful for inspecting large collections of
k-mers, but was hard to maintain and evolve.

In ~2016 @ctb's lab more or less switched over to developing
sourmash, which was initially built on a similar tech stack to khmer
(Python & C++).

At some point, @luizirber rewrote the sourmash C++ code into Rust.

This forced @ctb to learn Rust to maintain sourmash.

@ctb then decided he liked Rust an awful lot, and missed some of the
khmer functionality.

And, voila! oxli was born.

## Authors

* C. Titus Brown (@ctb), ctbrown@ucdavis.edu
* Adam Taranto (@Adamtaranto)
