<a href="https://opensource.org/licenses/BSD-3-Clause">
  <img src="https://img.shields.io/badge/License-BSD_3--Clause-blue.svg" align="left" height="20"/>
</a> 

<a href="https://gitpod.io/#https://github.com/oxli-bio/oxli">
  <img src="https://gitpod.io/button/open-in-gitpod.svg" align="right" height="35"/>
</a>

<br>

# oxli

oxli is a powerful Rust library with a simple Python interface for counting k-mers
in genomic sequencing data.

Use oxli to bring fast kmer counting and comparison operations to you Python projects.

This library is written on top of the
[sourmash](https://sourmash.readthedocs.io/)
[rust library](https://sourmash.readthedocs.io/), and the underlying
code for dealing with sequence data is well tested.

## Installation

### Quick setup

oxli is
[available on conda-forge for Linux, Mac OS X, and Windows](https://github.com/conda-forge/oxli-feedstock) for Python versions 3.10, 3.11, and 3.12:

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

See the [the Oxli Wiki](https://github.com/oxli-bio/oxli/wiki/Getting-Started) for documentation on the Python API.

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

Author: C. Titus Brown (@ctb), ctbrown@ucdavis.edu  
with with miscellaneous features by @Adamtaranto