# oxli

Author: C. Titus Brown (@ctb), ctbrown@ucdavis.edu

oxli is a simple Rust library + Python interface for counting k-mers
in genomic sequencing data.

## Documentation

Please see [the API documentation](doc/api.md).

## Is there anything I should know about oxli?

Two things -

First, oxli is channeling
[khmer](https://khmer.readthedocs.io/en/latest/), a package written by
@ctb and many others.  You shouldn't be too surprised to see useful
functionality from khmer making an appearance in oxli.

Second, it's written on top of the
[sourmash](https://sourmash.readthedocs.io/)
[rust library](https://sourmash.readthedocs.io/), and the underlying
code for dealing with sequence data is pretty well tested.

## What's the history here?

The history is a bit convoluted:

* the khmer package was useful for inspecting large collections of
  k-mers, but was hard to maintain and evolve.

* in ~2016 @ctb's lab more or less switched over to developing
  sourmash, which was initially built on a similar tech stack to khmer
  (Python & C++).
  
* at some point, @luizirber rewrote the sourmash C++ code into Rust.

* this forced @ctb to learn Rust to maintain sourmash.

* @ctb then decided he liked Rust an awful lot, and missed some of the
  khmer functionality.
  
* voila, oxli was born.

---

(Sep 2024)
