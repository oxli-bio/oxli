"""Fast k-mer counting for genomic sequencing data, powered by Rust.

This package exposes :class:`KmerCountTable`, implemented as a compiled
extension module. The pure-Python ``__init__`` simply re-exports the compiled
symbols so that ``import oxli`` continues to work in a mixed Rust/Python
(maturin) layout, while the accompanying ``__init__.pyi`` stub provides type
hints and documentation for editors and type checkers.
"""

from .oxli import KmerCountTable

__all__ = ['KmerCountTable']
