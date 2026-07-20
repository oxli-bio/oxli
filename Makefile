.PHONY: all install test bench bench-py wheel sdist

PYTHON ?= python

all:
	maturin develop --all-features

# Rust/criterion benchmarks (downloads the E. coli genome on first run).
bench:
	cargo bench

# Python/pytest-codspeed benchmarks over doc/example.fa. Requires `.[test]`.
bench-py:
	$(PYTHON) -m pytest src/python/benchmarks --codspeed

install:
	$(PYTHON) -m pip install -e .

test:
	$(PYTHON) -m pytest

wheel:
	$(PYTHON) -m maturin build -r --all-features

sdist:
	rm -f target/wheels/oxli-*.tar.gz
	$(PYTHON) -m maturin sdist
