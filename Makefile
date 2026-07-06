.PHONY: all install test wheel sdist

PYTHON ?= python

all:
	maturin develop --all-features

bench:
	cargo bench

install:
	$(PYTHON) -m pip install -e .

test:
	$(PYTHON) -m pytest

wheel:
	$(PYTHON) -m maturin build -r --all-features

sdist:
	rm -f target/wheels/oxli-*.tar.gz
	$(PYTHON) -m maturin sdist
