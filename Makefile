.PHONY: all install test wheel sdist

PYTHON ?= python

all:
	maturin develop

install:
	$(PYTHON) -m pip install -e .

test:
	$(PYTHON) -m pytest

wheel:
	$(PYTHON) -m maturin build -r

sdist:
	rm -f target/wheels/oxli-*.tar.gz
	$(PYTHON) -m maturin sdist
