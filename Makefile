PYTHON ?= python

all:
	maturin develop

install:
	$(PYTHON) -m pip install -e .

test:
	$(PYTHON) -m pytest

wheel:
	$(PYTHON) -m maturin build -r
