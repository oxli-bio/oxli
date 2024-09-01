# A simple example of the API

Import necessary modules:

```python
>>> import screed
>>> import oxli

```

Create a KmerCountTable with a k-mer size of 31:

```python
>>> counts = oxli.KmerCountTable(31)

```

Open a FASTA file and consume k-mers from all the sequences within:

```python
>>> for record in screed.open('example.fa'):
...    counts.consume(record.sequence)
349900

```

Get the count of `CGGAGGAAGCAAGAACAAAATATTTTTTCAT` in the data::

```python
>>> counts.get('CGGAGGAAGCAAGAACAAAATATTTTTTCAT')
1

```
