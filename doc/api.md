# A simple example of the API

Import necessary modules:

```python
>>> import screed # FASTA/FASTQ parsing
>>> import oxli   # this package!

```

## The basics

Create a KmerCountTable with a k-mer size of 31:

```python
>>> counts = oxli.KmerCountTable(ksize=31)

```

Open a FASTA file and consume k-mers from all the sequences within:

```python
>>> for record in screed.open('example.fa'):
...    counts.consume(record.sequence)
349900

```

Here, `consume` reports the total number of k-mers consumed.

Get the count of `CGGAGGAAGCAAGAACAAAATATTTTTTCAT` in the data:

```python
>>> counts.get('CGGAGGAAGCAAGAACAAAATATTTTTTCAT')
1

```

Reverse complement are of course handled:
```python
>>> counts.get('ATGAAAAAATATTTTGTTCTTGCTTCCTCCG')
1

```


## Handling bad k-mers in DNA sequence

You can fail on bad k-mers:

```python
>>> counts.consume('XXXCGGAGGAAGCAAGAACAAAATATTTTTTCATGGG', skip_bad_kmers=False)
Traceback (most recent call last):
...
ValueError: bad k-mer encountered at position 0

```

or allow them (which is default):

```python
>>> counts.consume('XXXCGGAGGAAGCAAGAACAAAATATTTTTTCATGGG', skip_bad_kmers=True)
4

```

If you allow bad k-mers, then all of the valid k-mers will be counted, and all of the bad k-mers will be skipped:
```
>>> counts.get("CGGAGGAAGCAAGAACAAAATATTTTTTCAT")
2

>>> counts.get("AGGAAGCAAGAACAAAATATTTTTTCATGGG")
1

```

@CTB note: right now, all k-mers up to the bad k-mers are
counted... oops. Fix this!
