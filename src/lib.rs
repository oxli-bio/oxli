use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
// use rayon::prelude::*;

use anyhow::{anyhow, Result};
use std::collections::HashMap;

// use sourmash::sketch::nodegraph::Nodegraph;
use sourmash::_hash_murmur;
use sourmash::encodings::HashFunctions;
use sourmash::signature::SeqToHashes;

#[pyclass]
struct KmerCountTable {
    counts: HashMap<u64, u64>,
    pub ksize: u8,
}

#[pymethods]
impl KmerCountTable {
    #[new]
    #[pyo3(signature = (ksize))]
    pub fn new(ksize: u8) -> Self {
        Self {
            counts: HashMap::new(),
            ksize,
        }
    }

    fn hash_kmer(&self, kmer: String) -> Result<u64> {
        if kmer.len() as u8 != self.ksize {
            Err(anyhow!("wrong ksize"))
        } else {
            // mut?
            let mut hashes = SeqToHashes::new(
                kmer.as_bytes(),
                self.ksize.into(),
                false,
                false,
                HashFunctions::Murmur64Dna,
                42,
            );

            let hashval = hashes.next().unwrap();
            Ok(hashval?)
        }
    }

    pub fn count_hash(&mut self, hashval: u64) -> u64 {
        let mut count: u64 = 1;
        if self.counts.contains_key(&hashval) {
            count = *self.counts.get(&hashval).unwrap();
            count = count + 1;
        }
        self.counts.insert(hashval, count);

        count
    }

    pub fn count(&mut self, kmer: String) -> PyResult<u64> {
        if kmer.len() as u8 != self.ksize {
            Err(PyValueError::new_err(
                "kmer size does not match count table ksize",
            ))
        } else {
            let hashval = _hash_murmur(kmer.as_bytes(), 42);
            let count = self.count_hash(hashval);
            Ok(count)
        }
    }

    pub fn get(&self, kmer: String) -> PyResult<u64> {
        if kmer.len() as u8 != self.ksize {
            Err(PyValueError::new_err(
                "kmer size does not match count table ksize",
            ))
        } else {
            let hashval = self.hash_kmer(kmer).unwrap();

            let count = match self.counts.get(&hashval) {
                Some(count) => count,
                None => &0,
            };
            Ok(*count)
        }
    }

    // Consume this DNA strnig. Return number of k-mers consumed.
    #[pyo3(signature = (seq, allow_bad_kmers=true))]
    pub fn consume(&mut self, seq: String, allow_bad_kmers: bool) -> PyResult<u64> {
        let hashes = SeqToHashes::new(
            seq.as_bytes(),
            self.ksize.into(),
            allow_bad_kmers,
            false,
            HashFunctions::Murmur64Dna,
            42,
        );

        let mut n = 0;
        for hash_value in hashes {
            // eprintln!("hash_value: {:?}", hash_value);
            match hash_value {
                Ok(0) => continue,
                Ok(x) => {
                    self.count_hash(x);
                    ()
                }
                Err(_) => {
                    let msg = format!("bad k-mer encountered at position {}", n);
                    return Err(PyValueError::new_err(msg));
                }
            }

            n += 1;
        }

        Ok(n)
    }
}

#[pymodule]
fn oxli(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<KmerCountTable>()?;
    Ok(())
}
