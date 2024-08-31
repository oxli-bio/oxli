use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
// use rayon::prelude::*;

use anyhow::{Result, Error, anyhow};
use std::collections::HashMap;

// use sourmash::sketch::nodegraph::Nodegraph;
use sourmash::_hash_murmur;
use sourmash::signature::SeqToHashes;
use sourmash::encodings::HashFunctions;


#[pyclass]
struct KmerCountTable {
    counts: HashMap<u64, u64>,
    pub ksize: u8,
}

#[pymethods]
impl KmerCountTable {
    #[new]
    pub fn new(ksize: u8) -> Self {
        Self { counts: HashMap::new(), ksize }
    }

    fn hash_kmer(&self, kmer: String) -> Result<u64> {
        if kmer.len() as u8 != self.ksize {
            Err(anyhow!("wrong ksize"))
        } else {
            // mut?
            let mut hashes = SeqToHashes::new(kmer.as_bytes(),
                                          self.ksize.into(),
                                          false,
                                          false,
                                          HashFunctions::Murmur64Dna,
                                          42);

            let mut hashval = hashes.next().unwrap();
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
            Err(PyValueError::new_err("kmer size does not match count table ksize"))
        } else {
            let hashval = _hash_murmur(kmer.as_bytes(), 42);
            let count = self.count_hash(hashval);
            Ok(count)
        }
    }

    pub fn get(&self, kmer: String) -> PyResult<u64> {
        if kmer.len() as u8 != self.ksize {
            Err(PyValueError::new_err("kmer size does not match count table ksize"))
        } else {
            let hashval = self.hash_kmer(kmer).unwrap();

            let count = match self.counts.get(&hashval) {
                Some(count) => count,
                None => &0
            };
            Ok(*count)
        }
    }

    // Consume this DNA strnig. Return number of k-mers consumed.
    pub fn consume(&mut self, seq: String) -> PyResult<u64> {
        let hashes = SeqToHashes::new(seq.as_bytes(),
                                      self.ksize.into(),
                                      false,
                                      false,
                                      HashFunctions::Murmur64Dna,
                                      42);

        let mut n = 0;
        for hash_value in hashes {
            match hash_value {
                Ok(0) => continue,
                Ok(x) => { self.count_hash(x); () }
                Err(err) => (),
            }
            n += 1;
        }

        Ok(n)
    }
}

#[pymodule]
fn oxli(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<KmerCountTable>()?;
    Ok(())
}
