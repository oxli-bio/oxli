// Standard library imports
use std::collections::{HashMap, HashSet};
use std::collections::hash_map::IntoIter;

// External crate imports
use anyhow::{anyhow, Result};
use log::debug;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use sourmash::encodings::HashFunctions;
use sourmash::signature::SeqToHashes;

// Set version variable
const VERSION: &str = env!("CARGO_PKG_VERSION");

#[pyclass]
struct KmerCountTable {
    counts: HashMap<u64, u64>,
    pub ksize: u8,
    version: String,
    consumed: u64,
}

#[pymethods]
impl KmerCountTable {
    #[new]
    #[pyo3(signature = (ksize))]
    pub fn new(ksize: u8) -> Self {
        Self {
            counts: HashMap::new(),
            ksize,
            version: VERSION.to_string(), // Initialize the version field
            consumed: 0,                  // Initialize the total sequence length tracker
        }
    }

    // TODO: Optionally store hash:kmer pair when counting a new kmer
    // Modify KmerCountTable to optionally store map of hash:kmer
    // Modify SeqToHashes to return canonical kmer & hash

    // TODO: Add function to get canonical kmer using hash key

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
        let count = self.counts.entry(hashval).or_insert(0);
        *count += 1;
        *count
    }

    pub fn count(&mut self, kmer: String) -> PyResult<u64> {
        if kmer.len() as u8 != self.ksize {
            Err(PyValueError::new_err(
                "kmer size does not match count table ksize",
            ))
        } else {
            self.consumed += kmer.len() as u64;
            let hashval = self.hash_kmer(kmer)?;
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
            debug!("get: hashval {}, count {}", hashval, count);
            Ok(*count)
        }
    }

    // Get the count for a specific hash value directly
    pub fn get_hash(&self, hashval: u64) -> u64 {
        // Return the count for the hash value, or 0 if it does not exist
        *self.counts.get(&hashval).unwrap_or(&0)
    }

    // Get counts for a list of hash keys and return an list of counts
    pub fn get_hash_array(&self, hash_keys: Vec<u64>) -> Vec<u64> {
        // Map each hash key to its count, defaulting to 0 if the key is not present
        hash_keys.iter().map(|&key| self.get_hash(key)).collect()
    }

    // TODO: Add method "drop"
    // remove kmer from table

    // TODO: Add method "drop_hash"
    // remove hash from table

    // TODO: Add "mincut". Remove counts below a minimum cutoff.

    // TODO: Add "maxcut". Remove counts above an maximum cutoff.

    // TODO: Serialize the KmerCountTable instance to a JSON string.

    // TODO: Compress JSON string with gzip and save to file

    // TODO: Static method to load KmerCountTable from serialized JSON. Yield new object.

    // TODO: Add method "dump"
    // Output tab delimited kmer:count pairs
    // Default sort by count
    // Option sort kmers lexicographically

    // TODO: Add method "dump_hash"
    // Output tab delimited hash:count pairs
    // Default sort by count
    // Option sort on keys

    // TODO: Add method "histo"
    // Output frequency counts

    // Getter for the 'hashes' attribute, returning all hash keys in the table
    #[getter]
    pub fn hashes(&self) -> Vec<u64> {
        // Collect and return all keys from the counts HashMap
        self.counts.keys().cloned().collect()
    }

    // Attribute to access the version of oxli that the table was created with
    #[getter]
    pub fn version(&self) -> &str {
        &self.version
    }

    // Attribute to access the total bases processed with count or consume.
    #[getter]
    pub fn consumed(&self) -> u64 {
        self.consumed
    }

    // Getter for the sum of all counts in the table.
    #[getter]
    pub fn sum_counts(&self) -> u64 {
        self.counts.values().sum()
    }

    // Consume this DNA string. Return number of k-mers consumed.
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

        // Update the total sequence consumed tracker
        self.consumed += seq.len() as u64;

        Ok(n)
    }

    // Helper method to get hash set of k-mers
    fn hash_set(&self) -> HashSet<u64> {
        self.counts.keys().cloned().collect()
    }

    // Set operation methods
    pub fn union(&self, other: &KmerCountTable) -> HashSet<u64> {
        self.hash_set().union(&other.hash_set()).cloned().collect()
    }

    pub fn intersection(&self, other: &KmerCountTable) -> HashSet<u64> {
        self.hash_set()
            .intersection(&other.hash_set())
            .cloned()
            .collect()
    }

    pub fn difference(&self, other: &KmerCountTable) -> HashSet<u64> {
        self.hash_set()
            .difference(&other.hash_set())
            .cloned()
            .collect()
    }

    pub fn symmetric_difference(&self, other: &KmerCountTable) -> HashSet<u64> {
        self.hash_set()
            .symmetric_difference(&other.hash_set())
            .cloned()
            .collect()
    }

    // Python dunder methods for set operations
    fn __or__(&self, other: &KmerCountTable) -> HashSet<u64> {
        self.union(other)
    }

    fn __and__(&self, other: &KmerCountTable) -> HashSet<u64> {
        self.intersection(other)
    }

    fn __sub__(&self, other: &KmerCountTable) -> HashSet<u64> {
        self.difference(other)
    }

    fn __xor__(&self, other: &KmerCountTable) -> HashSet<u64> {
        self.symmetric_difference(other)
    }

    // Python __iter__ method to return an iterator
    pub fn __iter__(slf: PyRef<Self>) -> KmerCountTableIterator {
        KmerCountTableIterator {
            inner: slf.counts.clone().into_iter(), // Clone the HashMap and convert to iterator
        }
    }

    // Python dunder method for __len__
    fn __len__(&self) -> usize {
        self.counts.len()
    }

    // Python dunder method for __getitem__
    fn __getitem__(&self, kmer: String) -> PyResult<u64> {
        self.get(kmer)
    }

    // Python dunder method for __setitem__
    pub fn __setitem__(&mut self, kmer: String, count: u64) -> PyResult<()> {
        // Calculate the hash for the k-mer
        let hashval = self.hash_kmer(kmer)?;
        // Set the count for the k-mer
        self.counts.insert(hashval, count);
        Ok(())
    }
}

// Iterator implementation for KmerCountTable
#[pyclass]
pub struct KmerCountTableIterator {
    inner: IntoIter<u64, u64>, // Now we own the iterator
}

#[pymethods]
impl KmerCountTableIterator {
    pub fn __next__(mut slf: PyRefMut<Self>) -> Option<(u64, u64)> {
        slf.inner.next()
    }
}

// Python module definition
#[pymodule]
fn oxli(m: &Bound<'_, PyModule>) -> PyResult<()> {
    env_logger::init();
    m.add_class::<KmerCountTable>()?;
    Ok(())
}
