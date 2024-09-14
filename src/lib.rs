// Standard library imports
use std::collections::{HashMap, HashSet};

// External crate imports
use anyhow::{anyhow, Result};
use log::debug;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
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
            let hashval = self.hash_kmer(kmer).unwrap();
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

    /// Drop a k-mer from the count table by its string representation
    pub fn drop(&mut self, kmer: String) -> PyResult<()> {
        // Ensure that the k-mer length matches the table's ksize
        if kmer.len() as u8 != self.ksize {
            // Return an error if the lengths do not match
            Err(PyValueError::new_err(
                "kmer size does not match count table ksize",
            ))
        } else {
            // Compute the hash of the k-mer using the same method used for counting
            let hashval = self.hash_kmer(kmer).unwrap();

            // Attempt to remove the k-mer's hash from the counts HashMap
            if self.counts.remove(&hashval).is_some() {
                // If the k-mer was successfully removed, return Ok
                debug!("K-mer with hashval {} removed from table", hashval);
                Ok(())
            } else {
                // If the k-mer was not found, return Ok without an error
                debug!("K-mer with hashval {} not found in table", hashval);
                Ok(())
            }
        }
    }

    /// Drop a k-mer from the count table by its hash value
    pub fn drop_hash(&mut self, hashval: u64) -> PyResult<()> {
        // Attempt to remove the hash value from the counts HashMap
        if self.counts.remove(&hashval).is_some() {
            // If the hash value was successfully removed, log and return Ok
            debug!("Hash value {} removed from table", hashval);
            Ok(())
        } else {
            // If the hash value was not found, log and return Ok without error
            debug!("Hash value {} not found in table", hashval);
            Ok(())
        }
    }

    /// Remove all k-mers with counts less than a given threshold
    pub fn mincut(&mut self, min_count: u64) -> PyResult<u64> {
        // Create a vector to store the keys (hashes) to be removed
        let mut to_remove = Vec::new();

        // Iterate over the HashMap and identify keys with counts less than the threshold
        for (&hash, &count) in self.counts.iter() {
            if count < min_count {
                to_remove.push(hash);
            }
        }

        // Remove the identified keys from the counts HashMap
        for &hash in &to_remove {
            self.counts.remove(&hash);
        }

        // Return the number of k-mers removed
        Ok(to_remove.len() as u64)
    }

    /// Remove all k-mers with counts greater than a given threshold
    pub fn maxcut(&mut self, max_count: u64) -> PyResult<u64> {
        // Create a vector to store the keys (hashes) to be removed
        let mut to_remove = Vec::new();

        // Iterate over the HashMap and identify keys with counts greater than the threshold
        for (&hash, &count) in self.counts.iter() {
            if count > max_count {
                to_remove.push(hash);
            }
        }

        // Remove the identified keys from the counts HashMap
        for &hash in &to_remove {
            self.counts.remove(&hash);
        }

        // Return the number of k-mers removed
        Ok(to_remove.len() as u64)
    }

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

    // TODO: Getter for the version attribute
    // Store oxli version when instance is created

    // TODO: Getter for the consumed seq len attribute
    // Update tracker when DNA is processed with count() or consume()

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

    // Python dunder method for __iter__

    // Python dunder method for __next__

    // Python dunder method for __len__

    // Python dunder method for __getitem__

    // Python dunder method for __setitem__
}

// Python module definition
#[pymodule]
fn oxli(m: &Bound<'_, PyModule>) -> PyResult<()> {
    env_logger::init();
    m.add_class::<KmerCountTable>()?;
    Ok(())
}
