// Standard library imports
use std::collections::hash_map::IntoIter;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
//use std::path::Path;

// External crate imports
use anyhow::{anyhow, Result};
use log::debug;
use niffler::compression::Format;
use niffler::get_writer;
use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;
use pyo3::PyResult;
use serde::{Deserialize, Serialize};
use sourmash::encodings::HashFunctions;
use sourmash::signature::SeqToHashes;

// Set version variable
const VERSION: &str = env!("CARGO_PKG_VERSION");

#[pyclass]
#[derive(Serialize, Deserialize, Debug)]
/// Basic KmerCountTable struct, mapping hashes to counts.
struct KmerCountTable {
    counts: HashMap<u64, u64>,
    pub ksize: u8,
    version: String,
    consumed: u64,
}

#[pymethods]
/// Methods on KmerCountTable.
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

    /// Turn a k-mer into a hashval.
    fn hash_kmer(&self, kmer: String) -> Result<u64> {
        if kmer.len() as u8 != self.ksize {
            Err(anyhow!("wrong ksize"))
        } else {
            let mut hashes = SeqToHashes::new(
                kmer.as_bytes(),
                self.ksize.into(),
                false,
                false,
                HashFunctions::Murmur64Dna,
                42,
            );

            let hashval = hashes.next().expect("error hashing this k-mer");
            Ok(hashval?)
        }
    }

    /// Increment the count of a hashval by 1.
    pub fn count_hash(&mut self, hashval: u64) -> u64 {
        let count = self.counts.entry(hashval).or_insert(0);
        *count += 1;
        *count
    }

    /// Increment the count of a k-mer by 1.
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

    /// Retrieve the count of a k-mer.
    pub fn get(&self, kmer: String) -> PyResult<u64> {
        if kmer.len() as u8 != self.ksize {
            Err(PyValueError::new_err(
                "kmer size does not match count table ksize",
            ))
        } else {
            let hashval = self.hash_kmer(kmer).expect("error hashing this k-mer");

            let count = match self.counts.get(&hashval) {
                Some(count) => count,
                None => &0,
            };
            debug!("get: hashval {}, count {}", hashval, count);
            Ok(*count)
        }
    }

    /// Get the count for a specific hash value directly
    pub fn get_hash(&self, hashval: u64) -> u64 {
        // Return the count for the hash value, or 0 if it does not exist
        *self.counts.get(&hashval).unwrap_or(&0)
    }

    /// Get counts for a list of hashvals and return a list of counts
    pub fn get_hash_array(&self, hash_keys: Vec<u64>) -> Vec<u64> {
        // Map each hash key to its count, defaulting to 0 if the key is not present
        hash_keys.iter().map(|&key| self.get_hash(key)).collect()
    }

    /// Drop a k-mer from the count table by its string representation
    pub fn drop(&mut self, kmer: String) -> PyResult<()> {
        // Compute the hash of the k-mer using the same method used for counting
        let hashval = self.hash_kmer(kmer)?;
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

    /// Serialize the KmerCountTable as a JSON string
    pub fn serialize_json(&self) -> Result<String> {
        serde_json::to_string(&self).map_err(|e| anyhow::anyhow!("Serialization error: {}", e))
    }

    /// Save the KmerCountTable to a compressed file using Niffler.
    pub fn save(&self, filepath: &str) -> PyResult<()> {
        // Open the file for writing
        let file = File::create(filepath).map_err(|e| PyIOError::new_err(e.to_string()))?;

        // Create a Gzipped writer with niffler, using the default compression level
        let writer = BufWriter::new(file);
        let mut writer = get_writer(Box::new(writer), Format::Gzip, niffler::level::Level::One)
            .map_err(|e| PyIOError::new_err(e.to_string()))?;

        // Serialize the KmerCountTable to JSON
        let json_data = self.serialize_json()?;

        // Write the serialized JSON to the compressed file
        writer
            .write_all(json_data.as_bytes())
            .map_err(|e| PyIOError::new_err(e.to_string()))?;

        Ok(())
    }

    #[staticmethod]
    /// Load a KmerCountTable from a compressed file using Niffler.
    pub fn load(filepath: &str) -> Result<KmerCountTable> {
        // Open the file for reading
        let file = File::open(filepath)?;

        // Use Niffler to get a reader that detects the compression format
        let reader = BufReader::new(file);
        let (mut reader, _format) = niffler::get_reader(Box::new(reader))?;

        // Read the decompressed data into a string
        let mut decompressed_data = String::new();
        reader.read_to_string(&mut decompressed_data)?;

        // Deserialize the JSON string to a KmerCountTable
        let loaded_table: KmerCountTable = serde_json::from_str(&decompressed_data)
            .map_err(|e| anyhow::anyhow!("Deserialization error: {}", e))?;

        // Check version compatibility and issue a warning if necessary
        if loaded_table.version != VERSION {
            eprintln!(
                "Version mismatch: loaded version is {}, but current version is {}",
                loaded_table.version, VERSION
            );
        }

        Ok(loaded_table)
    }

    /// Dump (hash,count) pairs, optional sorted by count or hash key.
    ///
    /// # Arguments
    /// * `file` - Optional file path to write the output. If not provided, returns a list of tuples.
    /// * `sortkeys` - Optional flag to sort by hash keys (default: False).
    /// * `sortcounts` - Sort on counts, secondary sort on keys. (default: False).
    #[pyo3(signature = (file=None, sortcounts=false, sortkeys=false))]
    pub fn dump(
        &self,
        file: Option<String>,
        sortcounts: bool,
        sortkeys: bool,
    ) -> PyResult<Vec<(u64, u64)>> {
        // Raise an error if both sortcounts and sortkeys are true
        if sortcounts && sortkeys {
            return Err(PyValueError::new_err(
                "Cannot sort by both counts and keys at the same time.",
            ));
        }

        // Collect hashes and counts
        let mut hash_count_pairs: Vec<(&u64, &u64)> = self.counts.iter().collect();

        // Handle sorting based on the flags
        if sortkeys {
            // Sort by hash keys if `sortkeys` is set to true
            hash_count_pairs.sort_by_key(|&(hash, _)| *hash);
        } else if sortcounts {
            // Sort by count, secondary sort by hash if `sortcounts` is true
            hash_count_pairs.sort_by(|&(hash1, count1), &(hash2, count2)| {
                count1.cmp(count2).then_with(|| hash1.cmp(hash2))
            });
        }
        // If both sortcounts and sortkeys are false, no sorting is done.

        // If a file is provided, write to the file
        if let Some(filepath) = file {
            let f = File::create(filepath)?;
            let mut writer = BufWriter::new(f);

            // Write each hash:count pair to the file
            for (hash, count) in hash_count_pairs {
                writeln!(writer, "{}\t{}", hash, count)?;
            }

            writer.flush()?; // Flush the buffer
            Ok(vec![]) // Return empty vector to Python
        } else {
            // Convert the vector of references to owned values
            let result: Vec<(u64, u64)> = hash_count_pairs
                .into_iter()
                .map(|(&hash, &count)| (hash, count))
                .collect();

            // Return the vector of (hash, count) tuples
            Ok(result)
        }
    }

    /// Calculates the frequency histogram for k-mer counts
    /// Returns a vector of tuples (frequency, count), where 'frequency' is
    /// the observed number of times a k-mer count occurred and 'count' is
    /// how many different k-mers have that frequency.
    /// If `zero` is True, include all frequencies from 0 to max observed count,
    /// even if no k-mers were observed for those frequencies.
    #[pyo3(signature = (zero=true))]
    pub fn histo(&self, zero: bool) -> Vec<(u64, u64)> {
        let mut freq_count: HashMap<u64, u64> = HashMap::new();

        // Step 1: Count the frequencies of observed k-mer counts
        for &count in self.counts.values() {
            *freq_count.entry(count).or_insert(0) += 1;
        }

        let mut histo_vec: Vec<(u64, u64)>;

        if zero {
            // Step 2 (optional): Include all frequencies from 0 to max_count
            let max_count = self.max();
            histo_vec = (0..=max_count)
                .map(|freq| (freq, *freq_count.get(&freq).unwrap_or(&0)))
                .collect();
        } else {
            // Step 2: Only include observed frequencies
            histo_vec = freq_count.into_iter().collect();
            histo_vec.sort_by_key(|&(frequency, _)| frequency);
        }

        histo_vec
    }

    /// Finds and returns the minimum count in the counts HashMap.
    /// Returns 0 if the HashMap is empty.
    #[getter]
    pub fn min(&self) -> u64 {
        // Check if the HashMap is empty, return 0 if true
        if self.counts.is_empty() {
            return 0;
        }

        // Iterate over the counts and find the minimum value
        *self.counts.values().min().unwrap_or(&0)
    }

    /// Finds and returns the maximum count in the counts HashMap.
    /// Returns 0 if the HashMap is empty.
    #[getter]
    pub fn max(&self) -> u64 {
        // Check if the HashMap is empty, return 0 if true
        if self.counts.is_empty() {
            return 0;
        }

        // Iterate over the counts and find the maximum value
        *self.counts.values().max().unwrap_or(&0)
    }

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

#[pyclass]
/// Iterator implementation for KmerCountTable
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
