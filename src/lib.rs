// Standard library imports
use std::collections::hash_map::IntoIter;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::sync::{
    atomic::{AtomicU64, Ordering},
    Mutex,
};
//use std::path::Path;

// External crate imports
use anyhow::{anyhow, Result};
use log::debug;
use niffler::compression::Format;
use niffler::get_writer;
use pyo3::exceptions::{PyIOError, PyKeyError, PyValueError};
use pyo3::prelude::*;
use pyo3::PyResult;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use sourmash::encodings::revcomp;
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
    store_kmers: bool, // Store hash:kmer mapping if true
    hash_to_kmer: Option<HashMap<u64, String>>,
}

#[pymethods]
impl KmerCountTable {
    /// Constructor for KmerCountTable
    #[new]
    #[pyo3(signature = (ksize, store_kmers=false))]
    pub fn new(ksize: u8, store_kmers: bool) -> Self {
        // Optional init HashMap for tracking hash:kmer pairs
        let hash_to_kmer = if store_kmers {
            Some(HashMap::new())
        } else {
            None
        };
        // Init new KmerCountTable
        Self {
            counts: HashMap::new(),
            ksize,
            version: VERSION.to_string(), // Initialize the version field
            consumed: 0,                  // Initialize the total sequence length tracker
            store_kmers,
            hash_to_kmer,
        }
    }

    /// Turn a k-mer into a hashval.
    pub fn hash_kmer(&self, kmer: String) -> Result<u64> {
        if (kmer.len() as u8) != self.ksize {
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

    /// Unhash function to retrieve the canonical kmer for a given hash
    pub fn unhash(&self, hash: u64) -> PyResult<String> {
        if self.store_kmers {
            if let Some(kmer) = self.hash_to_kmer.as_ref().unwrap().get(&hash) {
                Ok(kmer.clone())
            } else {
                // Raise KeyError if hash does not exist
                let msg = format!("Warning: Hash {} not found in table.", hash);
                Err(PyKeyError::new_err(msg))
            }
        } else {
            // Raise an error if store_kmers is false
            Err(PyValueError::new_err("K-mer storage is not enabled."))
        }
    }

    /// Increment the count of a hashval by 1.
    pub fn count_hash(&mut self, hashval: u64) -> u64 {
        let count = self.counts.entry(hashval).or_insert(0);
        *count += 1;
        *count
    }

    /// Return the canonical form of a k-mer: the lexicographically smaller of the k-mer or its reverse complement.
    fn canon(&self, kmer: &str) -> PyResult<String> {
        // Check if the k-mer length matches the table ksize
        if kmer.len() != self.ksize as usize {
            return Err(PyValueError::new_err(
                "kmer size does not match count table ksize",
            ));
        }

        // Convert k-mer to uppercase
        let kmer_upper = kmer.to_uppercase();

        // Ensure k-mer contains only valid DNA characters
        if !kmer_upper.chars().all(|c| "ATCG".contains(c)) {
            return Err(PyValueError::new_err("kmer contains invalid characters"));
        }

        // Compute the reverse complement
        let rev_comp: String = kmer_upper
            .chars()
            .rev()
            .map(|c| match c {
                'A' => 'T',
                'T' => 'A',
                'C' => 'G',
                'G' => 'C',
                _ => c, // This should not happen due to earlier validation
            })
            .collect();

        // Return the lexicographically smaller of kmer or its reverse complement
        if kmer_upper <= rev_comp {
            Ok(kmer_upper)
        } else {
            Ok(rev_comp)
        }
    }

    /// Increment the count of a k-mer by 1.
    pub fn count(&mut self, kmer: String) -> PyResult<u64> {
        if kmer.len() as u8 != self.ksize {
            Err(PyValueError::new_err(
                "kmer size does not match count table ksize",
            ))
        } else {
            let hashval = self.hash_kmer(kmer.clone())?; // Clone the kmer before passing it to hash_kmer
            let count = self.count_hash(hashval); // count with count_hash() function, return tally
            self.consumed += kmer.len() as u64; // Add kmer len to total consumed bases

            if self.store_kmers {
                // Get the canonical k-mer
                let canonical_kmer = self.canon(&kmer)?;
                // Optional: Store hash:kmer pair
                self.hash_to_kmer
                    .as_mut()
                    .unwrap()
                    .insert(hashval, canonical_kmer);
            }

            Ok(count) // Return the current total count for the hash
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

            let count = self.counts.get(&hashval).unwrap_or(&0);
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

    /// Dump (canonical_kmer,count) pairs, optional sorted by count or canonical kmer.
    ///
    /// # Arguments
    /// * `file` - Optional file path to write the output. If not provided, returns a list of tuples.
    /// * `sortkeys` - Optional flag to sort by canonical kmers (default: False).
    /// * `sortcounts` - Sort on counts, secondary sort on canonical kmers. (default: False).
    #[pyo3(signature = (file=None, sortcounts=false, sortkeys=false))]
    pub fn dump_kmers(
        &self,
        file: Option<String>,
        sortcounts: bool,
        sortkeys: bool,
    ) -> PyResult<Vec<(String, u64)>> {
        // Ensure that the hash:kmer mapping is stored
        if !self.store_kmers {
            return Err(PyValueError::new_err(
                "K-mer storage is disabled. No hash:kmer map is available.",
            ));
        }

        // Raise an error if both sortcounts and sortkeys are true
        if sortcounts && sortkeys {
            return Err(PyValueError::new_err(
                "Cannot sort by both counts and kmers at the same time.",
            ));
        }

        // Collect canonical k-mers and their counts, skipping those not found in the counts table
        let mut kmer_count_pairs: Vec<(&String, &u64)> = self
            .hash_to_kmer
            .as_ref()
            .unwrap()
            .par_iter() // Use rayon for parallel iteration
            .filter_map(|(&hash, kmer)| {
                // Use filter_map to only include (kmer, count) pairs where the count exists
                self.counts.get(&hash).map(|count| (kmer, count))
            })
            .collect();

        // Handle sorting based on the flags
        if sortkeys {
            // Sort by canonical kmer lexicographically
            kmer_count_pairs.par_sort_by_key(|&(kmer, _)| kmer.clone());
        } else if sortcounts {
            // Sort by count, secondary sort by kmer
            kmer_count_pairs.par_sort_by(|&(kmer1, count1), &(kmer2, count2)| {
                count1.cmp(count2).then_with(|| kmer1.cmp(kmer2))
            });
        }
        // If both sortcounts and sortkeys are false, no sorting is done.

        // If a file is provided, write to the file
        if let Some(filepath) = file {
            let f = File::create(filepath)?;
            let mut writer = BufWriter::new(f);

            // Write each kmer:count pair to the file
            for (kmer, count) in kmer_count_pairs {
                writeln!(writer, "{}\t{}", kmer, count)?;
            }

            writer.flush()?; // Ensure all data is written to the file
            Ok(vec![]) // Return an empty vector when writing to a file
        } else {
            // Convert the vector of references to owned values
            let result: Vec<(String, u64)> = kmer_count_pairs
                .into_par_iter() // Use rayon for parallel conversion
                .map(|(kmer, &count)| (kmer.clone(), count))
                .collect();

            // Return the vector of (kmer, count) tuples
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

    // Consume this DNA string. Return total number of k-mers consumed.
    // If "skip_bad_kmers = true" then ignore kmers with non-DNA characters
    // else if "false" consume kmers until a bad kmer in encountered, then
    // exit with error.
    #[pyo3(signature = (seq, skip_bad_kmers=true))]
    pub fn consume(&mut self, seq: &str, skip_bad_kmers: bool) -> PyResult<u64> {
        // Incoming seq len
        let new_len = seq.len();
        // Init tally for consumed kmers
        let mut n = 0;
        // If store_kmers is true, then count & log hash:kmer pairs
        if self.store_kmers {
            let hash_to_kmer = self.hash_to_kmer.as_mut().unwrap();

            // Create an iterator for (canonical_kmer, hash) pairs
            let iter = KmersAndHashesIter::new(seq, self.ksize as usize, skip_bad_kmers);

            // Iterate over the k-mers and their hashes
            for result in iter {
                match result {
                    Ok((kmer, hash)) => {
                        if hash != 0 {
                            // Insert hash:kmer pair into the hashmap
                            hash_to_kmer.insert(hash, kmer.clone());
                            // Increment the count for the hash
                            *self.counts.entry(hash).or_insert(0) += 1;
                            // Tally kmers added
                            n += 1;
                        }
                    }
                    Err(e) => return Err(e),
                }
            }
        } else {
            // Else, hash and count kmers as usual
            let hashes = SeqToHashes::new(
                seq.as_bytes(),
                self.ksize.into(),
                skip_bad_kmers,
                false,
                HashFunctions::Murmur64Dna,
                42,
            );

            for hash_value in hashes {
                // eprintln!("hash_value: {:?}", hash_value);
                match hash_value {
                    Ok(0) => continue,
                    Ok(x) => {
                        self.count_hash(x);
                    }
                    Err(_) => {
                        let msg = format!("bad k-mer encountered at position {}", n);
                        return Err(PyValueError::new_err(msg));
                    }
                }

                n += 1;
            }
        }

        // Update the total sequence consumed tracker
        self.consumed += new_len as u64;

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

    #[pyo3(signature = (seq, skip_bad_kmers=true))]
    pub fn kmers_and_hashes(
        &self,
        seq: &str,
        skip_bad_kmers: bool,
    ) -> PyResult<Vec<(String, u64)>> {
        let mut v: Vec<(String, u64)> = vec![];

        // Create the iterator
        let iter = KmersAndHashesIter::new(seq, self.ksize as usize, skip_bad_kmers);

        // Collect the k-mers and their hashes
        for result in iter {
            match result {
                Ok((kmer, hash)) => v.push((kmer, hash)),
                Err(e) => return Err(e),
            }
        }

        Ok(v)
    }

    /// Calculates the Jaccard Similarity Coefficient between two KmerCountTable objects.
    /// # Returns
    /// The Jaccard Similarity Coefficient between the two tables as a float value between 0 and 1.
    pub fn jaccard(&self, other: &KmerCountTable) -> f64 {
        // Get the intersection of the two k-mer sets.
        let intersection_size = self.intersection(other).len();

        // Get the union of the two k-mer sets.
        let union_size = self.union(other).len();

        // Handle the case where the union is empty (both sets are empty).
        if union_size == 0 {
            return 1.0; // By convention, two empty sets are considered identical.
        }

        // Calculate and return the Jaccard similarity as a ratio of intersection to union.
        intersection_size as f64 / union_size as f64
    }

    /// Cosine similarity between two `KmerCountTable` objects.
    /// # Returns
    /// The cosine similarity between the two tables as a float value between 0 and 1.
    pub fn cosine(&self, other: &KmerCountTable) -> f64 {
        // Early return if either table is empty.
        if self.counts.is_empty() || other.counts.is_empty() {
            return 0.0;
        }

        // Calculate the dot product in parallel.
        let dot_product: u64 = self
            .counts
            .par_iter()
            .filter_map(|(&hash, &count1)| {
                // Only include in the dot product if both tables have the k-mer.
                other.counts.get(&hash).map(|&count2| count1 * count2)
            })
            .sum();

        // Calculate magnitudes in parallel for both tables.
        let magnitude_self: f64 = self
            .counts
            .par_iter()
            .map(|(_, v)| (*v as f64).powi(2)) // Access the value, square it
            .sum::<f64>()
            .sqrt();

        let magnitude_other: f64 = other
            .counts
            .par_iter()
            .map(|(_, v)| (*v as f64).powi(2)) // Access the value, square it
            .sum::<f64>()
            .sqrt();

        // If either magnitude is zero (no k-mers), return 0 to avoid division by zero.
        if magnitude_self == 0.0 || magnitude_other == 0.0 {
            return 0.0;
        }

        // Calculate and return cosine similarity.
        dot_product as f64 / (magnitude_self * magnitude_other)
    }

    /// Add counts from another KmerCountTable to this one.
    ///
    /// # Arguments
    ///
    /// * `other` - The KmerCountTable to add from
    ///
    /// # Returns
    ///
    /// Returns a PyResult with a tuple containing:
    /// * The number of k-mer counts added
    /// * The number of new keys added
    #[pyo3(signature = (other))]
    pub fn add(&mut self, other: &KmerCountTable) -> PyResult<(u64, u64)> {
        if self.ksize != other.ksize {
            return Err(PyValueError::new_err(
                "KmerCountTables must have the same ksize",
            ));
        }

        let total_counts_added = AtomicU64::new(0);
        let new_keys_added = AtomicU64::new(0);
        let counts_mutex = Mutex::new(&mut self.counts);

        // Use thread-local storage to collect updates
        let updates: Vec<_> = other
            .counts
            .par_iter()
            .map(|(&hash, &count)| (hash, count))
            .collect();

        // Apply updates in parallel
        updates.par_iter().for_each(|(hash, count)| {
            let mut counts_lock = counts_mutex.lock().unwrap();
            let current_count = counts_lock.entry(*hash).or_insert(0);
            if *current_count == 0 {
                new_keys_added.fetch_add(1, Ordering::Relaxed);
            }
            *current_count += count;
            total_counts_added.fetch_add(*count, Ordering::Relaxed);
        });

        self.consumed += other.consumed;

        if self.store_kmers {
            if other.store_kmers {
                let hash_to_kmer_mutex = Mutex::new(self.hash_to_kmer.as_mut().unwrap());

                other
                    .hash_to_kmer
                    .as_ref()
                    .unwrap()
                    .par_iter()
                    .for_each(|(&hash, kmer)| {
                        let mut hash_to_kmer_lock = hash_to_kmer_mutex.lock().unwrap();
                        hash_to_kmer_lock
                            .entry(hash)
                            .or_insert_with(|| kmer.clone());
                    });
            } else {
                eprintln!("Warning: Incoming table does not store k-mers, but target table does. K-mer information for new hashes will be missing.");
            }
        }

        let total_added = total_counts_added.load(Ordering::Relaxed);
        let new_keys = new_keys_added.load(Ordering::Relaxed);

        println!("Added {} k-mer counts to the table", total_added);
        println!("Added {} new keys to the table", new_keys);

        Ok((total_added, new_keys))
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

pub struct KmersAndHashesIter {
    seq: String,          // The sequence to iterate over
    seq_rc: String,       // reverse complement sequence
    ksize: usize,         // K-mer size
    pos: usize,           // Current position in the sequence
    end: usize,           // The end position for k-mer extraction
    hasher: SeqToHashes,  // Iterator for generating hashes
    skip_bad_kmers: bool, // Flag to skip bad k-mers
}

impl KmersAndHashesIter {
    pub fn new(seq: &str, ksize: usize, skip_bad_kmers: bool) -> Self {
        let seq = seq.to_ascii_uppercase(); // Ensure uppercase for uniformity
        let seqb = seq.as_bytes().to_vec(); // Convert to bytes for hashing
        let seqb_rc = revcomp(&seqb);
        let seq_rc = std::str::from_utf8(&seqb_rc)
            .expect("invalid utf-8 sequence for rev comp")
            .to_string();

        let end = seq.len() - ksize + 1; // Calculate the endpoint for k-mer extraction
        let hasher = SeqToHashes::new(
            &seqb,
            ksize,
            true,  // Set force to true, bad kmers will emit hash=0 instead of killing process
            false, // Other flags, e.g., reverse complement
            HashFunctions::Murmur64Dna,
            42, // Seed for hashing
        );

        Self {
            seq,
            seq_rc,
            ksize,
            pos: 0, // Start at the beginning of the sequence
            end,
            hasher,
            skip_bad_kmers,
        }
    }
}

impl Iterator for KmersAndHashesIter {
    type Item = PyResult<(String, u64)>;

    fn next(&mut self) -> Option<Self::Item> {
        // Check if we've reached the end of the sequence
        if self.pos >= self.end {
            return None;
        }

        let start = self.pos;
        let ksize = self.ksize;
        let rpos = self.end - start - 1;

        // Extract the current k-mer and its reverse complement
        let substr = &self.seq[start..start + ksize];
        let substr_rc = &self.seq_rc[rpos..rpos + ksize];

        // Get the next hash value from the hasher
        let hashval = self.hasher.next().expect("should not run out of hashes");

        // Increment position for the next k-mer
        self.pos += 1;

        // Handle hash value logic
        if let Ok(hashval) = hashval {
            // Good kmer, all is well, store canonical k-mer and hashval;
            if hashval > 0 {
                // Select the canonical k-mer (lexicographically smaller between forward and reverse complement)
                let canonical_kmer = if substr < substr_rc {
                    substr
                } else {
                    substr_rc
                };
                // If valid hash, return (canonical_kmer,hashval) tuple
                Some(Ok((canonical_kmer.to_string(), hashval)))
            } else {
                // If the hash is 0, handle based on `skip_bad_kmers`
                // Prepare msg identifying bad kmer
                let msg = format!("bad k-mer at position {}: {}", start + 1, substr);
                if self.skip_bad_kmers {
                    // Print a message and skip adding the bad k-mer to the result
                    eprintln!("{}", msg);
                    self.next() // Recursively call `next()` to skip this k-mer
                } else {
                    // If skip_bad_kmer is false, return an empty string and 0, but still print a message
                    eprintln!("{}", msg);
                    Some(Ok(("".to_string(), 0)))
                }
            }
        } else {
            // If error raised by SeqToHashes
            let msg = format!("bad k-mer at position {}: {}", start + 1, substr);
            Some(Err(PyValueError::new_err(msg)))
        }
    }
}

// Python module definition
#[pymodule]
fn oxli(m: &Bound<'_, PyModule>) -> PyResult<()> {
    env_logger::init();
    m.add_class::<KmerCountTable>()?;
    Ok(())
}
