// Standard library imports
use std::cmp::max;
use std::collections::hash_map::IntoIter;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Cursor, Read, Write};
//use std::path::Path;

// External crate imports
use anyhow::{anyhow, Result};
use log::debug;
use niffler::compression::Format;
use niffler::get_writer;
use nohash_hasher::BuildNoHashHasher;
use pyo3::exceptions::{PyIOError, PyKeyError, PyValueError};
use pyo3::prelude::*;
use pyo3::PyResult;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use sourmash::encodings::revcomp;
use sourmash::encodings::HashFunctions;
use sourmash::signature::SeqToHashes;
// use sourmash::_hash_murmur;
// use sourmash::sketch::nodegraph::Nodegraph;

extern crate needletail;
use needletail::{parse_fastx_file, Sequence};

// Set version variable
const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Magic prefix identifying oxli's binary (bincode) save format. The trailing
/// byte is the on-disk format version. It is written *before* the gzip stream so
/// `load` can distinguish a new binary file from a legacy gzip-JSON one (whose
/// first bytes are the gzip magic `1f 8b`) or an uncompressed JSON file (`{`).
const SAVE_MAGIC: &[u8] = b"OXLIBIN\x01";

/// A `HashMap` keyed by 64-bit k-mer hashes using an identity hasher.
///
/// The keys are MurmurHash64 values (uniformly distributed), so the default
/// SipHash re-hash is wasted work; `BuildNoHashHasher` uses the key directly.
/// `serde` serializes this exactly like a default-hasher `HashMap` (entries
/// only), so the on-disk format is unchanged and old tables still load.
type IntMap<V> = HashMap<u64, V, BuildNoHashHasher<u64>>;

/// Error raised by the pure-Rust counting engine.
///
/// The engine runs inside `Python::detach`, which forbids holding any
/// GIL-bound value (a `PyErr` transitively references Python objects). Errors
/// are therefore reported with this plain Rust type and converted to the
/// appropriate `PyErr` by the caller *after* the GIL is re-acquired, preserving
/// the exact exception messages the Python API has always produced.
enum SeqError {
    /// A non-DNA k-mer was hit at this 0-based k-mer index (`skip_bad_kmers=false`).
    BadKmer(u64),
    /// The input bytes were not valid UTF-8 (only possible on the store-k-mers path).
    NonUtf8,
}

impl SeqError {
    /// Convert to the Python exception the API historically raised.
    fn into_pyerr(self) -> PyErr {
        match self {
            SeqError::BadKmer(pos) => {
                PyValueError::new_err(format!("bad k-mer encountered at position {}", pos))
            }
            SeqError::NonUtf8 => {
                PyValueError::new_err("sequence contains invalid (non-UTF-8) bytes")
            }
        }
    }
}

/// A canonical k-mer stored 2 bits per base (`A=00, C=01, G=10, T=11`).
///
/// Replaces the previous `String` value in the `store_kmers` map: DNA needs only
/// two bits per base, so this uses ~4x less heap than the UTF-8 string and skips
/// per-window UTF-8 handling. `nbases` records the k-mer length so a partial
/// final byte can be decoded unambiguously.
///
/// It (de)serializes transparently *as its decoded string*, so the on-disk JSON
/// (`{hash: "ACGT..."}`) is byte-for-byte identical to older oxli versions and
/// tables saved before this change still load.
#[derive(Debug, Clone, PartialEq, Eq)]
struct PackedKmer {
    nbases: u8,
    bits: Box<[u8]>,
}

impl PackedKmer {
    /// Pack an ASCII, upper-case `ACGT` k-mer. Non-ACGT bytes encode as `A`;
    /// callers only ever pass validated canonical k-mers, so this is not hit in
    /// practice.
    fn encode(kmer: &[u8]) -> Self {
        let nbases = kmer.len();
        let mut bits = vec![0u8; nbases.div_ceil(4)].into_boxed_slice();
        for (i, &base) in kmer.iter().enumerate() {
            let code: u8 = match base {
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => 0, // A (and, defensively, anything unexpected)
            };
            bits[i / 4] |= code << ((i % 4) * 2);
        }
        PackedKmer {
            nbases: nbases as u8,
            bits,
        }
    }

    /// Decode back to the upper-case `ACGT` string that was packed.
    fn decode(&self) -> String {
        const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
        let mut out = Vec::with_capacity(self.nbases as usize);
        for i in 0..self.nbases as usize {
            let code = (self.bits[i / 4] >> ((i % 4) * 2)) & 0b11;
            out.push(BASES[code as usize]);
        }
        // SAFETY: every pushed byte comes from `BASES`, i.e. one of A/C/G/T,
        // so `out` is always valid UTF-8. Skips the redundant validation scan
        // that `String::from_utf8` would run on this hot decode path.
        unsafe { String::from_utf8_unchecked(out) }
    }
}

impl Serialize for PackedKmer {
    fn serialize<S: serde::Serializer>(
        &self,
        serializer: S,
    ) -> std::result::Result<S::Ok, S::Error> {
        serializer.serialize_str(&self.decode())
    }
}

impl<'de> Deserialize<'de> for PackedKmer {
    fn deserialize<D: serde::Deserializer<'de>>(
        deserializer: D,
    ) -> std::result::Result<Self, D::Error> {
        let s = String::deserialize(deserializer)?;
        Ok(PackedKmer::encode(s.as_bytes()))
    }
}

#[pyclass]
#[derive(Serialize, Deserialize, Debug)]
/// Basic KmerCountTable struct, mapping hashes to counts.
pub struct KmerCountTable {
    counts: IntMap<u64>,
    /// K-mer size the table counts; exposed to Python as a read-only attribute.
    #[pyo3(get)]
    pub ksize: u8,
    version: String,
    consumed: u64,
    store_kmers: bool, // Store hash:kmer mapping if true
    hash_to_kmer: Option<IntMap<PackedKmer>>,
}

#[pymethods]
impl KmerCountTable {
    /// Constructor for KmerCountTable
    #[new]
    #[pyo3(signature = (ksize, store_kmers=false))]
    pub fn new(ksize: u8, store_kmers: bool) -> Self {
        // Optional init HashMap for tracking hash:kmer pairs
        let hash_to_kmer = if store_kmers {
            Some(IntMap::default())
        } else {
            None
        };
        // Init new KmerCountTable
        Self {
            counts: IntMap::default(),
            ksize,
            version: VERSION.to_string(), // Initialize the version field
            consumed: 0,                  // Initialize the total sequence length tracker
            store_kmers,
            hash_to_kmer,
        }
    }

    /// Turn a k-mer into a hashval.
    pub fn hash_kmer(&self, kmer: &str) -> Result<u64> {
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
            )?;

            let hashval = hashes.next().expect("error hashing this k-mer");
            Ok(hashval?)
        }
    }

    /// Unhash function to retrieve the canonical kmer for a given hash
    pub fn unhash(&self, hash: u64) -> PyResult<String> {
        if self.store_kmers {
            if let Some(kmer) = self.hash_to_kmer.as_ref().unwrap().get(&hash) {
                Ok(kmer.decode())
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
    pub fn count(&mut self, kmer: &str) -> PyResult<u64> {
        if kmer.len() as u8 != self.ksize {
            Err(PyValueError::new_err(
                "kmer size does not match count table ksize",
            ))
        } else {
            let hashval = self.hash_kmer(kmer)?;
            let count = self.count_hash(hashval); // count with count_hash() function, return tally
            self.consumed += kmer.len() as u64; // Add kmer len to total consumed bases

            if self.store_kmers {
                // Get the canonical k-mer and store it 2-bit packed.
                let canonical_kmer = self.canon(kmer)?;
                self.hash_to_kmer
                    .as_mut()
                    .unwrap()
                    .insert(hashval, PackedKmer::encode(canonical_kmer.as_bytes()));
            }

            Ok(count) // Return the current total count for the hash
        }
    }

    /// Retrieve the count of a k-mer.
    pub fn get(&self, kmer: &str) -> PyResult<u64> {
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
    pub fn drop(&mut self, kmer: &str) -> PyResult<()> {
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

    /// Save the KmerCountTable to a gzip-compressed binary file.
    ///
    /// The format is `SAVE_MAGIC` followed by a gzipped `bincode` payload, which
    /// is both smaller and faster to read/write than the previous gzip-JSON.
    /// [`load`](Self::load) still reads tables written by older versions, but
    /// files written here are not readable by them.
    pub fn save(&self, filepath: &str) -> PyResult<()> {
        // Serialize to the compact binary encoding first.
        let payload = bincode::serialize(self)
            .map_err(|e| PyValueError::new_err(format!("Serialization error: {}", e)))?;

        // Write the format magic straight to the file, then gzip the payload
        // after it. (Writing the magic to the raw file first — rather than to a
        // buffered writer that is then handed to niffler — guarantees it lands
        // ahead of the gzip stream.)
        let mut file = File::create(filepath).map_err(|e| PyIOError::new_err(e.to_string()))?;
        file.write_all(SAVE_MAGIC)
            .map_err(|e| PyIOError::new_err(e.to_string()))?;

        let mut gz = get_writer(Box::new(file), Format::Gzip, niffler::level::Level::One)
            .map_err(|e| PyIOError::new_err(e.to_string()))?;
        gz.write_all(&payload)
            .map_err(|e| PyIOError::new_err(e.to_string()))?;
        gz.flush().map_err(|e| PyIOError::new_err(e.to_string()))?;

        Ok(())
    }

    #[staticmethod]
    /// Load a KmerCountTable saved with [`save`](Self::save).
    ///
    /// The format is auto-detected: files carrying `SAVE_MAGIC` are read as
    /// gzipped bincode, while anything else falls back to the legacy path
    /// (niffler-decompressed JSON), so gzip-JSON tables written by older oxli
    /// versions still load.
    pub fn load(filepath: &str) -> Result<KmerCountTable> {
        // Read the whole file so we can sniff the format prefix.
        let mut file = File::open(filepath)?;
        let mut raw = Vec::new();
        file.read_to_end(&mut raw)?;

        let loaded_table: KmerCountTable = if raw.starts_with(SAVE_MAGIC) {
            // New format: gzipped bincode after the magic. Stream the
            // decompressed bytes straight into bincode (buffered) so the whole
            // uncompressed payload is never materialized in memory at once.
            let (reader, _format) =
                niffler::get_reader(Box::new(Cursor::new(&raw[SAVE_MAGIC.len()..])))?;
            bincode::deserialize_from(std::io::BufReader::new(reader))
                .map_err(|e| anyhow::anyhow!("Deserialization error: {}", e))?
        } else {
            // Legacy format: niffler auto-detects gzip/plain, then parse JSON.
            let (mut reader, _format) = niffler::get_reader(Box::new(Cursor::new(raw)))?;
            let mut decompressed_data = String::new();
            reader.read_to_string(&mut decompressed_data)?;
            serde_json::from_str(&decompressed_data)
                .map_err(|e| anyhow::anyhow!("Deserialization error: {}", e))?
        };

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

        // Collect (decoded canonical k-mer, count) pairs, skipping hashes absent
        // from the counts table. Decoding runs in parallel via Rayon.
        let mut kmer_count_pairs: Vec<(String, u64)> = self
            .hash_to_kmer
            .as_ref()
            .unwrap()
            .par_iter()
            .filter_map(|(&hash, kmer)| self.counts.get(&hash).map(|&count| (kmer.decode(), count)))
            .collect();

        // Handle sorting based on the flags
        if sortkeys {
            // Sort by canonical kmer lexicographically
            kmer_count_pairs.par_sort_by(|a, b| a.0.cmp(&b.0));
        } else if sortcounts {
            // Sort by count, secondary sort by kmer
            kmer_count_pairs.par_sort_by(|a, b| a.1.cmp(&b.1).then_with(|| a.0.cmp(&b.0)));
        }
        // If both sortcounts and sortkeys are false, no sorting is done.

        // If a file is provided, write to the file
        if let Some(filepath) = file {
            let f = File::create(filepath)?;
            let mut writer = BufWriter::new(f);

            // Write each kmer:count pair to the file
            for (kmer, count) in &kmer_count_pairs {
                writeln!(writer, "{}\t{}", kmer, count)?;
            }

            writer.flush()?; // Ensure all data is written to the file
            Ok(vec![]) // Return an empty vector when writing to a file
        } else {
            // Return the vector of (kmer, count) tuples
            Ok(kmer_count_pairs)
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

    /// Internal engine shared by `consume` and `consume_file`.
    ///
    /// Counts the canonical k-mers of the raw sequence bytes `seq` and returns
    /// the number of k-mers consumed. When `store_kmers` is enabled on the
    /// table, the hash -> canonical-k-mer mapping is populated as a side effect
    /// (driven by `KmersAndHashesIter`); otherwise the faster `SeqToHashes`
    /// path is used. `skip_bad_kmers` controls whether non-DNA k-mers are
    /// silently skipped (true) or raise an error (false).
    #[pyo3(signature = (seq, skip_bad_kmers=true))]
    fn consume_bytes(&mut self, py: Python<'_>, seq: &[u8], skip_bad_kmers: bool) -> PyResult<u64> {
        // Release the GIL: the counting engine touches only Rust-owned data, so
        // other Python threads can run while it works.
        py.detach(|| self.count_seq_bytes(seq, skip_bad_kmers))
            .map_err(SeqError::into_pyerr)
    }

    // Consume this DNA string. Return total number of k-mers consumed.
    // If "skip_bad_kmers = true" then ignore kmers with non-DNA characters
    // else if "false" consume kmers until a bad kmer is encountered, then
    // exit with error.
    #[pyo3(signature = (seq, skip_bad_kmers=true))]
    pub fn consume(&mut self, py: Python<'_>, seq: &str, skip_bad_kmers: bool) -> PyResult<u64> {
        // Thin GIL-releasing wrapper over the shared byte-level engine.
        py.detach(|| self.count_seq_bytes(seq.as_bytes(), skip_bad_kmers))
            .map_err(SeqError::into_pyerr)
    }

    /// Consume all sequences from a FASTA/FASTQ file. Return total k-mers consumed.
    ///
    /// The file is parsed with needletail, which auto-detects the format and
    /// transparently decompresses gzip/bzip2/xz inputs. Each record is
    /// normalized (upper-cased, line breaks stripped) and its k-mers counted
    /// via the shared engine, so `store_kmers` is honoured when enabled.
    /// `skip_bad_kmers` behaves as in `consume`.
    ///
    /// The GIL is released around the k-mer counting of each record (needletail's
    /// reader is not `Send`, so parsing itself stays on the calling thread).
    #[pyo3(signature = (filename, skip_bad_kmers=true))]
    pub fn consume_file(
        &mut self,
        py: Python<'_>,
        filename: &str,
        skip_bad_kmers: bool,
    ) -> PyResult<u64> {
        // Total k-mers consumed across all records.
        let mut n: u64 = 0;
        // Number of records processed (for logging).
        let mut n_records: u64 = 0;

        // Open the file; needletail auto-detects the format and compression.
        let mut reader = parse_fastx_file(filename)
            .map_err(|e| PyIOError::new_err(format!("failed to open '{}': {}", filename, e)))?;

        // Iterate over each sequence record in the file.
        while let Some(record) = reader.next() {
            let record = record.map_err(|e| {
                PyValueError::new_err(format!("invalid record in '{}': {}", filename, e))
            })?;

            // Normalise the sequence (upper-case, strip newlines/whitespace),
            // then count its k-mers with the GIL released, propagating any
            // bad-k-mer error.
            let normseq = record.normalize(false);
            let seq_bytes: &[u8] = normseq.as_ref();
            n += py
                .detach(|| self.count_seq_bytes(seq_bytes, skip_bad_kmers))
                .map_err(SeqError::into_pyerr)?;
            n_records += 1;
        }

        debug!(
            "consume_file: processed {} record(s), {} k-mer(s) from '{}'",
            n_records, n, filename
        );

        Ok(n)
    }

    /// Consume a DNA string in parallel by splitting it into overlapping chunks,
    /// processing each chunk concurrently using Rayon, and merging the results.
    ///
    /// Each chunk overlaps its neighbours by `ksize - 1` bases so that no k-mer
    /// spanning a chunk boundary is missed.  After all chunks have been processed
    /// the per-chunk tables are merged into `self` using a fast serial merge.
    ///
    /// # Arguments
    /// * `seq`            - The DNA sequence to consume.
    /// * `chunk_size`     - Target number of k-mers per chunk (clamped to at
    ///                      least `ksize`).  Defaults to 50 000.
    /// * `skip_bad_kmers` - If `true`, k-mers containing non-DNA characters are
    ///                      silently skipped.  If `false`, the first bad k-mer
    ///                      raises an error.  Defaults to `true`.
    ///
    /// # Returns
    /// The total number of k-mers consumed (identical to what `consume` would
    /// return for the same sequence).
    #[pyo3(signature = (seq, chunk_size=50_000, skip_bad_kmers=true))]
    pub fn parallel_consume(
        &mut self,
        py: Python<'_>,
        seq: &str,
        chunk_size: usize,
        skip_bad_kmers: bool,
    ) -> PyResult<u64> {
        let ksize = self.ksize as usize;
        let seq_len = seq.len();

        // Nothing to do for sequences shorter than k.
        if seq_len < ksize {
            self.consumed += seq_len as u64;
            return Ok(0);
        }

        // Clamp chunk_size so it is always >= ksize.
        let chunk_size = max(chunk_size, ksize);
        let this_ksize = self.ksize;
        let store_kmers = self.store_kmers;
        let seq_bytes = seq.as_bytes();

        // All counting (and the Rayon fan-out) runs with the GIL released.
        py.detach(|| -> Result<u64, SeqError> {
            // For short sequences that fit in a single chunk, count serially.
            if seq_len <= chunk_size {
                return self.count_seq_bytes(seq_bytes, skip_bad_kmers);
            }

            // Build a list of (start, end) byte-index pairs for each chunk.
            // Adjacent chunks overlap by (ksize - 1) bases so that every k-mer
            // crossing a chunk boundary appears in exactly one chunk.
            let mut coord_pairs: Vec<(usize, usize)> = Vec::new();
            let mut start = 0;
            while start < seq_len {
                let end = (start + chunk_size + ksize - 1).min(seq_len);
                coord_pairs.push((start, end));
                if end == seq_len {
                    break;
                }
                start += chunk_size;
            }

            // Process chunks in parallel: each chunk produces a local table.
            let chunk_results: Vec<Result<(KmerCountTable, u64), SeqError>> = coord_pairs
                .into_par_iter()
                .map(|(start, end)| {
                    let mut t = KmerCountTable::new(this_ksize, store_kmers);
                    let n = t.count_seq_bytes(&seq_bytes[start..end], skip_bad_kmers)?;
                    Ok((t, n))
                })
                .collect();

            // Merge the per-chunk tables into self and accumulate the k-mer count.
            let mut total_n: u64 = 0;
            for result in chunk_results {
                let (t, n) = result?;
                self._merge(t);
                total_n += n;
            }

            // Record the total bases processed (full sequence, counted once).
            self.consumed += seq_len as u64;

            Ok(total_n)
        })
        .map_err(SeqError::into_pyerr)
    }

    // Helper method to get hash set of k-mers
    fn hash_set(&self) -> HashSet<u64> {
        self.counts.keys().cloned().collect()
    }

    // Set operation methods.
    //
    // These iterate the count maps directly instead of first materializing two
    // intermediate `HashSet`s (as the old `self.hash_set()` / `other.hash_set()`
    // approach did), halving the allocation and copying work.
    pub fn union(&self, other: &KmerCountTable) -> HashSet<u64> {
        let mut result: HashSet<u64> =
            HashSet::with_capacity(self.counts.len() + other.counts.len());
        result.extend(self.counts.keys().copied());
        result.extend(other.counts.keys().copied());
        result
    }

    pub fn intersection(&self, other: &KmerCountTable) -> HashSet<u64> {
        // Probe the larger map while iterating the smaller one.
        let (small, large) = if self.counts.len() <= other.counts.len() {
            (&self.counts, &other.counts)
        } else {
            (&other.counts, &self.counts)
        };
        small
            .keys()
            .filter(|k| large.contains_key(k))
            .copied()
            .collect()
    }

    pub fn difference(&self, other: &KmerCountTable) -> HashSet<u64> {
        self.counts
            .keys()
            .filter(|k| !other.counts.contains_key(k))
            .copied()
            .collect()
    }

    pub fn symmetric_difference(&self, other: &KmerCountTable) -> HashSet<u64> {
        let mut result: HashSet<u64> = self
            .counts
            .keys()
            .filter(|k| !other.counts.contains_key(k))
            .copied()
            .collect();
        result.extend(
            other
                .counts
                .keys()
                .filter(|k| !self.counts.contains_key(k))
                .copied(),
        );
        result
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
    fn __getitem__(&self, kmer: &str) -> PyResult<u64> {
        self.get(kmer)
    }

    // Python dunder method for __setitem__
    pub fn __setitem__(&mut self, kmer: &str, count: u64) -> PyResult<()> {
        // Calculate the hash for the k-mer
        let hashval = self.hash_kmer(kmer)?;
        // Set the count for the k-mer
        self.counts.insert(hashval, count);
        Ok(())
    }

    #[pyo3(signature = (seq, skip_bad_kmers=true))]
    pub fn kmers_and_hashes(
        &self,
        py: Python<'_>,
        seq: &str,
        skip_bad_kmers: bool,
    ) -> PyResult<Vec<(String, u64)>> {
        let ksize = self.ksize as usize;
        // Build the (canonical_kmer, hash) list with the GIL released.
        py.detach(|| -> Result<Vec<(String, u64)>, SeqError> {
            KmersAndHashesIter::new(seq, ksize, skip_bad_kmers).collect()
        })
        .map_err(SeqError::into_pyerr)
    }

    /// Calculates the Jaccard Similarity Coefficient between two KmerCountTable objects.
    /// # Returns
    /// The Jaccard Similarity Coefficient between the two tables as a float value between 0 and 1.
    pub fn jaccard(&self, other: &KmerCountTable) -> f64 {
        // Single pass: count the intersection by probing the larger map while
        // iterating the smaller, then derive the union size arithmetically
        // (|A| + |B| - |A ∩ B|). Allocates no intermediate sets.
        let (small, large) = if self.counts.len() <= other.counts.len() {
            (&self.counts, &other.counts)
        } else {
            (&other.counts, &self.counts)
        };
        let intersection_size = small.keys().filter(|k| large.contains_key(k)).count();
        let union_size = self.counts.len() + other.counts.len() - intersection_size;

        // Two empty sets are identical by convention.
        if union_size == 0 {
            return 1.0;
        }

        intersection_size as f64 / union_size as f64
    }

    /// Cosine similarity between two `KmerCountTable` objects.
    /// # Returns
    /// The cosine similarity between the two tables as a float value between 0 and 1.
    pub fn cosine(&self, py: Python<'_>, other: &KmerCountTable) -> f64 {
        // Early return if either table is empty.
        if self.counts.is_empty() || other.counts.is_empty() {
            return 0.0;
        }

        // The dot product and magnitudes are Rayon-parallel; release the GIL so
        // other Python threads can run during the computation.
        py.detach(|| {
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
        })
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

        // A serial merge. The previous version fanned out with Rayon but guarded
        // the single shared counts map with one global `Mutex`, so every update
        // contended on the same lock — slower than a straight serial merge and
        // not worth the thread-dispatch overhead. Merging is a cheap hash-probe
        // per entry, so we do it directly.
        let mut total_added: u64 = 0;
        let mut new_keys: u64 = 0;

        self.counts.reserve(other.counts.len());
        for (&hash, &count) in &other.counts {
            let current = self.counts.entry(hash).or_insert(0);
            if *current == 0 {
                new_keys += 1;
            }
            *current += count;
            total_added += count;
        }

        self.consumed += other.consumed;

        if self.store_kmers {
            if other.store_kmers {
                let my_map = self.hash_to_kmer.as_mut().unwrap();
                let other_map = other.hash_to_kmer.as_ref().unwrap();
                my_map.reserve(other_map.len());
                for (&hash, kmer) in other_map {
                    my_map.entry(hash).or_insert_with(|| kmer.clone());
                }
            } else {
                // Kept on stderr (not the logger) so it always surfaces, matching
                // the historical behaviour the test-suite pins.
                eprintln!(
                    "Warning: Incoming table does not store k-mers, but target table does. \
                     K-mer information for new hashes will be missing."
                );
            }
        }

        debug!(
            "add: {} k-mer counts merged, {} new keys added",
            total_added, new_keys
        );

        Ok((total_added, new_keys))
    }
}

// Private (non-Python-visible) methods
impl KmerCountTable {
    /// Pure-Rust counting engine shared by `consume`, `consume_bytes`,
    /// `consume_file`, and `parallel_consume`.
    ///
    /// Counts the canonical k-mers of the raw sequence bytes `seq` and returns
    /// the number of k-mers consumed. When `store_kmers` is enabled the
    /// hash -> canonical-k-mer mapping is populated as a side effect (driven by
    /// `KmersAndHashesIter`); otherwise the faster `SeqToHashes` path is used.
    /// `skip_bad_kmers` controls whether non-DNA k-mers are silently skipped
    /// (true) or raise an error (false).
    ///
    /// Holds no GIL-bound state, so callers run it inside
    /// `Python::detach`; errors are the plain-Rust [`SeqError`].
    fn count_seq_bytes(&mut self, seq: &[u8], skip_bad_kmers: bool) -> Result<u64, SeqError> {
        // Raw bytes processed (added to the `consumed` tracker below).
        let new_len = seq.len();
        // Running tally of k-mers consumed.
        let mut n: u64 = 0;

        // Pre-size the counts map by the number of k-mer windows (an upper bound
        // on the distinct keys this call can add), capped so a single very long
        // sequence cannot request a pathologically large allocation. This avoids
        // repeated rehashing as the map fills.
        let est_new = new_len
            .saturating_sub(self.ksize as usize)
            .saturating_add(1)
            .min(1 << 20);
        self.counts.reserve(est_new);

        // If store_kmers is true, count & record hash:kmer pairs.
        if self.store_kmers {
            // KmersAndHashesIter works on &str, so the bytes must be valid UTF-8.
            // needletail's `normalize` yields ASCII, so this only fails on
            // genuinely malformed input.
            let seq_str = std::str::from_utf8(seq).map_err(|_| SeqError::NonUtf8)?;
            self.hash_to_kmer.as_mut().unwrap().reserve(est_new);

            // Iterate over (canonical_kmer, hash) pairs.
            let iter = KmersAndHashesIter::new(seq_str, self.ksize as usize, skip_bad_kmers);
            for result in iter {
                let (kmer, hash) = result?;
                if hash != 0 {
                    self.hash_to_kmer
                        .as_mut()
                        .unwrap()
                        .insert(hash, PackedKmer::encode(kmer.as_bytes()));
                    *self.counts.entry(hash).or_insert(0) += 1;
                    n += 1;
                }
            }
        } else {
            // Fast path: hash and count k-mers directly.
            let hashes = SeqToHashes::new(
                seq,
                self.ksize.into(),
                skip_bad_kmers,
                false,
                HashFunctions::Murmur64Dna,
                42,
            )
            .expect("Failed to create SeqToHashes");

            for hash_value in hashes {
                match hash_value {
                    Ok(0) => continue,
                    Ok(x) => {
                        self.count_hash(x);
                    }
                    Err(_) => return Err(SeqError::BadKmer(n)),
                }
                n += 1;
            }
        }

        // Update the total sequence consumed tracker.
        self.consumed += new_len as u64;

        Ok(n)
    }

    /// Merge `other` into `self` by summing counts for shared hashes and
    /// inserting new hashes.  `consumed` is intentionally **not** propagated
    /// because callers that split sequences into chunks track `consumed`
    /// at the top level.
    fn _merge(&mut self, other: KmerCountTable) {
        for (hashval, count) in other.counts {
            *self.counts.entry(hashval).or_insert(0) += count;
        }

        if self.store_kmers {
            if let Some(other_map) = other.hash_to_kmer {
                let my_map = self.hash_to_kmer.as_mut().unwrap();
                for (hash, kmer) in other_map {
                    my_map.entry(hash).or_insert(kmer);
                }
            }
        }
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

struct KmersAndHashesIter {
    seq: String,          // The sequence to iterate over
    seq_rc: String,       // reverse complement sequence
    ksize: usize,         // K-mer size
    pos: usize,           // Current position in the sequence
    end: usize,           // The end position for k-mer extraction
    hasher: SeqToHashes,  // Iterator for generating hashes
    skip_bad_kmers: bool, // Flag to skip bad k-mers
}

impl KmersAndHashesIter {
    fn new(seq: &str, ksize: usize, skip_bad_kmers: bool) -> Self {
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
        )
        .expect("Failed to create SeqToHashes");

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
    type Item = Result<(String, u64), SeqError>;

    fn next(&mut self) -> Option<Self::Item> {
        // Loop (rather than recurse) so a long run of bad k-mers — e.g. a
        // stretch of `N`s — cannot overflow the stack when skipping.
        loop {
            // Stop once every window has been visited.
            if self.pos >= self.end {
                return None;
            }

            let start = self.pos;
            let ksize = self.ksize;
            let rpos = self.end - start - 1;

            // Extract the current k-mer and its reverse complement.
            let substr = &self.seq[start..start + ksize];
            let substr_rc = &self.seq_rc[rpos..rpos + ksize];

            // Get the next hash value from the hasher, advancing position.
            let hashval = self.hasher.next().expect("should not run out of hashes");
            self.pos += 1;

            match hashval {
                // Good k-mer: return its canonical form (lexicographically
                // smaller of the forward and reverse-complement windows).
                Ok(h) if h > 0 => {
                    let canonical_kmer = if substr < substr_rc {
                        substr
                    } else {
                        substr_rc
                    };
                    return Some(Ok((canonical_kmer.to_string(), h)));
                }
                // Bad k-mer (hash 0): warn, then skip or emit a ("", 0) sentinel.
                Ok(_) => {
                    eprintln!("bad k-mer at position {}: {}", start + 1, substr);
                    if self.skip_bad_kmers {
                        continue; // advance to the next window
                    }
                    return Some(Ok((String::new(), 0)));
                }
                // Error raised by SeqToHashes. With `force = true` (set in `new`)
                // bad k-mers surface as hash 0 above rather than here, so this
                // branch is effectively unreachable, but we map it faithfully.
                Err(_) => return Some(Err(SeqError::BadKmer(start as u64))),
            }
        }
    }
}

/// Shared module setup for both the GIL and free-threaded module entry points.
fn register_oxli(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // `try_init` (unlike `init`) does not panic if a logger is already set,
    // which can happen when the extension is imported more than once.
    let _ = env_logger::try_init();
    m.add_class::<KmerCountTable>()?;
    Ok(())
}

// Free-threaded (no-GIL) build: declare the module safe to import without
// re-enabling the GIL. This is sound because all mutable state lives on
// `#[pyclass]` instances, whose method borrows PyO3 serializes per object; the
// module holds no shared mutable state.
//
// The `gil_used = false` option emits a `Py_mod_gil` module slot, which is not
// available under the abi3 (limited API) floor the GIL wheels are built against
// — so it is applied only here, to the free-threaded build (`Py_GIL_DISABLED`),
// which is never an abi3 build. Applying it unconditionally makes the abi3 wheel
// segfault on import.
#[cfg(Py_GIL_DISABLED)]
#[pymodule(gil_used = false)]
fn oxli(m: &Bound<'_, PyModule>) -> PyResult<()> {
    register_oxli(m)
}

// Standard (GIL) build, including abi3 wheels.
#[cfg(not(Py_GIL_DISABLED))]
#[pymodule]
fn oxli(m: &Bound<'_, PyModule>) -> PyResult<()> {
    register_oxli(m)
}
