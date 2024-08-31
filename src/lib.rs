use pyo3::prelude::*;

use std::collections::HashMap;

use sourmash::sketch::nodegraph::Nodegraph;
use sourmash::_hash_murmur;

#[pyclass]
struct KmerCountTable {
    counts: HashMap<u64, usize>,
}

#[pymethods]
impl KmerCountTable {
    #[new]
    pub fn new() -> Self {
        Self { counts: HashMap::new() }
    }

    pub fn count(&mut self, kmer: String) -> PyResult<usize> {
        let hashval = _hash_murmur(kmer.as_bytes(), 42);

        let mut count: usize = 1;
        if self.counts.contains_key(&hashval) {
            count = *self.counts.get(&hashval).unwrap();
            count = count + 1;
        }
        self.counts.insert(hashval, count);

        Ok(count)
    }

    pub fn get(&self, kmer: String) -> PyResult<usize> {
        let hashval = _hash_murmur(kmer.as_bytes(), 42);

        let count = match self.counts.get(&hashval) {
            Some(count) => count,
            None => &(0 as usize)
        };
        Ok(*count)
    }
}

/// Formats the sum of two numbers as string.
#[pyfunction]
fn sum_as_string(a: String) -> PyResult<usize> {
    let mut ng: Nodegraph = Nodegraph::with_tables(23, 6, 3);

    let hashval = _hash_murmur(a.as_bytes(), 42);
    ng.count(hashval);
    Ok(ng.get(hashval))
}

/// A Python module implemented in Rust.
#[pymodule]
fn oxli(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    m.add_class::<KmerCountTable>()?;
    Ok(())
}
