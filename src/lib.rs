use pyo3::prelude::*;

use sourmash::sketch::nodegraph::Nodegraph;
use sourmash::_hash_murmur;

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
    Ok(())
}
