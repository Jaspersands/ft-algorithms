//! Python bindings: `ftalgo._ftsim`. The thin Python layer in `python/ftalgo/engine.py` turns
//! the byte buffers into numpy arrays.

use crate::batch::{final_state, sample, Backend};
use crate::noise::{FaultPlan, SiteTree};
use crate::program::Program;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyDict};

fn err(e: String) -> PyErr {
    PyValueError::new_err(e)
}

#[pyclass(name = "Program", module = "ftalgo._ftsim", frozen)]
pub struct PyProgram {
    program: Program,
    sites: SiteTree,
}

#[pymethods]
impl PyProgram {
    #[new]
    fn new(text: &str) -> PyResult<Self> {
        let program = Program::parse(text).map_err(err)?;
        let sites = SiteTree::build(&program);
        Ok(PyProgram { program, sites })
    }

    #[getter]
    fn num_qubits(&self) -> u32 {
        self.program.num_qubits
    }

    #[getter]
    fn num_measurements(&self) -> u64 {
        self.program.num_measurements
    }

    #[getter]
    fn num_sites(&self) -> u64 {
        self.sites.num_sites()
    }

    #[getter]
    fn expected_faults(&self) -> f64 {
        self.sites.expected_faults()
    }

    #[getter]
    fn classes(&self) -> Vec<String> {
        self.program.classes.clone()
    }

    /// [P(K = 0), …, P(K = kmax), P(K > kmax)] for the number K of static faults.
    fn fault_count_distribution(&self, kmax: usize) -> PyResult<Vec<f64>> {
        if kmax > 10_000 {
            return Err(err("kmax at most 10000".into()));
        }
        Ok(self.sites.fault_count_distribution(kmax))
    }

    /// Returns a dict of raw buffers: records (shots × bytes_per_shot), and with log_faults the
    /// fault arrays (u64 shot, u64 site, u16 class, u32 pauli; native endian).
    #[pyo3(signature = (shots, *, seed, faults=None, backend="auto", threads=0, first_shot=0, log_faults=false))]
    #[allow(clippy::too_many_arguments)]
    fn sample<'py>(
        &self,
        py: Python<'py>,
        shots: u64,
        seed: u64,
        faults: Option<Bound<'py, PyAny>>,
        backend: &str,
        threads: usize,
        first_shot: u64,
        log_faults: bool,
    ) -> PyResult<Bound<'py, PyDict>> {
        let plan = match faults {
            None => FaultPlan::Plain,
            Some(f) => {
                if let Ok(s) = f.extract::<String>() {
                    match s.as_str() {
                        "plain" => FaultPlan::Plain,
                        "none" => FaultPlan::None,
                        _ => return Err(err(format!("faults must be 'plain', 'none' or an integer, not '{s}'"))),
                    }
                } else if let Ok(k) = f.extract::<usize>() {
                    FaultPlan::Exactly(k)
                } else {
                    return Err(err("faults must be 'plain', 'none' or a non-negative integer".into()));
                }
            }
        };
        let backend = Backend::parse(backend).map_err(err)?;
        let (p, t) = (&self.program, &self.sites);
        let r = py.detach(|| sample(p, t, shots, seed, first_shot, plan, backend, threads, log_faults)).map_err(err)?;
        let d = PyDict::new(py);
        d.set_item("shots", r.shots)?;
        d.set_item("bytes_per_shot", r.bytes_per_shot)?;
        d.set_item("records", PyBytes::new(py, &r.records))?;
        if log_faults {
            let b64 = |v: &[u64]| v.iter().flat_map(|x| x.to_ne_bytes()).collect::<Vec<u8>>();
            d.set_item("fault_shot", PyBytes::new(py, &b64(&r.fault_shot)))?;
            d.set_item("fault_site", PyBytes::new(py, &b64(&r.fault_site)))?;
            let c: Vec<u8> = r.fault_class.iter().flat_map(|x| x.to_ne_bytes()).collect();
            d.set_item("fault_class", PyBytes::new(py, &c))?;
            let q: Vec<u8> = r.fault_pauli.iter().flat_map(|x| x.to_ne_bytes()).collect();
            d.set_item("fault_pauli", PyBytes::new(py, &q))?;
        }
        Ok(d)
    }

    /// One noiseless shot: (indices as decimal strings, re, im, records).
    #[pyo3(signature = (*, seed=0, backend="auto"))]
    fn final_state(&self, py: Python<'_>, seed: u64, backend: &str) -> PyResult<(Vec<u128>, Vec<f64>, Vec<f64>, Vec<bool>)> {
        let backend = Backend::parse(backend).map_err(err)?;
        let p = &self.program;
        let (amps, recs) = py.detach(|| final_state(p, seed, backend)).map_err(err)?;
        let idx = amps.iter().map(|a| a.0).collect();
        let re = amps.iter().map(|a| a.1.re).collect();
        let im = amps.iter().map(|a| a.1.im).collect();
        Ok((idx, re, im, recs))
    }
}

#[pymodule]
fn _ftsim(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyProgram>()?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
