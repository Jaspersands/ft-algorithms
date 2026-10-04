//! Many shots: seeded per shot, optionally in parallel, records bit-packed.

use crate::dense::{Dense, MAX_DENSE_QUBITS};
use crate::exec::{run_shot, RunApply, ShotOutput};
use crate::noise::{FaultPlan, SiteTree};
use crate::program::Program;
use crate::rng::shot_rng;
use crate::sparse::Sparse;
use crate::state::State;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Backend {
    /// Dense up to `AUTO_DENSE_QUBITS` qubits, sparse above.
    Auto,
    Dense,
    Sparse,
}

pub const AUTO_DENSE_QUBITS: u32 = 12;

impl Backend {
    pub fn parse(s: &str) -> Result<Backend, String> {
        match s {
            "auto" => Ok(Backend::Auto),
            "dense" => Ok(Backend::Dense),
            "sparse" => Ok(Backend::Sparse),
            _ => Err(format!("unknown backend '{s}' (want auto, dense or sparse)")),
        }
    }

    pub fn resolve(self, num_qubits: u32) -> Result<Backend, String> {
        match self {
            Backend::Auto => Ok(if num_qubits <= AUTO_DENSE_QUBITS { Backend::Dense } else { Backend::Sparse }),
            Backend::Dense if num_qubits > MAX_DENSE_QUBITS => {
                Err(format!("{num_qubits} qubits is too many for the dense backend (at most {MAX_DENSE_QUBITS})"))
            }
            b => Ok(b),
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct BatchResult {
    pub shots: u64,
    pub bytes_per_shot: usize,
    /// shots × bytes_per_shot; record j of a shot is bit j % 8 of byte j / 8.
    pub records: Vec<u8>,
    pub fault_shot: Vec<u64>,
    pub fault_site: Vec<u64>,
    pub fault_class: Vec<u16>,
    pub fault_pauli: Vec<u32>,
}

fn pack(records: &[bool], out: &mut [u8]) {
    for (j, &b) in records.iter().enumerate() {
        if b {
            out[j / 8] |= 1 << (j % 8);
        }
    }
}

fn run_range<S: State + RunApply>(p: &Program, sites: &SiteTree, seed: u64, plan: FaultPlan, shots: std::ops::Range<u64>, log: bool) -> Result<Vec<ShotOutput>, String> {
    let mut out = Vec::with_capacity((shots.end - shots.start) as usize);
    for s in shots {
        let mut rng = shot_rng(seed, s);
        let (o, _) = run_shot::<S>(p, sites, plan, &mut rng, log)?;
        out.push(o);
    }
    Ok(out)
}

/// Runs shots `first_shot .. first_shot + shots`; shot s uses RNG stream (seed, s), so any split
/// of a job into calls, and any thread count, gives the same records.
#[allow(clippy::too_many_arguments)]
pub fn sample(
    p: &Program,
    sites: &SiteTree,
    shots: u64,
    seed: u64,
    first_shot: u64,
    plan: FaultPlan,
    backend: Backend,
    threads: usize,
    log: bool,
) -> Result<BatchResult, String> {
    let backend = backend.resolve(p.num_qubits)?;
    if let FaultPlan::Exactly(k) = plan {
        if k as u64 > sites.num_sites() {
            return Err(format!("cannot place {k} faults: the program has {} noise sites", sites.num_sites()));
        }
    }
    let run = |r: std::ops::Range<u64>| -> Result<Vec<ShotOutput>, String> {
        match backend {
            Backend::Dense => run_range::<Dense>(p, sites, seed, plan, r, log),
            _ => run_range::<Sparse>(p, sites, seed, plan, r, log),
        }
    };
    let end = first_shot.checked_add(shots).ok_or("shot range overflows")?;
    let outputs: Vec<ShotOutput> = run_parallel(first_shot..end, threads, &run)?;

    let bps = (p.num_measurements as usize).div_ceil(8);
    let mut res = BatchResult { shots, bytes_per_shot: bps, records: vec![0u8; bps * shots as usize], ..Default::default() };
    for (i, o) in outputs.iter().enumerate() {
        pack(&o.records, &mut res.records[i * bps..(i + 1) * bps]);
        for f in &o.faults {
            res.fault_shot.push(first_shot + i as u64);
            res.fault_site.push(f.site);
            res.fault_class.push(f.class);
            res.fault_pauli.push(f.pauli);
        }
    }
    Ok(res)
}

#[cfg(feature = "parallel")]
fn run_parallel<F>(range: std::ops::Range<u64>, threads: usize, run: &F) -> Result<Vec<ShotOutput>, String>
where
    F: Fn(std::ops::Range<u64>) -> Result<Vec<ShotOutput>, String> + Sync,
{
    use rayon::prelude::*;
    let n = range.end - range.start;
    let workers = if threads == 0 { rayon::current_num_threads() } else { threads };
    if workers <= 1 || n < 2 {
        return run(range);
    }
    // Several chunks per worker for balance; chunk boundaries never affect results.
    let chunk = n.div_ceil(workers as u64 * 8).max(1);
    let chunks: Vec<std::ops::Range<u64>> =
        (0..n.div_ceil(chunk)).map(|c| range.start + c * chunk..(range.start + (c + 1) * chunk).min(range.end)).collect();
    let pool = rayon::ThreadPoolBuilder::new().num_threads(workers).build().map_err(|e| e.to_string())?;
    let parts: Vec<Result<Vec<ShotOutput>, String>> = pool.install(|| chunks.into_par_iter().map(run).collect());
    let mut out = Vec::with_capacity(n as usize);
    for p in parts {
        out.extend(p?);
    }
    Ok(out)
}

#[cfg(not(feature = "parallel"))]
fn run_parallel<F>(range: std::ops::Range<u64>, _threads: usize, run: &F) -> Result<Vec<ShotOutput>, String>
where
    F: Fn(std::ops::Range<u64>) -> Result<Vec<ShotOutput>, String>,
{
    run(range)
}

/// The final state of one noiseless shot (records drawn with stream (seed, 0)).
pub fn final_state(p: &Program, seed: u64, backend: Backend) -> Result<(Vec<(u128, crate::C64)>, Vec<bool>), String> {
    let sites = SiteTree::build(p);
    let mut rng = shot_rng(seed, 0);
    match backend.resolve(p.num_qubits)? {
        Backend::Dense => {
            let (o, s) = run_shot::<Dense>(p, &sites, FaultPlan::None, &mut rng, false)?;
            Ok((s.amplitudes(), o.records))
        }
        _ => {
            let (o, s) = run_shot::<Sparse>(p, &sites, FaultPlan::None, &mut rng, false)?;
            Ok((s.amplitudes(), o.records))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn noisy_program() -> Program {
        Program::parse(
            "H 0 1 2 3\nREPEAT 20 {\nCX 0 1\nT 1\nCCX 1 2 3\nH 2\nDEPOLARIZE1(0.05) 0 1 2 3\nDEPOLARIZE2(0.02) 0 3\nM 0\n\
             TABLE rec[-1] {S 2} {\nH 2\nX_ERROR(0.1) 1\n}\n}\nM 1 2 3",
        )
        .unwrap()
    }

    #[test]
    fn independent_of_threads_and_splits() {
        let p = noisy_program();
        let t = SiteTree::build(&p);
        let a = sample(&p, &t, 300, 42, 0, FaultPlan::Plain, Backend::Sparse, 1, true).unwrap();
        for threads in [2, 7] {
            let b = sample(&p, &t, 300, 42, 0, FaultPlan::Plain, Backend::Sparse, threads, true).unwrap();
            assert_eq!(a.records, b.records);
            assert_eq!(a.fault_site, b.fault_site);
            assert_eq!(a.fault_pauli, b.fault_pauli);
        }
        let first = sample(&p, &t, 100, 42, 0, FaultPlan::Plain, Backend::Sparse, 3, false).unwrap();
        let rest = sample(&p, &t, 200, 42, 100, FaultPlan::Plain, Backend::Sparse, 3, false).unwrap();
        let mut joined = first.records.clone();
        joined.extend(rest.records);
        assert_eq!(joined, a.records);
    }

    #[test]
    fn dense_equals_sparse() {
        let p = noisy_program();
        let t = SiteTree::build(&p);
        for plan in [FaultPlan::None, FaultPlan::Plain, FaultPlan::Exactly(3)] {
            let a = sample(&p, &t, 200, 7, 0, plan, Backend::Sparse, 0, true).unwrap();
            let b = sample(&p, &t, 200, 7, 0, plan, Backend::Dense, 0, true).unwrap();
            assert_eq!(a.records, b.records, "{plan:?}");
            assert_eq!(a.fault_site, b.fault_site);
        }
    }

    #[test]
    fn packing_and_errors() {
        let p = Program::parse("X 0 2\nM 0 1 2 0 1 2 0 1 2").unwrap();
        let t = SiteTree::build(&p);
        let r = sample(&p, &t, 2, 0, 0, FaultPlan::None, Backend::Auto, 0, false).unwrap();
        assert_eq!(r.bytes_per_shot, 2);
        assert_eq!(r.records, vec![0b0110_1101, 0b1, 0b0110_1101, 0b1]);
        assert!(sample(&p, &t, 1, 0, 0, FaultPlan::Exactly(1), Backend::Auto, 0, false).is_err());
        let big = Program::parse("H 40").unwrap();
        assert!(Backend::Dense.resolve(big.num_qubits).is_err());
    }
}
