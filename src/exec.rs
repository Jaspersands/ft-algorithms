//! Running one shot: walk the program, apply gates, measure, and inject faults.
//!
//! RNG consumption order (part of the determinism contract, identical on every backend):
//! the static fault set is drawn first; then, in execution order, one f64 per measurement or
//! reset, one f64 per fired static site (to choose its Pauli), and one f64 per dynamic site
//! executed (to decide whether and how it fires).
//!
//! Runs of permutation-and-phase gates (X, Y, Z, S, T, CX, CY, CZ, SWAP, CCX, CCZ and fault
//! Paulis) are buffered and handed to the backend together, which lets the sparse backend apply a
//! whole run to each basis state in one pass.

use crate::gates::{Basis, Gate};
use crate::noise::{FaultPlan, SiteTree};
use crate::program::{Block, Op, Program};
use crate::rng::Xoshiro;
use crate::state::State;

/// A fault that fired. `site` is the static site index, or `u64::MAX` for a dynamic site.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Fault {
    pub site: u64,
    pub class: u16,
    /// Packed Pauli, two bits per target of the group (see `Noise`).
    pub pauli: u32,
}

#[derive(Clone, Debug, Default)]
pub struct ShotOutput {
    pub records: Vec<bool>,
    pub faults: Vec<Fault>,
}

/// One gate of a buffered run, targets padded to three.
pub type RunGate = (Gate, [u32; 3]);

#[inline]
fn is_permutation(g: Gate) -> bool {
    !matches!(g, Gate::H | Gate::SqrtX | Gate::SqrtXdg)
}

struct Runner<'a, S: State> {
    state: S,
    rng: &'a mut Xoshiro,
    records: Vec<bool>,
    fired: Vec<u64>,
    next: usize,
    site: u64,
    class: u16,
    dynamic_depth: u32,
    noisy: bool,
    log: bool,
    faults: Vec<Fault>,
    run: Vec<RunGate>,
}

impl<'a, S: State + RunApply> Runner<'a, S> {
    #[inline]
    fn push_gate(&mut self, g: Gate, t: &[u32]) {
        if g == Gate::I {
            return;
        }
        if is_permutation(g) {
            let mut a = [0u32; 3];
            a[..t.len()].copy_from_slice(t);
            self.run.push((g, a));
        } else {
            self.flush();
            self.state.apply(g, t);
        }
    }

    #[inline]
    fn flush(&mut self) {
        if !self.run.is_empty() {
            self.state.apply_run(&self.run);
            self.run.clear();
        }
    }

    fn push_pauli(&mut self, qs: &[u32], pauli: u32) {
        for (i, &q) in qs.iter().enumerate() {
            let g = match (pauli >> (2 * i)) & 3 {
                1 => Gate::X,
                2 => Gate::Y,
                3 => Gate::Z,
                _ => continue,
            };
            self.push_gate(g, &[q]);
        }
    }

    /// Measures `q` in Z, drawing one f64; outcomes of (numerically) zero probability never occur.
    fn measure_z(&mut self, q: u32) -> bool {
        self.flush();
        let p1 = self.state.prob_one(q);
        let u = self.rng.next_f64();
        let outcome = if p1 < 1e-13 {
            false
        } else if p1 > 1.0 - 1e-13 {
            true
        } else {
            u < p1
        };
        self.state.collapse(q, outcome);
        outcome
    }

    fn to_z(&mut self, b: Basis, q: u32) {
        match b {
            Basis::Z => {}
            Basis::X => self.push_gate(Gate::H, &[q]),
            Basis::Y => {
                self.push_gate(Gate::Sdg, &[q]);
                self.push_gate(Gate::H, &[q]);
            }
        }
    }

    fn from_z(&mut self, b: Basis, q: u32) {
        match b {
            Basis::Z => {}
            Basis::X => self.push_gate(Gate::H, &[q]),
            Basis::Y => {
                self.push_gate(Gate::H, &[q]);
                self.push_gate(Gate::S, &[q]);
            }
        }
    }

    fn block(&mut self, block: &Block) {
        for op in block {
            match op {
                Op::Gate(g, qs) => {
                    for t in qs.chunks(g.arity()) {
                        self.push_gate(*g, t);
                    }
                }
                Op::Measure(b, qs) => {
                    for &q in qs {
                        self.to_z(*b, q);
                        let m = self.measure_z(q);
                        self.records.push(m);
                        self.from_z(*b, q);
                    }
                }
                Op::Reset(b, qs) => {
                    for &q in qs {
                        if self.measure_z(q) {
                            self.push_gate(Gate::X, &[q]);
                        }
                        self.from_z(*b, q);
                    }
                }
                Op::Noise(n, qs) => {
                    let a = n.arity();
                    for t in qs.chunks(a) {
                        if self.dynamic_depth > 0 {
                            if self.noisy {
                                let u = self.rng.next_f64();
                                let total = n.total();
                                if u < total {
                                    let p = n.choose(u / total);
                                    self.push_pauli(t, p);
                                    if self.log {
                                        self.faults.push(Fault { site: u64::MAX, class: self.class, pauli: p });
                                    }
                                }
                            }
                        } else {
                            if self.next < self.fired.len() && self.fired[self.next] == self.site {
                                self.next += 1;
                                let p = n.choose(self.rng.next_f64());
                                self.push_pauli(t, p);
                                if self.log {
                                    self.faults.push(Fault { site: self.site, class: self.class, pauli: p });
                                }
                            }
                            self.site += 1;
                        }
                    }
                }
                Op::Mark(c) => self.class = *c,
                Op::Repeat(k, body) => {
                    for _ in 0..*k {
                        self.block(body);
                    }
                }
                Op::Table(offs, cases) => {
                    let len = self.records.len();
                    let mut idx = 0usize;
                    for (i, &k) in offs.iter().enumerate() {
                        if self.records[len - k as usize] {
                            idx |= 1 << i;
                        }
                    }
                    self.dynamic_depth += 1;
                    self.block(&cases[idx]);
                    self.dynamic_depth -= 1;
                }
            }
        }
    }
}

/// Backends apply a buffered run of permutation-and-phase gates; the default applies them one
/// by one.
pub trait RunApply: State {
    fn apply_run(&mut self, run: &[RunGate]) {
        for (g, t) in run {
            self.apply(*g, &t[..g.arity()]);
        }
    }
}

/// Runs one shot of `p` on a fresh state of backend `S`.
pub fn run_shot<S: State + RunApply>(p: &Program, sites: &SiteTree, plan: FaultPlan, rng: &mut Xoshiro, log: bool) -> Result<(ShotOutput, S), String> {
    let fired = match plan {
        FaultPlan::None => vec![],
        FaultPlan::Plain => sites.sample_plain(rng),
        FaultPlan::Exactly(k) => sites.sample_exactly(k, rng)?,
    };
    let mut r = Runner {
        state: S::new(p.num_qubits.max(1)),
        rng,
        records: Vec::with_capacity(p.num_measurements as usize),
        fired,
        next: 0,
        site: 0,
        class: 0,
        dynamic_depth: 0,
        noisy: plan != FaultPlan::None,
        log,
        faults: vec![],
        run: Vec::with_capacity(256),
    };
    r.block(&p.body);
    r.flush();
    Ok((ShotOutput { records: r.records, faults: r.faults }, r.state))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dense::Dense;
    use crate::rng::shot_rng;
    use crate::sparse::Sparse;

    fn shot<S: State + RunApply>(text: &str, seed: u64, plan: FaultPlan) -> ShotOutput {
        let p = Program::parse(text).unwrap();
        let t = SiteTree::build(&p);
        run_shot::<S>(&p, &t, plan, &mut shot_rng(seed, 0), true).unwrap().0
    }

    fn both(text: &str, seed: u64, plan: FaultPlan) -> ShotOutput {
        let a = shot::<Dense>(text, seed, plan);
        let b = shot::<Sparse>(text, seed, plan);
        assert_eq!(a.records, b.records, "dense and sparse disagree on {text:?}");
        assert_eq!(a.faults, b.faults);
        a
    }

    #[test]
    fn bell_pairs_agree() {
        for s in 0..50 {
            let o = both("H 0\nCX 0 1\nM 0 1", s, FaultPlan::None);
            assert_eq!(o.records[0], o.records[1]);
        }
    }

    #[test]
    fn repeat_and_reset() {
        assert_eq!(both("REPEAT 3 {\nX 0\n}\nM 0", 0, FaultPlan::None).records, vec![true]);
        assert_eq!(both("X 0\nR 0\nM 0", 0, FaultPlan::None).records, vec![false]);
        for s in 0..20 {
            assert_eq!(both("H 0\nRX 0\nMX 0", s, FaultPlan::None).records, vec![false]);
            assert_eq!(both("H 0\nRY 0\nMY 0", s, FaultPlan::None).records, vec![false]);
        }
        // |−⟩ measured in X gives 1, S|+⟩ = |+i⟩ measured in Y gives 0, S_DAG|+⟩ gives 1.
        assert_eq!(both("X 0\nH 0\nMX 0", 1, FaultPlan::None).records, vec![true]);
        assert_eq!(both("H 0\nS 0\nMY 0", 1, FaultPlan::None).records, vec![false]);
        assert_eq!(both("H 0\nS_DAG 0\nMY 0", 1, FaultPlan::None).records, vec![true]);
    }

    #[test]
    fn teleportation_with_feed_forward() {
        // Teleport T|+⟩ from qubit 0 to qubit 2, undo the T and check X = +1.
        let text = "H 0\nT 0\nH 1\nCX 1 2\nCX 0 1\nH 0\nM 0 1\n\
                    TABLE rec[-1] rec[-2] {I 2} {X 2} {Z 2} {\nX 2\nZ 2\n}\nT_DAG 2\nMX 2";
        for s in 0..100 {
            let o = both(text, s, FaultPlan::None);
            assert!(!o.records[2], "seed {s}");
        }
    }

    #[test]
    fn static_and_dynamic_faults() {
        // A certain-ish error on a static site is logged with its site index and class.
        let o = both("MARK(a)\nX_ERROR(0.999) 0\nM 0", 2, FaultPlan::Plain);
        assert_eq!(o.records, vec![true]);
        assert_eq!(o.faults, vec![Fault { site: 0, class: 1, pauli: 1 }]);
        // Exactly(1) on two sites puts the fault on one of them.
        let o = both("X_ERROR(0.01) 0\nX_ERROR(0.01) 1\nM 0 1", 3, FaultPlan::Exactly(1));
        assert_eq!(o.records.iter().filter(|&&b| b).count(), 1);
        // Dynamic noise fires in Plain mode and is logged with site u64::MAX.
        let o = both("M 0\nTABLE rec[-1] {X_ERROR(0.999) 1} {I 1}\nM 1", 4, FaultPlan::Plain);
        assert_eq!(o.records, vec![false, true]);
        assert_eq!(o.faults[0].site, u64::MAX);
        // …and not in None mode.
        assert_eq!(both("M 0\nTABLE rec[-1] {X_ERROR(0.999) 1} {I 1}\nM 1", 4, FaultPlan::None).records, vec![false, false]);
    }

    #[test]
    fn pauli_channel_2_conditional_frequencies() {
        // Fire PAULI_CHANNEL_2 with every outcome and read which Pauli happened via Bell-basis-free
        // logging: count packed Paulis.
        let probs: Vec<f64> = (1..=15).map(|i| i as f64 / 1000.0).collect();
        let args: Vec<String> = probs.iter().map(|p| p.to_string()).collect();
        let text = format!("PAULI_CHANNEL_2({}) 0 1", args.join(","));
        let p = Program::parse(&text).unwrap();
        let t = SiteTree::build(&p);
        let total: f64 = probs.iter().sum();
        let mut counts = [0u64; 16];
        let shots = 200_000u64;
        for s in 0..shots {
            let (o, _) = run_shot::<Sparse>(&p, &t, FaultPlan::Exactly(1), &mut shot_rng(9, s), true).unwrap();
            counts[o.faults[0].pauli as usize] += 1;
        }
        let n = crate::gates::Noise::P2([0.0; 15]);
        for j in 0..15 {
            let m = shots as f64 * probs[j] / total;
            let got = counts[n.pauli(j) as usize] as f64;
            assert!((got - m).abs() < 5.0 * m.sqrt(), "outcome {j}: {got} vs {m}");
        }
    }

    #[test]
    fn x_error_rate() {
        let p = Program::parse("X_ERROR(0.2) 0\nM 0").unwrap();
        let t = SiteTree::build(&p);
        let shots = 200_000u64;
        let ones = (0..shots)
            .filter(|&s| run_shot::<Dense>(&p, &t, FaultPlan::Plain, &mut shot_rng(10, s), false).unwrap().0.records[0])
            .count() as f64;
        let m = shots as f64 * 0.2;
        assert!((ones - m).abs() < 5.0 * (m * 0.8).sqrt(), "{ones}");
    }
}
