//! The state-vector interface both backends implement. Qubit q is bit q of a basis index.

use crate::complex::C64;
use crate::gates::Gate;

pub trait State {
    /// |0…0⟩ on `num_qubits` qubits.
    fn new(num_qubits: u32) -> Self
    where
        Self: Sized;
    /// Applies `g` to one target group (`t.len() == g.arity()`).
    fn apply(&mut self, g: Gate, t: &[u32]);
    /// Applies the Pauli `p` (1 = X, 2 = Y, 3 = Z; 0 does nothing) to qubit `q`.
    fn pauli(&mut self, q: u32, p: u32) {
        match p & 3 {
            1 => self.apply(Gate::X, &[q]),
            2 => self.apply(Gate::Y, &[q]),
            3 => self.apply(Gate::Z, &[q]),
            _ => {}
        }
    }
    /// Probability that a Z measurement of `q` gives 1.
    fn prob_one(&self, q: u32) -> f64;
    /// Projects `q` onto `outcome` and renormalizes. The outcome must have non-zero probability.
    fn collapse(&mut self, q: u32, outcome: bool);
    /// Non-zero amplitudes sorted by basis index.
    fn amplitudes(&self) -> Vec<(u128, C64)>;
    /// Number of stored non-zero amplitudes.
    fn nnz(&self) -> usize;
}

#[cfg(test)]
pub(crate) mod testutil {
    //! A slow, obviously correct reference: full unitaries built by Kronecker products.
    use super::*;
    use crate::rng::Xoshiro;

    pub type Mat = Vec<Vec<C64>>;

    /// The full 2^n × 2^n unitary of `g` on targets `t`, by its action on basis states
    /// written out from the gate definitions.
    pub fn full_unitary(n: u32, g: Gate, t: &[u32]) -> Mat {
        let dim = 1usize << n;
        let mut u = vec![vec![C64::ZERO; dim]; dim];
        for col in 0..dim {
            // Column `col` = image of basis state |col⟩.
            let bit = |q: u32| (col >> q) & 1;
            match g.arity() {
                1 => {
                    let m = g.matrix1();
                    let q = t[0];
                    let b = bit(q);
                    let c0 = col & !(1 << q);
                    let c1 = col | (1 << q);
                    u[c0][col] += m[b];
                    u[c1][col] += m[2 + b];
                }
                _ => {
                    let (row, phase) = match g {
                        Gate::CX => (if bit(t[0]) == 1 { col ^ (1 << t[1]) } else { col }, C64::ONE),
                        Gate::CZ => (col, if bit(t[0]) & bit(t[1]) == 1 { -C64::ONE } else { C64::ONE }),
                        Gate::CY => {
                            if bit(t[0]) == 1 {
                                let ph = if bit(t[1]) == 0 { C64::I } else { -C64::I };
                                (col ^ (1 << t[1]), ph)
                            } else {
                                (col, C64::ONE)
                            }
                        }
                        Gate::Swap => {
                            let (a, b) = (bit(t[0]), bit(t[1]));
                            let mut r = col & !(1 << t[0]) & !(1 << t[1]);
                            r |= b << t[0];
                            r |= a << t[1];
                            (r, C64::ONE)
                        }
                        Gate::CCX => (if bit(t[0]) & bit(t[1]) == 1 { col ^ (1 << t[2]) } else { col }, C64::ONE),
                        Gate::CCZ => (col, if bit(t[0]) & bit(t[1]) & bit(t[2]) == 1 { -C64::ONE } else { C64::ONE }),
                        _ => unreachable!(),
                    };
                    u[row][col] += phase;
                }
            }
        }
        u
    }

    pub fn mat_vec(m: &Mat, v: &[C64]) -> Vec<C64> {
        m.iter().map(|row| row.iter().zip(v).fold(C64::ZERO, |acc, (a, b)| acc + *a * *b)).collect()
    }

    pub const ALL_GATES: [Gate; 17] = [
        Gate::I, Gate::X, Gate::Y, Gate::Z, Gate::H, Gate::S, Gate::Sdg, Gate::T, Gate::Tdg,
        Gate::SqrtX, Gate::SqrtXdg, Gate::CX, Gate::CY, Gate::CZ, Gate::Swap, Gate::CCX, Gate::CCZ,
    ];

    /// A random gate on distinct random qubits of an n-qubit register.
    pub fn random_gate(rng: &mut Xoshiro, n: u32) -> (Gate, Vec<u32>) {
        let g = ALL_GATES[(rng.next_u64() % ALL_GATES.len() as u64) as usize];
        let mut t: Vec<u32> = Vec::new();
        while t.len() < g.arity() {
            let q = (rng.next_u64() % n as u64) as u32;
            if !t.contains(&q) {
                t.push(q);
            }
        }
        (g, t)
    }

    pub fn dense_vec<S: State>(s: &S, n: u32) -> Vec<C64> {
        let mut v = vec![C64::ZERO; 1 << n];
        for (k, a) in s.amplitudes() {
            v[k as usize] = a;
        }
        v
    }

    pub fn assert_close(a: &[C64], b: &[C64], tol: f64) {
        assert_eq!(a.len(), b.len());
        for (i, (x, y)) in a.iter().zip(b).enumerate() {
            assert!((*x - *y).norm_sqr().sqrt() < tol, "index {i}: {x:?} vs {y:?}");
        }
    }

    /// Every gate, on every target choice of a 3-qubit register, from every basis state and
    /// from a random superposition, matches the reference unitary.
    pub fn check_backend<S: State>() {
        let n = 3u32;
        let mut rng = Xoshiro::from_seed(11);
        for g in ALL_GATES {
            let mut targets: Vec<Vec<u32>> = Vec::new();
            for a in 0..n {
                for b in 0..n {
                    for c in 0..n {
                        let t: Vec<u32> = [a, b, c][..g.arity()].to_vec();
                        let mut uniq = t.clone();
                        uniq.sort();
                        uniq.dedup();
                        if uniq.len() == t.len() && !targets.contains(&t) {
                            targets.push(t);
                        }
                    }
                }
            }
            for t in targets {
                let u = full_unitary(n, g, &t);
                // A random state prepared by a fixed gate sequence on both sides.
                let mut s = S::new(n);
                let mut prep = Vec::new();
                for _ in 0..12 {
                    let (pg, pt) = random_gate(&mut rng, n);
                    s.apply(pg, &pt);
                    prep.push((pg, pt));
                }
                let mut v = vec![C64::ZERO; 1 << n];
                v[0] = C64::ONE;
                for (pg, pt) in &prep {
                    v = mat_vec(&full_unitary(n, *pg, pt), &v);
                }
                assert_close(&dense_vec(&s, n), &v, 1e-12);
                s.apply(g, &t);
                let want = mat_vec(&u, &v);
                assert_close(&dense_vec(&s, n), &want, 1e-12);
            }
        }
    }

    pub fn check_measurement<S: State>() {
        let mut s = S::new(2);
        s.apply(Gate::H, &[0]);
        s.apply(Gate::CX, &[0, 1]);
        assert!((s.prob_one(0) - 0.5).abs() < 1e-12);
        assert!((s.prob_one(1) - 0.5).abs() < 1e-12);
        s.collapse(0, true);
        assert!((s.prob_one(1) - 1.0).abs() < 1e-12);
        let amps = s.amplitudes();
        assert_eq!(amps.len(), 1);
        assert_eq!(amps[0].0, 3);
        assert!((amps[0].1.norm_sqr() - 1.0).abs() < 1e-12);
        let mut t = S::new(1);
        t.apply(Gate::H, &[0]);
        t.apply(Gate::T, &[0]);
        t.collapse(0, false);
        assert!((t.amplitudes()[0].1.norm_sqr() - 1.0).abs() < 1e-12);
    }
}
