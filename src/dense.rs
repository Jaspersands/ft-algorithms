//! A dense state vector: 2^n amplitudes, qubit q = bit q of the index.

use crate::complex::C64;
use crate::gates::Gate;
use crate::state::State;

/// The largest register the dense backend accepts (2^30 amplitudes = 16 GiB).
pub const MAX_DENSE_QUBITS: u32 = 30;

#[derive(Clone, Debug)]
pub struct Dense {
    pub n: u32,
    pub amp: Vec<C64>,
}

impl Dense {
    #[inline]
    fn for_pairs(&mut self, q: u32, mut f: impl FnMut(&mut C64, &mut C64)) {
        let m = 1usize << q;
        let len = self.amp.len();
        let mut base = 0;
        while base < len {
            for i in base..base + m {
                let (lo, hi) = self.amp.split_at_mut(i + m);
                f(&mut lo[i], &mut hi[0]);
            }
            base += 2 * m;
        }
    }

    fn matrix(&mut self, q: u32, u: [C64; 4]) {
        self.for_pairs(q, |a, b| {
            let (x, y) = (*a, *b);
            *a = u[0] * x + u[1] * y;
            *b = u[2] * x + u[3] * y;
        });
    }

    fn diag(&mut self, mask: usize, phase: C64) {
        for (i, a) in self.amp.iter_mut().enumerate() {
            if i & mask == mask {
                *a *= phase;
            }
        }
    }

    /// Swaps amplitude pairs (i, i ^ flip) for every i with `ctrl` bits set and bit `t` clear.
    fn controlled_flip(&mut self, ctrl: usize, t: u32) {
        let tm = 1usize << t;
        for i in 0..self.amp.len() {
            if i & ctrl == ctrl && i & tm == 0 {
                self.amp.swap(i, i | tm);
            }
        }
    }
}

impl State for Dense {
    fn new(num_qubits: u32) -> Dense {
        assert!(num_qubits <= MAX_DENSE_QUBITS, "dense backend limited to {MAX_DENSE_QUBITS} qubits");
        let mut amp = vec![C64::ZERO; 1usize << num_qubits];
        amp[0] = C64::ONE;
        Dense { n: num_qubits, amp }
    }

    fn apply(&mut self, g: Gate, t: &[u32]) {
        match g {
            Gate::I => {}
            Gate::X => self.controlled_flip(0, t[0]),
            Gate::Z | Gate::S | Gate::Sdg | Gate::T | Gate::Tdg => {
                self.diag(1 << t[0], g.diag_phase().unwrap());
            }
            Gate::Y => self.for_pairs(t[0], |a, b| {
                let (x, y) = (*a, *b);
                *a = -(y.mul_i());
                *b = x.mul_i();
            }),
            Gate::H | Gate::SqrtX | Gate::SqrtXdg => self.matrix(t[0], g.matrix1()),
            Gate::CX => self.controlled_flip(1 << t[0], t[1]),
            Gate::CCX => self.controlled_flip((1 << t[0]) | (1 << t[1]), t[2]),
            Gate::CZ => self.diag((1 << t[0]) | (1 << t[1]), -C64::ONE),
            Gate::CCZ => self.diag((1 << t[0]) | (1 << t[1]) | (1 << t[2]), -C64::ONE),
            Gate::CY => {
                let c = 1usize << t[0];
                let tm = 1usize << t[1];
                for i in 0..self.amp.len() {
                    if i & c != 0 && i & tm == 0 {
                        let (x, y) = (self.amp[i], self.amp[i | tm]);
                        self.amp[i] = -(y.mul_i());
                        self.amp[i | tm] = x.mul_i();
                    }
                }
            }
            Gate::Swap => {
                let (a, b) = (1usize << t[0], 1usize << t[1]);
                for i in 0..self.amp.len() {
                    if i & a != 0 && i & b == 0 {
                        self.amp.swap(i, i ^ a ^ b);
                    }
                }
            }
        }
    }

    fn prob_one(&self, q: u32) -> f64 {
        let m = 1usize << q;
        self.amp.iter().enumerate().filter(|(i, _)| i & m != 0).map(|(_, a)| a.norm_sqr()).sum()
    }

    fn collapse(&mut self, q: u32, outcome: bool) {
        let m = 1usize << q;
        let mut norm = 0.0;
        for (i, a) in self.amp.iter_mut().enumerate() {
            if (i & m != 0) != outcome {
                *a = C64::ZERO;
            } else {
                norm += a.norm_sqr();
            }
        }
        let s = 1.0 / norm.sqrt();
        for a in self.amp.iter_mut() {
            *a = a.scale(s);
        }
    }

    fn amplitudes(&self) -> Vec<(u128, C64)> {
        self.amp.iter().enumerate().filter(|(_, a)| a.norm_sqr() > 0.0).map(|(i, a)| (i as u128, *a)).collect()
    }

    fn nnz(&self) -> usize {
        self.amp.iter().filter(|a| a.norm_sqr() > 0.0).count()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::state::testutil;

    #[test]
    fn gates_match_reference() {
        testutil::check_backend::<Dense>();
    }

    #[test]
    fn measurement() {
        testutil::check_measurement::<Dense>();
    }
}
