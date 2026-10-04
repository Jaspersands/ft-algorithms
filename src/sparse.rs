//! A sparse state: only non-zero amplitudes, keyed by basis index. Permutation and diagonal
//! gates (X, CX, CCX, SWAP, Z, S, T, CZ, CCZ, Y, CY) never grow it; H and √X can at most
//! double it, and amplitudes that cancel are dropped.

use crate::complex::C64;
use crate::gates::Gate;
use crate::state::State;
use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};

/// A fast multiply-xorshift hasher for u128 keys (no HashDoS concern: keys are basis states).
#[derive(Default, Clone, Copy)]
pub struct KeyHasher(u64);

impl Hasher for KeyHasher {
    #[inline]
    fn finish(&self) -> u64 {
        let mut z = self.0;
        z = (z ^ (z >> 32)).wrapping_mul(0xD6E8_FEB8_6659_FD93);
        z ^ (z >> 32)
    }
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 = (self.0.rotate_left(8) ^ b as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
        }
    }
    #[inline]
    fn write_u128(&mut self, k: u128) {
        let lo = k as u64;
        let hi = (k >> 64) as u64;
        self.0 = (lo ^ hi.rotate_left(29)).wrapping_mul(0x9E37_79B9_7F4A_7C15);
    }
}

pub type Map = HashMap<u128, C64, BuildHasherDefault<KeyHasher>>;

/// Amplitudes with |a|² below this are treated as exact cancellations and dropped.
pub const DROP: f64 = 1e-26;

#[derive(Clone, Debug)]
pub struct Sparse {
    pub n: u32,
    pub amp: Map,
    scratch: Map,
}

impl Sparse {
    /// Rebuilds the map through a basis permutation with phases.
    #[inline]
    fn remap(&mut self, f: impl Fn(u128) -> (u128, C64)) {
        self.scratch.clear();
        self.scratch.reserve(self.amp.len());
        for (&k, &a) in self.amp.iter() {
            let (k2, ph) = f(k);
            self.scratch.insert(k2, a * ph);
        }
        std::mem::swap(&mut self.amp, &mut self.scratch);
    }

    #[inline]
    fn diag(&mut self, mask: u128, phase: C64) {
        for (k, a) in self.amp.iter_mut() {
            if k & mask == mask {
                *a *= phase;
            }
        }
    }

    fn branch(&mut self, q: u32, u: [C64; 4]) {
        let m = 1u128 << q;
        self.scratch.clear();
        self.scratch.reserve(self.amp.len() * 2);
        for (&k, &a) in self.amp.iter() {
            let b = ((k >> q) & 1) as usize;
            let k0 = k & !m;
            let k1 = k | m;
            let c0 = u[b] * a;
            let c1 = u[2 + b] * a;
            if c0.norm_sqr() > 0.0 {
                *self.scratch.entry(k0).or_insert(C64::ZERO) += c0;
            }
            if c1.norm_sqr() > 0.0 {
                *self.scratch.entry(k1).or_insert(C64::ZERO) += c1;
            }
        }
        self.scratch.retain(|_, a| a.norm_sqr() > DROP);
        std::mem::swap(&mut self.amp, &mut self.scratch);
    }
}

impl State for Sparse {
    fn new(num_qubits: u32) -> Sparse {
        assert!(num_qubits <= 128, "sparse backend limited to 128 qubits");
        let mut amp = Map::default();
        amp.insert(0, C64::ONE);
        Sparse { n: num_qubits, amp, scratch: Map::default() }
    }

    fn apply(&mut self, g: Gate, t: &[u32]) {
        match g {
            Gate::I => {}
            Gate::X => {
                let m = 1u128 << t[0];
                self.remap(|k| (k ^ m, C64::ONE));
            }
            Gate::Y => {
                let q = t[0];
                let m = 1u128 << q;
                // Y|0⟩ = i|1⟩, Y|1⟩ = −i|0⟩.
                self.remap(|k| (k ^ m, if k & m == 0 { C64::I } else { -C64::I }));
            }
            Gate::Z | Gate::S | Gate::Sdg | Gate::T | Gate::Tdg => self.diag(1 << t[0], g.diag_phase().unwrap()),
            Gate::H | Gate::SqrtX | Gate::SqrtXdg => self.branch(t[0], g.matrix1()),
            Gate::CX => {
                let (c, m) = (1u128 << t[0], 1u128 << t[1]);
                self.remap(|k| (if k & c != 0 { k ^ m } else { k }, C64::ONE));
            }
            Gate::CY => {
                let (c, m) = (1u128 << t[0], 1u128 << t[1]);
                self.remap(|k| {
                    if k & c != 0 {
                        (k ^ m, if k & m == 0 { C64::I } else { -C64::I })
                    } else {
                        (k, C64::ONE)
                    }
                });
            }
            Gate::CCX => {
                let c = (1u128 << t[0]) | (1u128 << t[1]);
                let m = 1u128 << t[2];
                self.remap(|k| (if k & c == c { k ^ m } else { k }, C64::ONE));
            }
            Gate::Swap => {
                let (a, b) = (1u128 << t[0], 1u128 << t[1]);
                self.remap(|k| (if ((k & a != 0) as u8) ^ ((k & b != 0) as u8) == 1 { k ^ a ^ b } else { k }, C64::ONE));
            }
            Gate::CZ => self.diag((1 << t[0]) | (1 << t[1]), -C64::ONE),
            Gate::CCZ => self.diag((1 << t[0]) | (1 << t[1]) | (1 << t[2]), -C64::ONE),
        }
    }

    fn prob_one(&self, q: u32) -> f64 {
        let m = 1u128 << q;
        self.amp.iter().filter(|(k, _)| *k & m != 0).map(|(_, a)| a.norm_sqr()).sum()
    }

    fn collapse(&mut self, q: u32, outcome: bool) {
        let m = 1u128 << q;
        self.amp.retain(|k, _| (k & m != 0) == outcome);
        let norm: f64 = self.amp.values().map(|a| a.norm_sqr()).sum();
        let s = 1.0 / norm.sqrt();
        for a in self.amp.values_mut() {
            *a = a.scale(s);
        }
    }

    fn amplitudes(&self) -> Vec<(u128, C64)> {
        let mut v: Vec<(u128, C64)> = self.amp.iter().map(|(k, a)| (*k, *a)).collect();
        v.sort_by_key(|x| x.0);
        v
    }

    fn nnz(&self) -> usize {
        self.amp.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dense::Dense;
    use crate::rng::Xoshiro;
    use crate::state::testutil;

    #[test]
    fn gates_match_reference() {
        testutil::check_backend::<Sparse>();
    }

    #[test]
    fn measurement() {
        testutil::check_measurement::<Sparse>();
    }

    #[test]
    fn agrees_with_dense_on_random_circuits() {
        let n = 5;
        for seed in 0..200 {
            let mut rng = Xoshiro::from_seed(seed);
            let mut d = Dense::new(n);
            let mut s = Sparse::new(n);
            for _ in 0..60 {
                let (g, t) = testutil::random_gate(&mut rng, n);
                d.apply(g, &t);
                s.apply(g, &t);
            }
            testutil::assert_close(&testutil::dense_vec(&s, n), &d.amp, 1e-12);
        }
    }

    #[test]
    fn high_qubits_and_cancellation() {
        let mut s = Sparse::new(101);
        s.apply(Gate::H, &[100]);
        s.apply(Gate::CX, &[100, 3]);
        assert_eq!(s.nnz(), 2);
        s.apply(Gate::CX, &[100, 3]);
        s.apply(Gate::H, &[100]);
        assert_eq!(s.nnz(), 1, "H·H must cancel back to one amplitude");
        assert_eq!(s.amplitudes()[0].0, 0);
    }
}
