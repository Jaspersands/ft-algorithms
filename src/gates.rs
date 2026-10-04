//! The gate set and the noise channels of the program format.

use crate::complex::C64;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Gate {
    I,
    X,
    Y,
    Z,
    H,
    S,
    Sdg,
    T,
    Tdg,
    SqrtX,
    SqrtXdg,
    CX,
    CY,
    CZ,
    Swap,
    CCX,
    CCZ,
}

impl Gate {
    pub fn from_name(name: &str) -> Option<Gate> {
        Some(match name {
            "I" => Gate::I,
            "X" => Gate::X,
            "Y" => Gate::Y,
            "Z" => Gate::Z,
            "H" => Gate::H,
            "S" => Gate::S,
            "S_DAG" => Gate::Sdg,
            "T" => Gate::T,
            "T_DAG" => Gate::Tdg,
            "SQRT_X" => Gate::SqrtX,
            "SQRT_X_DAG" => Gate::SqrtXdg,
            "CX" | "CNOT" => Gate::CX,
            "CY" => Gate::CY,
            "CZ" => Gate::CZ,
            "SWAP" => Gate::Swap,
            "CCX" | "TOFFOLI" => Gate::CCX,
            "CCZ" => Gate::CCZ,
            _ => return None,
        })
    }

    pub fn name(self) -> &'static str {
        match self {
            Gate::I => "I",
            Gate::X => "X",
            Gate::Y => "Y",
            Gate::Z => "Z",
            Gate::H => "H",
            Gate::S => "S",
            Gate::Sdg => "S_DAG",
            Gate::T => "T",
            Gate::Tdg => "T_DAG",
            Gate::SqrtX => "SQRT_X",
            Gate::SqrtXdg => "SQRT_X_DAG",
            Gate::CX => "CX",
            Gate::CY => "CY",
            Gate::CZ => "CZ",
            Gate::Swap => "SWAP",
            Gate::CCX => "CCX",
            Gate::CCZ => "CCZ",
        }
    }

    pub fn arity(self) -> usize {
        match self {
            Gate::CX | Gate::CY | Gate::CZ | Gate::Swap => 2,
            Gate::CCX | Gate::CCZ => 3,
            _ => 1,
        }
    }

    /// The 2×2 matrix of a single-qubit gate, row-major [m00, m01, m10, m11].
    pub fn matrix1(self) -> [C64; 4] {
        let h = std::f64::consts::FRAC_1_SQRT_2;
        let z = C64::ZERO;
        let o = C64::ONE;
        let i = C64::I;
        match self {
            Gate::I => [o, z, z, o],
            Gate::X => [z, o, o, z],
            Gate::Y => [z, -i, i, z],
            Gate::Z => [o, z, z, -o],
            Gate::H => [C64::new(h, 0.0), C64::new(h, 0.0), C64::new(h, 0.0), C64::new(-h, 0.0)],
            Gate::S => [o, z, z, i],
            Gate::Sdg => [o, z, z, -i],
            Gate::T => [o, z, z, C64::new(h, h)],
            Gate::Tdg => [o, z, z, C64::new(h, -h)],
            Gate::SqrtX => [C64::new(0.5, 0.5), C64::new(0.5, -0.5), C64::new(0.5, -0.5), C64::new(0.5, 0.5)],
            Gate::SqrtXdg => [C64::new(0.5, -0.5), C64::new(0.5, 0.5), C64::new(0.5, 0.5), C64::new(0.5, -0.5)],
            _ => panic!("matrix1 of a multi-qubit gate"),
        }
    }

    /// The phase on |1⟩ of a diagonal single-qubit gate (|0⟩ keeps phase 1).
    pub fn diag_phase(self) -> Option<C64> {
        let h = std::f64::consts::FRAC_1_SQRT_2;
        match self {
            Gate::Z => Some(-C64::ONE),
            Gate::S => Some(C64::I),
            Gate::Sdg => Some(-C64::I),
            Gate::T => Some(C64::new(h, h)),
            Gate::Tdg => Some(C64::new(h, -h)),
            _ => None,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Basis {
    Z,
    X,
    Y,
}

/// A Pauli channel acting on one target group. Paulis are packed two bits per qubit
/// (0 = I, 1 = X, 2 = Y, 3 = Z), the group's first qubit in the lowest bits.
#[derive(Clone, Debug, PartialEq)]
pub enum Noise {
    /// [p_X, p_Y, p_Z] on one qubit.
    P1([f64; 3]),
    /// Stim's PAULI_CHANNEL_2 order: index j = 4·a + b − 1 for Paulis a (first qubit), b (second).
    P2([f64; 15]),
    /// Z-type errors on three qubits: entry m − 1 applies Z on qubit k of the group when m has bit k.
    Z3([f64; 7]),
}

impl Noise {
    pub fn arity(&self) -> usize {
        match self {
            Noise::P1(_) => 1,
            Noise::P2(_) => 2,
            Noise::Z3(_) => 3,
        }
    }

    pub fn probs(&self) -> &[f64] {
        match self {
            Noise::P1(p) => p,
            Noise::P2(p) => p,
            Noise::Z3(p) => p,
        }
    }

    pub fn total(&self) -> f64 {
        self.probs().iter().sum()
    }

    /// The packed Pauli for outcome index `j` (0-based, into `probs()`).
    pub fn pauli(&self, j: usize) -> u32 {
        match self {
            Noise::P1(_) => (j + 1) as u32,
            Noise::P2(_) => {
                let k = j + 1;
                let a = (k / 4) as u32;
                let b = (k % 4) as u32;
                a | (b << 2)
            }
            Noise::Z3(_) => {
                let m = j + 1;
                let mut p = 0u32;
                for q in 0..3 {
                    if m >> q & 1 == 1 {
                        p |= 3 << (2 * q);
                    }
                }
                p
            }
        }
    }

    /// Chooses the Pauli of a fired site from the conditional distribution probs / total.
    pub fn choose(&self, u: f64) -> u32 {
        let probs = self.probs();
        let total = self.total();
        let mut acc = 0.0;
        let target = u * total;
        let mut last = 0;
        for (j, &p) in probs.iter().enumerate() {
            if p > 0.0 {
                last = j;
                acc += p;
                if target < acc {
                    return self.pauli(j);
                }
            }
        }
        self.pauli(last)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn pauli_packing() {
        assert_eq!(Noise::P1([0.0; 3]).pauli(0), 1);
        assert_eq!(Noise::P1([0.0; 3]).pauli(2), 3);
        // IX: first qubit I, second X.
        assert_eq!(Noise::P2([0.0; 15]).pauli(0), 1 << 2);
        // XI.
        assert_eq!(Noise::P2([0.0; 15]).pauli(3), 1);
        // ZZ.
        assert_eq!(Noise::P2([0.0; 15]).pauli(14), 3 | 3 << 2);
        // Z on qubits 0 and 2.
        assert_eq!(Noise::Z3([0.0; 7]).pauli(4), 3 | 3 << 4);
    }
    #[test]
    fn choose_respects_zeros() {
        let n = Noise::P1([0.0, 0.0, 0.3]);
        for u in [0.0, 0.5, 0.999999] {
            assert_eq!(n.choose(u), 3);
        }
    }
}
