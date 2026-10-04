//! `ftsim`: a logical-level simulator for fault-tolerant quantum programs.
//!
//! Programs are Clifford + T + Toffoli circuits in a Stim-like text format, with Pauli noise
//! channels, measurements, `REPEAT` blocks and classical feed-forward (`TABLE`). Shots run on a
//! sparse or dense state vector; noise is sampled independently per site or stratified by the
//! exact number of faults.

pub mod complex;
pub mod gates;
pub mod program;
pub mod rng;

pub use complex::C64;
pub use gates::{Basis, Gate, Noise};
pub use program::{Block, Op, Program};
