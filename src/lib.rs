//! `ftsim`: a logical-level simulator for fault-tolerant quantum programs.
//!
//! Programs are Clifford + T + Toffoli circuits in a Stim-like text format, with Pauli noise
//! channels, measurements, `REPEAT` blocks and classical feed-forward (`TABLE`). Shots run on a
//! sparse or dense state vector; noise is sampled independently per site or stratified by the
//! exact number of faults.

pub mod batch;
pub mod complex;
pub mod dense;
pub mod exec;
pub mod gates;
pub mod noise;
pub mod program;
#[cfg(feature = "python")]
pub mod python;
pub mod rng;
pub mod sparse;
pub mod state;
#[cfg(target_arch = "wasm32")]
pub mod wasm;

pub use complex::C64;
pub use gates::{Basis, Gate, Noise};
pub use dense::Dense;
pub use program::{Block, Op, Program};
pub use sparse::Sparse;
pub use state::State;
