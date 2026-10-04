# Phase 1: `ftsim` logical-level simulator — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A Rust state-vector simulator for logical Clifford+T+Toffoli programs with Pauli-channel
Monte Carlo (plain and stratified by fault count), classical feed-forward, sparse and dense
backends, exposed to Python as `ftalgo._ftsim`, and cross-checked against Qiskit Aer.

**Architecture:** Text program → parsed `Program` (nested blocks: ops, `REPEAT`, `TABLE`) →
a static noise-site tree (weights and hazards per block, REPEAT handled by multiplication) →
per-shot execution on a `State` trait with `Dense` and `Sparse` implementations. Shots run in
parallel with per-shot RNG streams derived from (seed, shot), so output is independent of the
thread count.

**Tech Stack:** Rust 1.90 (edition 2021), pyo3 0.29 (abi3-py39, optional feature `python`),
rayon (optional feature `parallel`), maturin 1.15, Python 3.13 venv at `.venv`, numpy,
qiskit 2.5 + qiskit-aer 0.17 (tests only), pytest + hypothesis.

## Global Constraints

- No dependency beyond pyo3 and rayon in Rust (own RNG, complex type, hasher) so the crate also
  builds for `wasm32-unknown-unknown` without features.
- Determinism: same (program, seed, shots, mode) ⇒ identical records for any thread count.
- Host toolchain is x86_64 (Rosetta); build the Python module with
  `maturin build --profile python --target aarch64-apple-darwin` and install the wheel
  (the `python` profile is release without strip, which avoids the macOS "mis-aligned LINKEDIT"
  load error).
- Errors are `Result<_, String>` inside Rust and `ValueError` in Python; no panics on bad input.
- Commit after each task with the attribution trailer.

## File structure

```
Cargo.toml                 crate ftsim; features python, parallel
pyproject.toml             maturin, module ftalgo._ftsim, python-source = "python"
src/lib.rs                 module wiring, public API
src/complex.rs             C64 (re, im) arithmetic
src/rng.rs                 SplitMix64, Xoshiro256++ ; shot_rng(seed, shot)
src/program.rs             AST (Op, Gate, Noise, Block), text parser, validation, stats
src/gates.rs               gate kinds: 1q matrices, permutation/diagonal classification
src/state.rs               trait State
src/dense.rs               Dense state vector
src/sparse.rs              Sparse (u128 → C64) state with fast hasher
src/noise.rs               static site tree: counts, weights, hazards, P(k), samplers
src/exec.rs                run one shot: walk program, apply gates, inject faults, TABLE
src/batch.rs               many shots (rayon), record packing, fault logs
src/python.rs              pyo3 bindings (feature python)
python/ftalgo/__init__.py  re-export Program; version
tests/ (Python)            test_engine.py, test_aer_crosscheck.py
```

## Program format (reference for every task)

One instruction per line, `#` comments. Targets are qubit indices (non-negative integers).

- Gates: `I X Y Z H S S_DAG T T_DAG SQRT_X SQRT_X_DAG` (1q, any number of targets, applied to
  each), `CX CY CZ SWAP` (pairs), `CCX CCZ` (triples).
- Measurement/reset: `M MX MY` (each target appends one record bit), `R RX`.
- Noise: `PAULI_CHANNEL_1(px,py,pz) q...`, `PAULI_CHANNEL_2(p_IX,p_IY,p_IZ,p_XI,p_XX,p_XY,p_XZ,
  p_YI,p_YX,p_YY,p_YZ,p_ZI,p_ZX,p_ZY,p_ZZ) a b ...` (Stim order: first letter acts on the first
  qubit), `Z_CHANNEL_3(p1..p7) a b c ...` where index m (1..7) applies Z on a if m&1, b if m&2,
  c if m&4. Sugar: `X_ERROR(p)`, `Z_ERROR(p)`, `DEPOLARIZE1(p)`, `DEPOLARIZE2(p)`.
- Blocks: `REPEAT n {` … `}`; `TABLE rec[-a] rec[-b] ... {` … `} {` … `}` … with exactly 2^k
  case blocks, chosen by index Σ bit_i·2^i where bit_i is the record named by the i-th target.
  Every case must contain the same number of measurements. Noise inside a TABLE is "dynamic".
- `MARK(name)` sets the fault class of subsequent noise sites (default class `unmarked`).
  `TICK` is accepted and ignored by the simulator.

---

### Task 1: Crate scaffold, complex numbers, RNG, parser

**Files:** Create `Cargo.toml`, `src/lib.rs`, `src/complex.rs`, `src/rng.rs`, `src/program.rs`,
`src/gates.rs`.

**Interfaces — Produces:**
- `pub struct C64 { pub re: f64, pub im: f64 }` with `Add, Sub, Mul, Neg, conj(), norm_sqr(),
  scale(f64), from_polar(r, θ)`, consts `ZERO, ONE, I`.
- `pub fn shot_rng(seed: u64, shot: u64) -> Xoshiro` ; `Xoshiro::next_u64()`, `next_f64()` in [0,1).
- `pub enum Gate { I,X,Y,Z,H,S,Sdg,T,Tdg,SqrtX,SqrtXdg,CX,CY,CZ,Swap,CCX,CCZ }`, `Gate::arity()`.
- `pub enum Noise { P1([f64;3]), P2([f64;15]), Z3([f64;7]) }` with `total()`.
- `pub enum Op { Gate(Gate, Vec<u32>), Measure(Basis, Vec<u32>), Reset(Basis, Vec<u32>),
  Noise(Noise, Vec<u32>), Mark(u16), Repeat(u64, Block), Table(Vec<u32> /*rec offsets k≥1*/, Vec<Block>) }`
  where `Basis { Z, X, Y }` and `Block = Vec<Op>`.
- `pub struct Program { pub body: Block, pub num_qubits: u32, pub num_measurements: u64,
  pub classes: Vec<String> }`, `Program::parse(&str) -> Result<Program, String>`.
- Validation: arity, distinct targets within one gate application, probabilities in [0,1] with
  total ≤ 1, TABLE case count 2^k and equal measurement counts, rec offsets refer to existing
  records (offset ≤ records so far, counted statically), REPEAT count ≥ 1.

- [ ] Step 1: Rust unit tests in `program.rs`: parse every gate/noise/block form; reject each
  invalid form above with a message naming the line; `num_measurements` counts REPEAT×body and
  one case per TABLE; `M 0 1` then `TABLE rec[-1] rec[-2] {X 2} {Y 2} {Z 2} {H 2}` parses.
- [ ] Step 2: `cargo test` → fails (missing code).
- [ ] Step 3: implement; RNG test: `shot_rng(1,0)` ≠ `shot_rng(1,1)`, mean of 1e5 `next_f64` in
  (0.49, 0.51).
- [ ] Step 4: `cargo test` passes.
- [ ] Step 5: commit `feat(ftsim): program format, parser, RNG`.

### Task 2: Dense backend

**Files:** Create `src/state.rs`, `src/dense.rs`.

**Interfaces — Produces:**
```rust
pub trait State {
    fn new(num_qubits: u32) -> Self where Self: Sized;      // |0…0⟩
    fn apply(&mut self, g: Gate, t: &[u32]);
    fn pauli(&mut self, q: u32, p: u8);                     // p: 1=X 2=Y 3=Z (Stim xz-bits order: I X Y Z)
    fn prob_one(&self, q: u32) -> f64;                       // Z basis
    fn collapse(&mut self, q: u32, outcome: bool);           // project + renormalize
    fn amplitudes(&self) -> Vec<(u128, C64)>;                // non-zero, sorted by index
    fn nnz(&self) -> usize;
}
```
Dense stores `Vec<C64>` of length 2^n, qubit q = bit q of the index (little-endian).
Fast paths: X/CX/CCX/SWAP as index permutations, Z/S/Sdg/T/Tdg/CZ/CCZ diagonal, H/Y/SqrtX/SqrtXdg
general 2×2 on pairs. Matrices: H = [[1,1],[1,−1]]/√2, S = diag(1,i), T = diag(1,e^{iπ/4}),
SqrtX = ½[[1+i,1−i],[1−i,1+i]], CY = control-Y.

- [ ] Step 1: tests: each gate on each computational basis input of up to 3 qubits equals the
  explicit matrix product (build the full 8×8 matrix by Kronecker products in the test);
  measurement probabilities of H|0⟩ = ½; collapse renormalizes.
- [ ] Step 2: run, fail. Step 3: implement. Step 4: pass. Step 5: commit
  `feat(ftsim): dense state vector`.

### Task 3: Sparse backend and dense = sparse

**Files:** Create `src/sparse.rs`.

Sparse stores `HashMap<u128, C64, FxBuild>` (own 64-bit multiply-xor hasher). Permutation gates
remap keys; diagonal gates multiply in place; branching gates build a new map, accumulating,
and drop entries with `|a|² < 1e-26`. Supports up to 128 qubits.

- [ ] Step 1: test: 200 random circuits (5 qubits, 60 gates each from all gates, fixed seeds)
  give the same amplitudes on Dense and Sparse to 1e-12; sparse handles qubit index 100.
- [ ] Steps 2–5 as usual; commit `feat(ftsim): sparse state, dense/sparse agreement`.

### Task 4: Executor — measurement, reset, REPEAT, TABLE

**Files:** Create `src/exec.rs`.

**Interfaces — Produces:**
```rust
pub struct ShotOutput { pub records: Vec<bool>, pub faults: Vec<Fault> }
pub struct Fault { pub site: u64 /* static site index, or u64::MAX for dynamic */, pub class: u16, pub pauli: u32 /* packed 2 bits per target */ }
pub fn run_shot<S: State>(p: &Program, sites: &SiteTree, plan: &FaultPlan, rng: &mut Xoshiro) -> ShotOutput
```
`M` samples with `prob_one`, collapses, appends the bit. `MX`: H, M, H. `MY`: S_DAG, H, M, H, S.
`R`: measure (not recorded), X if 1. `RX`: R then H. TABLE reads records by offset from the
current end. Noise handling arrives in Task 5; here `FaultPlan::none()` means no faults.

- [ ] Step 1: tests: Bell pair `H 0`, `CX 0 1`, `M 0 1` gives equal bits; `REPEAT 3 { X 0 }`, `M 0`
  gives 1; TABLE teleportation of |+⟩ through a Bell pair with X/Z corrections then `MX` on the
  output is always 0; `R` after `X` gives |0⟩; all on both backends.
- [ ] Steps 2–5; commit `feat(ftsim): executor with feed-forward`.

### Task 5: Noise sites, sampling modes, exact fault-count distribution, fault logs

**Files:** Create `src/noise.rs`; modify `src/exec.rs`.

Static sites are noise instructions outside any TABLE, one site per target group (P1: per qubit,
P2: per pair, Z3: per triple), in program order with REPEAT bodies repeated. For each block the
tree stores, per item, cumulative site count, cumulative odds weight w = q/(1−q) and cumulative
hazard h = −ln(1−q); a REPEAT item has count n·c_body, weight n·W_body, hazard n·H_body.

**Interfaces — Produces:**
```rust
pub struct SiteTree { /* … */ }
impl SiteTree {
    pub fn build(p: &Program) -> SiteTree;
    pub fn num_sites(&self) -> u64;
    pub fn fault_count_distribution(&self, kmax: usize) -> Vec<f64>;   // exact P(K=k), k≤kmax, plus P(K>kmax) last
    pub fn sample_plain(&self, rng) -> Vec<u64>;                          // sorted fired sites
    pub fn sample_exactly(&self, k: usize, rng) -> Vec<u64>;              // sorted, distinct
}
pub enum FaultPlan { Plain, Exactly(usize), None }
```
- P(k): truncated generating polynomial Π(1−q + q z); REPEAT by exponentiation by squaring,
  truncated at kmax+1 (the tail term = 1 − Σ).
- `sample_plain`: walk with exponential(1) hazard budgets: next fired site after g is the first
  index where cumulative hazard exceeds H(g) + E, found by descending the tree. Exact for
  independent Bernoulli sites.
- `sample_exactly(k)`: draw k sites i.i.d. ∝ w (descend by weight), reject and redraw the whole
  set if any duplicate. This samples sets with probability ∝ Π w_i, exactly the conditional law
  given K = k.
- Executor: a counter of static sites; when it equals the next fired site, choose the Pauli from
  the channel's conditional distribution (probabilities / total) and apply it, logging a Fault.
  Dynamic sites (inside TABLE cases) fire by a Bernoulli draw each time they execute in both
  Plain and Exactly modes (site = u64::MAX in logs); in mode None, nothing fires.

- [ ] Step 1: tests:
  - `num_sites` of `REPEAT 4 { DEPOLARIZE1(0.1) 0 1 }` is 8.
  - P(k) for 8 sites at q=0.1 equals Binomial(8, 0.1) to 1e-12; with mixed q equals a direct DP.
  - Plain mode on `X_ERROR(0.2) 0`, `M 0`: mean of records over 200k shots within 4σ of 0.2.
  - Plain on 1000 sites with q=1e-3 under REPEAT: mean fired count within 4σ of 1.0; each site's
    empirical frequency uniform (chi-square p > 1e-4).
  - Exactly(2) on mixed weights: pair frequencies ∝ w_i w_j (chi-square).
  - PAULI_CHANNEL_2 conditional choice frequencies match the given vector (chi-square).
- [ ] Steps 2–5; commit `feat(ftsim): noise sites, stratified sampling, exact P(k)`.

### Task 6: Batches, determinism, parallel shots

**Files:** Create `src/batch.rs`.

```rust
pub enum Backend { Auto, Dense, Sparse }   // Auto: Dense if num_qubits ≤ 20, else Sparse
pub struct BatchResult { pub records: Vec<u8> /* shots × ceil(m/8), little-endian bits */,
                         pub fault_shot: Vec<u64>, pub fault_site: Vec<u64>, pub fault_class: Vec<u16>, pub fault_pauli: Vec<u32> }
pub fn sample(p: &Program, shots: u64, seed: u64, plan: FaultPlan, backend: Backend, threads: usize, log_faults: bool) -> Result<BatchResult, String>
```
Each shot uses `shot_rng(seed, shot)`. With feature `parallel`, shots are split into chunks
processed by a rayon pool of `threads` (0 = all cores); results concatenated in shot order.

- [ ] Step 1: tests: records identical for threads 1, 2 and 7 on a noisy random program;
  identical between Dense and Sparse backends for the same seed (same RNG consumption order is
  part of the contract: one f64 per measurement, faults drawn before the shot executes, dynamic
  draws in execution order).
- [ ] Steps 2–5; commit `feat(ftsim): seeded parallel batches`.

### Task 7: Python bindings and package

**Files:** Create `src/python.rs`, `pyproject.toml`, `python/ftalgo/__init__.py`,
`tests/test_engine.py`, `tools/build.sh`.

Python API (module `ftalgo._ftsim`, re-exported as `ftalgo.Program`):
```python
p = ftalgo.Program(text)
p.num_qubits; p.num_measurements; p.num_sites; p.classes  # list[str]
p.fault_count_distribution(kmax) -> list[float]
p.sample(shots, *, seed, faults=None|"plain"|int, backend="auto", threads=0, log_faults=False)
    -> numpy bool array (shots, num_measurements)              if not log_faults
    -> (records, dict(shot=, site=, cls=, pauli=) numpy arrays) if log_faults
p.statevector(*, seed=0, backend="dense") -> (np.ndarray complex128 of length 2**n, records)  # one noiseless shot
```
`tools/build.sh` runs `maturin build --profile python --target aarch64-apple-darwin --out dist`
and `pip install --force-reinstall --no-deps dist/ftalgo-*.whl` into `.venv`.

- [ ] Step 1: `tests/test_engine.py`: Bell sampling; seed determinism across threads; errors
  raise ValueError; statevector of `H 0`, `T 0` equals [1, e^{iπ/4}]/√2.
- [ ] Steps 2–5; commit `feat: Python bindings ftalgo._ftsim`.

### Task 8: Cross-checks against Qiskit Aer

**Files:** Create `python/ftalgo/qiskit_bridge.py`, `tests/test_aer_crosscheck.py`.

`to_qiskit(text) -> QuantumCircuit` for gate-only programs (noise ignored) and
`exact_distribution(text) -> dict[bitstring, prob]` for programs of gates, Pauli channels and
final `M`, using Aer's density-matrix method with each channel as a `pauli_error` instruction.

- [ ] Step 1: tests:
  - 100 random 6-qubit, 80-gate circuits: our statevector equals `qiskit.quantum_info.Statevector`
    to 1e-10 (qubit order: both little-endian).
  - 10 random noisy 4-qubit programs (P1, P2, Z3 channels, q up to 0.05): our 200k-shot
    distribution vs Aer's exact one, chi-square p-value > 1e-4 each; and the fault-stratified
    estimate Σ_k P(k)·P̂(outcome | k) agrees too.
- [ ] Steps 2–5; commit `test: cross-check against Qiskit Aer`.
