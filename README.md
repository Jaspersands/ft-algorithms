# Fault-tolerant algorithms, run end to end

Shor's algorithm and molecular phase estimation, compiled to Clifford + T + Toffoli and run
**end to end** on a simulated surface-code quantum computer whose every logical error rate was
**measured** by circuit-level simulation. The results are success probabilities as functions of
code distance, physical error rate and magic-state source, with error budgets that say which
operations make the algorithm fail, and the same error model extrapolated to RSA-2048 and FeMoco.

- **Explainer site** with a live demo (Shor for N = 15 running in your browser on WebAssembly):
  `site/index.html` (serve `site/` with any static server).
- **Technical report**: [report/report.html](report/report.html) · [PDF](report/report.pdf).
  Every number in it is computed from the committed data by `tools/report.py`.

## What is in it

| layer | what | where |
|---|---|---|
| physical | rotated surface code, SD6 circuit noise; memory, lattice-surgery merges and CNOTs simulated and decoded (correlated matching) | [stabilizer-qec](https://pypi.org/project/stabilizer-qec/) |
| logical | each operation's Pauli channel, extracted with baselines pushed through the operation; fits in d and p; logical Y from a reference-qubit experiment | `python/ftalgo/calib/` |
| program | Shor with full modular exponentiation (Textbook Cuccaro & Modern Windowed + MBU); STO-3G chemistry from scratch; iterative QPE with Suzuki-4 steps; gridsynth rotations | `python/ftalgo/{shor,arith,arith_modern,phase,qpe,synth}.py`, `chem/` |
| machine | per-qubit clocks, measured idle channel for every gap, decoder reaction time, magic-state factories (cited) | `python/ftalgo/{arch,schedule}.py` |
| run | **ftsim**: Rust state-vector Monte Carlo (sparse and dense), feed-forward, seeded, stratified by the exact number of faults; Python module and WebAssembly | `src/` |

## Results

<!-- results:start -->
**Calibration**: 316 circuit-level experiments, 0.82 billion shots. Whole surgery experiments predicted from calibrated parts: 0.93–1.26× the measured failure rate (median 1.04×). Stim + PyMatching cross-check: 36 of 36 within 2.2σ.

**Shor (Textbook Cuccaro)** (p = 0.1%, cultivated magic states; peak probability = share of runs on one of the r ideal peaks):

| N | logical qubits | Toffolis | noiseless peak | d for 90% of it | physical qubits | run time |
|---:|---:|---:|---:|---:|---:|---:|
| 15 | 18 | 3,360 | 1.000 | 13 | 55.4 k | 225 ms |
| 21 | 21 | 6,250 | 0.790 | 13 | 58.1 k | 418 ms |
| 35 | 24 | 10,440 | 0.790 | 15 | 63.4 k | 808 ms |
| 77 | 27 | 16,170 | 0.774 | 15 | 67 k | 1.25 s |
| 143 | 30 | 23,680 | 0.774 | 15 | 70.8 k | 1.81 s |

**Shor (Modern Windowed + MBU)** (p = 0.1%, cultivated magic states; 2.7–3.0× Toffoli reduction):

| N | logical qubits | Toffolis | noiseless peak | d for 90% of it | physical qubits | run time |
|---:|---:|---:|---:|---:|---:|---:|
| 15 | 23 | 1,248 | 1.000 | 11 | 46.1 k | 107 ms |
| 21 | 27 | 2,390 | 0.793 | 13 | 47 k | 298 ms |
| 35 | 31 | 3,528 | 0.784 | 13 | 50.2 k | 447 ms |

**Phase estimation** (chemical accuracy, 1.6 mHa from FCI; p = 0.1%, cultivated):

| molecule | T gates per run | noiseless | d for 90% of it | run time |
|---|---:|---:|---:|---:|
| H2 | 1.79e7 | 0.92 | 19 | 7.2 min |
| HeH+ | 4.37e7 | 0.97 | 21 | 20.3 min |

**Cryptographic scale (RSA-2048, p = 0.1%, cultivated magic states, target E[faults] ≤ 0.1)**:

| Workload / Implementation | Logical qubits | Toffolis | Distance d | Physical qubits | Quantum run time |
|---|---:|---:|---:|---:|---:|
| Textbook Cuccaro (this work) | 6,150 | 3.44e11 | 35 | 32.5 M | 1.96 years (1.96 yr) |
| Modern Windowed + MBU (this work) | 4,120 | 8.61e10 | 35 | 43.2 M | 298 days (0.81 yr) |
| Gidney 2025 counts (our measured model) | 1,409 | 6.50e9 | 29 | 5.31 M | 4.63 days |
| Gidney 2025 published (assumed model) | 1,409 | 6.50e9 | 25 | 898 k | 4.96 days |

**Headline defensibility (d = 29 vs. d = 25)**: Compiling Gidney's 2025 counts under our measured SD6 model gives d = 29. At d = 25, our circuit-level simulations decoded by correlated matching measure an idle logical error rate of 4.30e-15 per round (1σ: [2.97e-15, 6.25e-15]), which is **4.3× higher** than Gidney's assumed 1.0e-15. Perturbing noise parameters (Λ ± 1σ, prefactor ± 1σ, factory ε_CCZ × 0.1/10×, correlated matching) keeps required distance at d ∈ {29, 31}. Feeding Gidney's assumed error model into our compiler reproduces his published d = 25.

**Biased noise & XZZX surface codes (RSA-2048 Modern, p = 0.1%, target E[faults] ≤ 0.01)**:

| Bias η = p_Z/p_X | XZZX distance (dX × dZ) | Tile qubits | Symmetric distance | Symmetric tile | XZZX physical | Symmetric physical | Qubit savings |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 25 × 21 | 1144 | 25 | 1352 | 19.1 M | 22.5 M | **1.18×** |
| 10 | 31 × 13 | 896 | 31 | 2048 | 14.9 M | 34.1 M | **2.29×** |
| 50 | 33 × 11 | 816 | 33 | 2312 | 13.6 M | 38.5 M | **2.83×** |
| 100 | 33 × 9 | 680 | 33 | 2312 | 11.3 M | 38.5 M | **3.40×** |
| 500 | 33 × 7 | 544 | 33 | 2312 | 9.06 M | 38.5 M | **4.25×** |
| 1000 | 33 × 7 | 544 | 33 | 2312 | 9.06 M | 38.5 M | **4.25×** |
<!-- results:end -->

## How it is checked

- **ftsim** against Qiskit Aer (100 random circuits' amplitudes; exact noisy distributions, plain and
  stratified), dense against sparse, WebAssembly against native.
- **Arithmetic** verified on every input by an independent bit-level simulator; Shor's noiseless output
  against the exact phase-estimation distribution; every synthesized rotation as a matrix.
- **Modern arithmetic** verified on all inputs: Gidney measurement-based carry uncomputation,
  classical feed-forward CZ fixups, and windowed modular exponentiation match exact algebra.
- **Sensitivity analysis**: parameter covariance, 1σ/2σ bands, and multi-factor tornado analysis.
- **Chemistry** against PySCF (integrals, HF and FCI to 1e-8 Ha) and Szabo & Ostlund.
- **Calibration** against Stim + PyMatching; the reference-qubit experiment against plain memories;
  **composition**: whole surgery experiments predicted from their parts vs measured.

```bash
cargo test --lib --features parallel     # engine
.venv/bin/python -m pytest -q tests      # everything else
node tools/wasm-smoke.mjs                # WebAssembly = native
```

## Reproduce

```bash
python3 -m venv .venv && .venv/bin/pip install maturin numpy scipy matplotlib pytest hypothesis \
    stabilizer-qec stim pymatching qiskit qiskit-aer pygridsynth pyscf
tools/build.sh                           # Rust engine → python/ftalgo/_ftsim.abi3.so
.venv/bin/python tools/calibrate.py      # circuit-level experiments (resumable, ~hours)
.venv/bin/python tools/fit_model.py      # data/calibration/model.json
.venv/bin/python tools/compose_check.py && .venv/bin/python tools/xcheck.py
.venv/bin/python tools/run_shor.py && .venv/bin/python tools/run_qpe.py && .venv/bin/python tools/scale.py
.venv/bin/python tools/site_data.py && .venv/bin/python tools/report.py && .venv/bin/python tools/readme_tables.py
cargo build --release --target wasm32-unknown-unknown --lib && cp target/wasm32-unknown-unknown/release/ftsim.wasm site/wasm/
```

Every dataset records its seed and the versions that produced it.

## Assumptions

Magic-state error rates and factory footprints (Litinski 2019; Gidney, Shutty & Jones 2024; Gidney
2025), the H and S operation model of the fast-block layout, 1 µs rounds and a 10 µs decoder reaction
time are cited inputs, not measured here. Algorithm runs at distances beyond the directly measured
range use extrapolated channels and are marked as such. See the report's Limitations.

MIT licensed.
