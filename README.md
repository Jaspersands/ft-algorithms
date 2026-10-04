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
| program | Shor with the full modular exponentiation (Toffoli arithmetic, semiclassical AQFT); STO-3G chemistry from scratch; iterative QPE with Suzuki-4 steps; gridsynth rotations | `python/ftalgo/{shor,arith,phase,qpe,synth}.py`, `chem/` |
| machine | per-qubit clocks, measured idle channel for every gap, decoder reaction time, magic-state factories (cited) | `python/ftalgo/{arch,schedule}.py` |
| run | **ftsim**: Rust state-vector Monte Carlo (sparse and dense), feed-forward, seeded, stratified by the exact number of faults; Python module and WebAssembly | `src/` |

## Results

<!-- results:start -->
**Calibration**: 421 circuit-level experiments, 1.21 billion shots. Whole surgery experiments predicted from calibrated parts: 0.93–1.26× the measured failure rate (median 1.04×). Stim + PyMatching cross-check: 36 of 36 within 2.2σ.

**Shor** (p = 0.1%, cultivated magic states; peak probability = share of runs on one of the r ideal peaks):

| N | logical qubits | Toffolis | noiseless peak | d for 90% of it | physical qubits | run time |
|---:|---:|---:|---:|---:|---:|---:|
| 15 | 18 | 3,360 | 1.000 | 13 | 55.4 k | 225 ms |
| 21 | 21 | 6,250 | 0.790 | 13 | 58.1 k | 418 ms |

**To scale** (≤ 0.1 expected faults, p = 0.1%): textbook-arithmetic RSA-2048 needs d = 35, 32.5 M physical qubits and 1.96 years; Gidney's 2025 counts need d = 29 under our measured error model (he assumed d = 25).
<!-- results:end -->

## How it is checked

- **ftsim** against Qiskit Aer (100 random circuits' amplitudes; exact noisy distributions, plain and
  stratified), dense against sparse, WebAssembly against native.
- **Arithmetic** verified on every input by an independent bit-level simulator; Shor's noiseless output
  against the exact phase-estimation distribution; every synthesized rotation as a matrix.
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
