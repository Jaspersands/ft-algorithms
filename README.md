# Fault-tolerant algorithms, run end to end

[![CI](https://github.com/Jaspersands/ft-algorithms/actions/workflows/ci.yml/badge.svg)](https://github.com/Jaspersands/ft-algorithms/actions/workflows/ci.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23792569.svg)](https://doi.org/10.5281/zenodo.23792569)

Shor's algorithm and molecular phase estimation, compiled to Clifford + T + Toffoli and run
**end to end** on a simulated surface-code quantum computer whose every logical error rate was
**measured** by circuit-level simulation. The results are success probabilities as functions of
code distance, physical error rate and magic-state source, with error budgets that say which
operations make the algorithm fail, and the same error model extrapolated to RSA-2048 and FeMoco.

- **Explainer site** with a live demo (Shor for N = 15 and H₂ phase estimation running in your
  browser on WebAssembly): [ft.jaspersands.com](https://ft.jaspersands.com) (source in `site/`).
- **Technical report**: [report/report.html](report/report.html) · [PDF](report/report.pdf).
  Every number in it is computed from the committed data by `tools/report.py`.

## What is in it

| layer | what | where |
|---|---|---|
| physical | rotated surface code, SD6 circuit noise; memory, lattice-surgery merges and CNOTs simulated and decoded (correlated matching) | [stabilizer-qec](https://pypi.org/project/stabilizer-qec/) |
| logical | each operation's Pauli channel, extracted with baselines pushed through the operation; fits in d and p; logical Y from a reference-qubit experiment | `python/ftalgo/calib/` |
| program | Shor with full modular exponentiation (textbook Cuccaro, and windowed with measurement-based uncomputation); STO-3G chemistry from scratch for H₂ and HeH⁺, a PySCF CASSCF active space for LiH; iterative QPE with Suzuki-4 steps; gridsynth rotations | `python/ftalgo/{shor,arith,arith_modern,phase,qpe,synth}.py`, `chem/` |
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

**Shor (Modern Windowed + MBU)** (p = 0.1%, cultivated magic states; 2.6–3.0× fewer Toffolis):

| N | logical qubits | Toffolis | noiseless peak | d for 90% of it | physical qubits | run time |
|---:|---:|---:|---:|---:|---:|---:|
| 15 | 23 | 1,248 | 1.000 | 11 | 46.1 k | 107 ms |
| 21 | 27 | 2,390 | 0.793 | 13 | 47 k | 298 ms |
| 35 | 31 | 3,528 | 0.784 | 13 | 50.2 k | 447 ms |

**Phase estimation** (p = 0.1%, cultivated magic states; success = energy within chemical accuracy, 1.6 mHa, of the exact ground state of the Hamiltonian simulated):

| molecule | system qubits | T gates per run | noiseless | d for 90% of it | run time |
|---|---:|---:|---:|---:|---:|
| H2 | 2 | 1.79e7 | 0.92 | 19 | 7.2 min |
| HeH+ | 2 | 4.37e7 | 0.97 | 21 | 20.3 min |
| LiH, CAS(2, 2) | 2 | 4.32e7 | 0.88 | 23 | 21.9 min |

H₂ and HeH⁺ in STO-3G are simulated whole (tapered to 2 qubits). LiH is a deliberately small active space: Li 1s² frozen, two electrons in two CASSCF orbitals, 9 Pauli terms on 2 qubits; its exact energy is 1.3 mHa above full STO-3G FCI at R = 1.595 Å.

**RSA-2048, this work's circuits** (p = 0.1%, cultivated magic states, E[faults] ≤ 0.1; operation counts exact, schedule fitted on compiled instances, error rates extrapolated):

| arithmetic | calibration | logical qubits | Toffolis | d | physical qubits | run time |
|---|---|---:|---:|---:|---:|---:|
| Textbook Cuccaro | separate measure and reset layers | 6,150 | 3.44e11 | 35 | 32.5 M | 1.96 years |
| Textbook Cuccaro | one measure-and-reset layer | 6,150 | 3.44e11 | 33 | 29 M | 1.85 years |
| Modern windowed + MBU | separate measure and reset layers | 8,199 | 8.61e10 | 35 | 43.2 M | 298 days |
| Modern windowed + MBU | one measure-and-reset layer | 8,199 | 8.61e10 | 33 | 38.5 M | 281 days |

**Published estimates under the measured model** (p = 0.1%; each at its authors' failure budget; distances only, since their qubit totals use layouts this model does not include):

| workload | failure budget | published d | d, measured model | d, recalibrated with one measure-and-reset layer |
|---|---:|---:|---:|---:|
| RSA-2048, Gidney 2025 (one shot) | 6.7% | 25 | 29 | 27 |
| FeMoco (THC), Lee et al. 2021 | 1.0% | 31 | 31 | 29 |

Gidney 2025 picks d = 25 because his projection of surface-code memory reaches 1e-15 errors per patch per round there. The measured model (SD6 noise, correlated matching) reaches it at d = 29 (d = 27–29 within the fit's 1σ band). At the measured distances its rates are already 3.4–6.5× his projection, and the gap widens with d because the fitted suppression per two units of distance (Λ = 10.1 and 10.9 for the two bases) is below the projection's 12.25; at d = 25 the fit extrapolates to 1.30e-14 per patch-round (1σ: 9.95e-15–1.70e-14). One shot of his algorithm through the estimator needs d = 29 at his budget (d = 27 at 90%); the same estimator with his projection in place of the measured model returns d = 25. The table below shows how much of that gap is the circuit and noise model rather than his design.

**Which memory model?** Memories at p = 0.1% under other circuits and noise models, each decoded by correlated matching and fitted the same way (per-round rates from d and 3d rounds), against the projection:

| memory model | Λ (Z, X) | measured to d | rate ÷ projection there | d, worse basis ≤ 1e-15 | d, either basis ≤ 1e-15 | d, one shot at his budget |
|---|---:|---:|---:|---:|---:|---:|
| SD6, as calibrated: separate measure and reset layers | 10.8, 10.5 | 9 | 5.9× | 27 | 29 | 27 |
| SD6, one measure-and-reset layer | 12.0, 11.8 | 9 | 3.2× | 27 | 27 | 27 |
| SI1000 on CZ, one measure-and-reset layer | 11.6, 12.3 | 9 | 2.1× | 27 | 27 | 27 |
| uniform depolarizing, Stim's generated circuit | 20.4, 15.9 | 9 | 0.3× | 23 | 23 | 23 |
| Gidney 2025 projection, 3.5^−d/40 | 12.25 | — | 1× | 25 | 25 | 25 |

The distance that reaches 1e-15 per patch-round at p = 0.1% runs from d = 23 to d = 29 across these circuit-level models, and Gidney's d = 25 lies inside that range. SI1000, the noise model of the yoked-surface-code paper's simulations (arXiv:2312.04522), has the projection's slope here (Λ = 11.6–12.3) but a 2.1× higher rate at d = 9, and needs d = 27. Gidney 2025 labels the projection uniform depolarizing; Stim's generated circuit with every parameter at p, one reading of that which leaves qubits idle during gate layers noiseless, is below it. The calibration's circuit is the most pessimistic of the four: it depolarizes idle qubits in every layer and measures and resets the ancillas in separate layers, so data qubits idle through two of them each round. Merging them alone lowers the d = 9 rates 1.7–2.0× and raises Λ to ≈ 12, so every result here leans pessimistic by about that. (The SD6 row is the calibration's own memories fitted like the others, over d ≤ 9; with all its memory lengths and d = 11 the calibrated model needs d = 29.)

**Recalibrated with one measure-and-reset layer.** Every calibration experiment was rerun on the same circuits with each reset layer folded into the measurement layer before it (`tools/calibrate.py --mr`; data/calibration/raw_mr/, model_mr.json). Idling then has Λ = 11.6–12.5, reaches 1e-15 per patch-round at d = 27 (d = 27 at either 1σ edge), and one shot of Gidney's algorithm needs d = 27 at his budget. The tables above show both calibrations.

The simulated algorithm runs were repeated on it too (16 workloads, data/results_mr/, run in parallel on GitHub Actions by `tools/grid.py`). The distance at which a run first reaches 90% of its noiseless success drops by one step (d − 2) for 9 of them and stays the same for 7; none needs more. Where the magic states set the ceiling (injected states, H₂ with 15-to-1), the ceiling is unchanged.

| workload | p | noiseless | d for 90% of it, as calibrated | one measure-and-reset layer |
|---|---:|---:|---:|---:|
| Shor, textbook, N = 15 | 0.1% | 1.00 | 13 | 11 |
| Shor, textbook, N = 15 | 0.2% | 1.00 | 19 | 17 |
| Shor, textbook, N = 21 | 0.1% | 0.79 | 13 | 13 |
| Shor, textbook, N = 21 | 0.2% | 0.79 | 21 | 19 |
| Shor, textbook, N = 35 | 0.1% | 0.79 | 15 | 13 |
| Shor, textbook, N = 35 | 0.2% | 0.79 | 21 | 21 |
| Shor, textbook, N = 77 | 0.1% | 0.77 | 15 | 15 |
| Shor, textbook, N = 77 | 0.2% | 0.77 | 23 | 21 |
| Shor, textbook, N = 143 | 0.1% | 0.77 | 15 | 15 |
| Shor, textbook, N = 143 | 0.2% | 0.77 | 23 | 21 |
| Shor, windowed + MBU, N = 15 | 0.1% | 1.00 | 11 | 11 |
| Shor, windowed + MBU, N = 21 | 0.1% | 0.79 | 13 | 13 |
| Shor, windowed + MBU, N = 35 | 0.1% | 0.78 | 13 | 13 |
| H2 | 0.1% | 0.92 | 19 | 17 |
| HeH+ | 0.1% | 0.97 | 21 | 19 |
| LiH, CAS(2, 2) | 0.1% | 0.88 | 23 | 21 |
<!-- results:end -->

## How it is checked

- **ftsim** against Qiskit Aer (100 random circuits' amplitudes; exact noisy distributions, plain and
  stratified), dense against sparse, WebAssembly against native.
- **Arithmetic** verified on every input by an independent bit-level simulator; Shor's noiseless output
  against the exact phase-estimation distribution; every synthesized rotation as a matrix.
- **Windowed arithmetic** verified on every input the same way, including the measurement-based
  carry uncomputation and its feed-forward CZ fixups.
- **Scaling**: the Toffoli and qubit counts used at RSA-2048 equal the compiled circuits' for every
  simulated N; the estimator, given Gidney 2025's own error projection and budget, returns his d = 25.
- **Chemistry** against PySCF (integrals, HF and FCI to 1e-8 Ha) and Szabo & Ostlund; the LiH
  active-space Hamiltonian's ground energy equals PySCF's CASSCF(2, 2) to 1e-9 Ha.
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
.venv/bin/python tools/noise_models.py && .venv/bin/python tools/noise_models.py --summary   # SI1000, uniform (~1 h)
.venv/bin/python tools/run_shor.py && .venv/bin/python tools/run_shor_modern.py
.venv/bin/python tools/run_qpe.py && .venv/bin/python tools/run_lih.py
.venv/bin/python tools/scale.py && .venv/bin/python tools/sensitivity.py
.venv/bin/python tools/site_data.py && .venv/bin/python tools/report.py && .venv/bin/python tools/readme_tables.py
cargo build --release --target wasm32-unknown-unknown --lib && cp target/wasm32-unknown-unknown/release/ftsim.wasm site/wasm/
```

Every dataset records its seed and the versions that produced it.

## Assumptions

Magic-state error rates and factory footprints (Litinski 2019; Gidney, Shutty & Jones 2024; Gidney
2025), the H and S operation model of the fast-block layout, 1 µs rounds and a 10 µs decoder reaction
time are cited inputs, not measured here. Logical channels are measured up to d = 11 (idling) and
d = 9 (surgery) at p = 0.1%; runs at larger distances use extrapolated channels and are marked as
such, and every RSA-2048 and FeMoco number is an extrapolation. The calibrated circuit is
stabilizer-qec's SD6 memory, which measures and resets ancillas in separate layers; with one
measure-and-reset layer the same noise errs about half as often at d = 9 (see "Which memory model?"),
so results lean pessimistic. See the report's Limitations.

MIT licensed. To cite it, see [CITATION.cff](CITATION.cff).
