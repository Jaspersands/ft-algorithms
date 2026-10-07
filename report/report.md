---
title: Factoring and chemistry on a simulated fault-tolerant quantum computer
subtitle: Algorithms run end to end on a surface-code machine whose logical error rates were measured by circuit-level simulation
author: Jasper Sands
date: 2026-10-08
version: 0.1
commit: 68dbe11
description: Technical report of the ft-algorithms project.
toc: true
abstract: |
  Shor's algorithm (full modular exponentiation in Toffoli arithmetic, semiclassical
  approximate QFT) and iterative phase estimation of the H₂ and HeH⁺ ground-state energies are
  run end to end on a simulated rotated-surface-code machine. Every logical error channel, idling,
  lattice-surgery merges and the surgery CNOT, was measured in 316 circuit-level experiments
  (0.82×10⁹ shots) with stabilizer-qec under SD6 noise and extracted with baselines
  propagated through each operation; magic states come from published factory models. Whole
  surgery experiments predicted from their parts agree with circuit-level measurements to
  0.93–1.26× (median 1.04×). At p = 0.1%, factoring 15
  reaches 90% of its noiseless peak probability at d = 13; phase estimation
  of H₂ needs 1.8×10^7^ T gates per run. Extrapolated with the same error model,
  textbook-arithmetic RSA-2048 needs d = 35 and 32.5 M
  physical qubits for 1.96 years; Gidney's 2025 algorithm needs
  d = 29 under our error model where he assumed 25.
---

# Introduction

Resource estimates for fault-tolerant algorithms usually multiply a gate count by an assumed
logical error rate. This work runs the algorithms instead. A logical-level simulator (ftsim,
written for this project) executes the actual circuits, with every operation followed by its
measured logical channel and every waiting patch accumulating measured idle error. The results
are success probabilities as functions of code distance, physical error rate and magic-state
source, with error budgets that say which operations cause the failures.

The approach rests on one fact and one approximation. The fact: under Pauli circuit noise, for a
Clifford circuit decoded into a Pauli-frame correction, the logical action of one surface-code
operation is exactly a Pauli channel, since every syndrome history ends in some logical Pauli. The
approximation: consecutive operations' channels are composed as independent. Section 4.3
measures how good that is.

Contributions:

- A calibrated logical model of a rotated surface code under SD6 noise: idle (including logical Y),
  Z⊗Z and X⊗X merges and the surgery CNOT, from 316 experiments, with fits in d and p.
- A composition test: whole lattice-surgery experiments predicted from parts at the logical level
  and compared with circuit-level measurement.
- Shor's algorithm with honest arithmetic (no use of the period anywhere in the circuit) run on
  that model for N up to 143, scored by peak probability because, for small N, uniformly random
  outcomes already "factor" with probability 1.00 (N = 15).
- Molecular phase estimation from first-principles Hamiltonians, to chemical accuracy.
- Extrapolation to RSA-2048 and FeMoco under the same error model, compared with published estimates.

# Methods

## Circuit-level calibration

All physical-level simulation uses stabilizer-qec 1.2 (the author's engine, which matches Stim and
PyMatching) with SD6 noise at p ∈ {0.1, 0.2, 0.3, 0.5}% and correlated matching. Experiments: memories
of d, 2d, 3d and 4d rounds in both bases (d = 3 … 11), reference-qubit memories (d ≤ 7), Z⊗Z and X⊗X
merges of d rounds with d rounds apart before and after, the surgery CNOT (two merges and an ancilla
patch, 4d rounds in all), and, for the composition test, three Z⊗Z measurements in a row and a
three-patch merge. Each ran to 500 failures or a 300-second cap; every record keeps its seed.

**Channels through fidelities.** A distribution over flip patterns has Pauli fidelities
f(s) = Σₓ P(x)(−1)^(s·x), which multiply under composition. The per-round idle channel is the slope
of log f over the memory lengths (preparation and readout cancel). An operation's own channel is the
measured one divided by its baseline, where the baseline (preparation plus d rounds before, d rounds
plus readout after, each the square root of a 2d-round memory's fidelity) is first pushed through the
ideal operation: a CNOT copies an X error on the control onto the target and a Z on the target onto the
control, and an error before a merge changes its outcome. Dividing naively gives negative
probabilities for "X on the control only". Uncertainties are bootstrap (300 resamples).

**Logical Y.** X- and Z-basis memories see X∪Y and Z∪Y. A reference-qubit experiment gives the full
channel: the patch is projected noiselessly (every Z check by a noiseless Pauli-product measurement,
which round 1 is compared with), entangled with a perfect qubit by a noiseless Z_R Z_L measurement,
run, and ended by a perfect round and noiseless X_R X_L and Z_R Z_L measurements. Without the initial
projection an X error before the first Z-check readout flips Z_R Z_L unseen.

**Fits.** Each component ε(d, p) is fitted per p as log ε = a + b(d + 1)/2 by weighted least squares
over points with at least ten events, and globally as A(p/p*)^((d+1)/2). Values beyond the largest
measured distance are extrapolations and are marked as such in every figure.

## The machine

| operation | rounds | logical channel |
| :--- | ---: | :--- |
| X, Y, Z | 0 | none (Pauli frame) |
| H | 0 | none: transversal, orientation tracked in the fast block † |
| S, S† | d | Z⊗Z merge channel (|Y⟩ resource); wrong outcome → Z |
| T, T† | d | merge channel + Z w.p. ε_T † + wrong outcome → Z w.p. ½ (twirled S) |
| CX, CZ | 2d | measured surgery-CNOT channel, idling included |
| CCX, CCZ | 2d | 3 merges + Z-type w.p. ε_CCZ † + wrong outcomes → twirled CZ; d rounds idle |
| M, R | 1 | one round of idling |
| idle | per round | measured (p_X, p_Y, p_Z) |
| feed-forward | 10 | decoder reaction time before a TABLE † |

† marks cited inputs. Magic states: injection (ε_T = p, conservative against Li 2015's ≈ 0.4p), 15-to-1
distillation (Litinski 2019, Table 1: 4.5×10⁻⁸ per T and 5.2×10⁻¹¹ per CCZ at p = 0.1%) and cultivation
(Gidney, Shutty and Jones 2024: 2×10⁻⁹ at p = 0.1%, 4×10⁻¹¹ at 0.05%; |CCZ⟩ by 8T-to-CCZ, 28 ε_T²). Physical
qubits: Litinski's fast block, 2n + ⌈√(8n)⌉ + 1 tiles of 2(d + 1)² qubits, plus factories sized by their
published qubit·rounds per state. One round is 1 µs.

The noise compiler gives every qubit a clock. An operation starts when its qubits are free (a TABLE
also waits for its records plus the reaction time); each waiting qubit receives the composed idle
channel for its gap; the operation's channel follows it with a class label for attribution. REPEAT
blocks are barriers, scheduled once; TABLE cases are padded to the longest, so the schedule does not
depend on outcomes. Channel composition is written with log1p/expm1: per-round rates reach 10⁻²⁰ at
large d and a naive 1 − (1 − 2q)ⁿ cancels to zero.

## The simulator

ftsim (Rust) runs programs in a Stim-like format: Clifford + T + Toffoli gates, Pauli channels on
one to three qubits, measurements, REPEAT, and TABLE (classical feed-forward choosing a block by
earlier records). Backends: a sparse state (hash map of basis states; runs of permutation and phase
gates fused per basis state), which handles 36-qubit Shor circuits in milliseconds per shot, and a
dense vector. Each shot's RNG derives from (seed, shot), so results do not depend on threads.

**Stratified sampling.** With K the number of static faults, P(success) = Σₖ P(K = k)·Sₖ. P(K = k) is
computed exactly (a truncated generating polynomial over the site tree, REPEAT by exponentiation), and Sₖ
is estimated from shots with exactly k faults, drawn by sampling k sites with weights q/(1 − q) and
rejecting duplicates, which is exactly the conditional law. Strata beyond the tail cutoff carry their
mass as an interval. The one-fault stratum, with faults logged, gives each operation class's harm.

## Algorithms

**Shor.** Cuccaro ripple-carry adders; modular addition of a classical constant with a sign flag
(Vedral–Barenco–Ekert, Beauregard); controlled modular multiplication with uncomputation by the inverse
and a controlled swap; 3n + 6 qubits and 2n(2n(10n + 12) + n) Toffolis. One control qubit is recycled
through 2n rounds of a semiclassical approximate QFT whose correction rotations, chosen by TABLE from
the previous ⌈log₂ 2n⌉ + 2 records, are synthesized by gridsynth to ε = 10⁻³/2n and verified as matrices.
The base is a = 2 for every N. Every multiplier is the generic circuit, even when its constant is 1.

**Chemistry.** STO-3G integrals of s Gaussians in closed form (Boys F₀), restricted Hartree–Fock with
DIIS, Jordan–Wigner with interleaved spins, and Z₂ tapering of both spin parities (Bravyi et al.).
Iterative phase estimation reads 10 bits of U = exp(−i(H − E_HF)τ) with τ = 2π, approximated by four
fourth-order Suzuki steps; controlled Pauli rotations use only uncontrolled synthesized Rz's, so their
global phases stay global.

# Results

## The calibrated machine

316 experiments, 0.82×10⁹ shots and 459,686 logical failures. Suppression
factors Λ (per two units of distance):

| p | idle (X∪Y) | idle (Z∪Y) | CNOT X on target | CNOT Z on control | Z⊗Z patch Z | Z⊗Z wrong outcome |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.1% | 10.93 | 10.06 | 7.91 | 7.10 | 7.12 | 12.30 |
| 0.2% | 4.86 | 4.72 | 3.42 | 3.43 | 3.17 | 6.33 |
| 0.3% | 2.80 | 2.85 | 2.01 | 2.01 | 1.98 | 2.97 |
| 0.5% | 1.55 | 1.62 | 1.13 | 1.11 | 1.11 | 1.58 |

![Per-round logical error of an idle patch (X∪Y). Points measured (1σ), lines per-p fits, dashed beyond the largest measured distance.](figures/idle.svg)

At p = 0.1% an idle patch errs 8.8×10^-7^ per round at d = 9, and the fit gives 4.3×10^-15^
at d = 25. Logical Y is 4% of the X∪Y rate at d = 3, 43× what
independent X and Z would give: Y errors are included as measured.

**Cross-check.** 36 experiments rerun with Stim and PyMatching (correlated) agree with
stabilizer-qec within statistics (largest |z| = 2.2).

## Composition

![Observable failure rates of whole surgery experiments: predicted from calibrated parts by ftsim vs measured at circuit level; dashed ±25%.](figures/composition.svg)

Across 48 experiment–noise–distance combinations the predicted failure rate is
0.93–1.26× the measured one, median 1.04×. The independent-composition
model is accurate to about 10% and leans pessimistic.

## Shor's algorithm

At p = 0.1% with cultivated magic states:

| N | n | logical qubits | Toffolis | T | noiseless peak | random peak | random factors | d for 90% | physical qubits | run time |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 15 | 4 | 18 | 3,360 | 223 | 1.000 | 0.0156 | 0.996 | 13 | 55.4 k | 225 ms |
| 21 | 5 | 21 | 6,250 | 317 | 0.790 | 0.0059 | 0.999 | 13 | 58.1 k | 418 ms |
| 35 | 6 | 24 | 10,440 | 425 | 0.790 | 0.0029 | 1.000 | 15 | 63.4 k | 808 ms |
| 77 | 7 | 27 | 16,170 | 519 | 0.774 | 0.0018 | 1.000 | 15 | 67 k | 1.25 s |
| 143 | 8 | 30 | 23,680 | 635 | 0.774 | 0.0009 | 1.000 | 15 | 70.8 k | 1.81 s |

![Left: peak probability vs d (cultivated states, p = 0.1%); dotted: noiseless. Right: the same runs scored by factors found; dotted: random-outcome baseline.](figures/shor.svg)

Factoring success is a poor score for small N. Classical post-processing that tries the
convergents' small multiples finds the factors from most outcomes, and the noisiest runs look best.
The peak probability falls to its random baseline (r/2^m) instead.

![N = 35: peak probability vs d for each magic-state source; solid p = 0.1%, dashed p = 0.2%.](figures/factories.svg)

| magic states | p | ε_T † | ε_CCZ † | best peak | d for 90% |
| :--- | ---: | ---: | ---: | ---: | ---: |
| cultivated | 0.1% | 2.0×10^-9^ | 1.1×10^-16^ | 0.794 (d = 17) | 15 |
| cultivated | 0.2% | 10.0×10^-8^ | 2.8×10^-13^ | 0.798 (d = 29) | 21 |
| 15-to-1 distilled | 0.1% | 4.5×10^-8^ | 5.2×10^-11^ | 0.799 (d = 19) | 15 |
| 15-to-1 distilled | 0.2% | 3.6×10^-7^ | 3.3×10^-9^ | 0.793 (d = 27) | 21 |
| injected | 0.1% | 1.0×10^-3^ | 2.8×10^-5^ | 0.612 (d = 17) | never |
| injected | 0.2% | 2.0×10^-3^ | 1.1×10^-4^ | 0.366 (d = 23) | never |

**Error budget** (N = 35, d = 15, cultivated, p = 0.1%):

| class | expected faults | share | score of one-fault runs |
| :--- | ---: | ---: | ---: |
| idle | 0.027 | 80.3% | 0.394 |
| CNOT | 5.3×10^-3^ | 15.9% | 0.564 |
| CCZ | 1.2×10^-3^ | 3.7% | 0.667 |
| T | 1.7×10^-5^ | 0.0% | — |
| S | 1.0×10^-5^ | 0.0% | — |
| measure | 3.7×10^-8^ | 0.0% | — |

Idling accounts for 80% of the expected faults. The textbook circuit runs one Toffoli at a
time (the schedule averages 5.13·d rounds per Toffoli with 3.2 patches busy)
while every other patch waits.

## Phase estimation

| molecule | logical qubits | terms | T gates | E_HF (Ha) | E_FCI (Ha) | noiseless | best | d for 90% | run time | physical qubits |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| H2 | 3 | 5 | 1.8×10^7^ | -1.116714 | -1.137276 | 0.92 | 0.91 (d = 23) | 19 | 7.2 min | 22.1 k |
| HeH+ | 3 | 9 | 4.4×10^7^ | -2.841836 | -2.851466 | 0.97 | 0.91 (d = 23) | 21 | 20.3 min | 22.4 k |

![Probability of chemical accuracy (1.6 mHa from FCI) vs d at p = 0.1%; dotted: noiseless.](figures/qpe.svg)

Each run uses ≈ 1.8×10^7^ T gates for H₂, so the magic-state error sets a floor
that no distance removes. With cultivated states (2×10⁻⁹) that floor is a few percent; with 15-to-1
distilled states (4.5×10⁻⁸) most runs fail.

![The answer itself: energies returned by 400 individual H₂ runs at several distances (cultivated states, p = 0.1%), as the error from the exact energy; shaded: chemical accuracy; grey: the same circuit without noise.](figures/qpe_hist.svg)

## To scale

Operation counts of the textbook circuit are exact; the schedule's shape (rounds per Toffoli in units
of d, 5.13; CNOTs per Toffoli, 2.28; patches busy, 3.18) is fitted on compiled
instances and is stable across them. The smallest odd d with at most 0.1 expected faults:

| workload | logical qubits | Toffolis | published d | published qubits | our d | our qubits | run time |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| RSA-2048, textbook arithmetic | 6,150 | 3.4×10^11^ | — | — | 35 | 32.5 M | 1.96 years |
| RSA-2048, Gidney 2025 | 1,409 | 6.5×10^9^ | 25 | 898 k | 29 | 5.31 M | 4.63 days |
| FeMoco (THC), Lee et al. 2021 | 2,142 | 5.3×10^9^ | 31 | ≈ 4 M | 29 | 8.01 M | 2.45 days |

![Textbook Shor from the simulated sizes to RSA-2048 (p = 0.1%).](figures/scaling.svg)

Under our measured error model, Gidney's 2025 RSA-2048 algorithm needs d = 29 where he assumed
25: our idle rate at d = 25 is 4.3×10^-15^ per round against his 10⁻¹⁵. Our physical-qubit totals use
plain patches throughout and so bound his yoked-storage design from above. The textbook arithmetic's
3.4×10¹¹ Toffolis exceed what 15-to-1 distillation's 5.2×10⁻¹¹ can support beyond n = 256.

# Validation

- ftsim: dense = sparse on random circuits; 100 random circuits equal Qiskit's statevector to 10⁻¹⁰;
  noisy distributions equal Aer's exact density matrices (χ², plain and stratified); WebAssembly = native.
- Arithmetic verified on every input by an independent bit-level simulator; Shor's noiseless outcomes
  match the exact distribution; synthesized rotations verified as matrices; QPE matches its exact
  Trotterized distribution.
- Chemistry: integrals, HF and FCI equal PySCF to 10⁻⁸ Ha; Szabo and Ostlund's H₂ reproduced; tapering
  preserves the spectrum.
- Calibration: Stim + PyMatching cross-check; reference-qubit experiment against plain memories;
  composition test.

# Limitations

- Composition is accurate to about 10% (Section 4.2), measured only for d ≤ 7.
- Most algorithm results lie beyond the directly measured distances (d ≤ 9 for the CNOT at p = 0.1%):
  extrapolated along measured exponential fits, marked in the figures.
- Magic-state error rates and footprints, the H and S operation model and the timing constants are
  cited inputs, not measured here.
- The arithmetic is textbook and sequential. Modern constructions (windowed arithmetic, measurement-based
  uncomputation, parallel Toffolis) would cut the cost by orders of magnitude; the scaling section
  separates the error model from the algorithm by also evaluating Gidney's 2025 counts.

# Reproduction

```
tools/build.sh                     # the engine (Rust → Python module)
python tools/calibrate.py          # circuit-level experiments (resumable)
python tools/fit_model.py          # data/calibration/model.json
python tools/compose_check.py; python tools/xcheck.py
python tools/run_shor.py; python tools/run_qpe.py; python tools/scale.py
python tools/site_data.py; python tools/report.py
```
