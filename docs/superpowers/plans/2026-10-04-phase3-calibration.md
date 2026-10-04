# Phase 3: calibration on `stabilizer-qec` — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Measured logical Pauli channels of a rotated surface-code patch (idle per round, full
X/Y/Z), of lattice-surgery Z⊗Z / X⊗X measurements and of the surgery CNOT, at SD6 noise, with
fits in d and p, a Stim + PyMatching cross-check, and composition checks.

**Architecture:** `ftalgo.calib.qec` samples and decodes circuits with stabilizer-qec (correlated
matching) and returns joint observable-flip counts. Channels are handled through their Pauli
fidelities f(s) = E[(−1)^{s·flips}], which multiply under composition: per-round channels come
from the ratio of two lengths, and an operation's own channel is the measured one divided by
its baseline (preparation, readout and idle rounds, measured by memory experiments of the
same length). `tools/calibrate.py` runs the grid and records raw counts in
`data/calibration/*.json`; `ftalgo.calib.model` fits and serves channels at any (d, p).

**Tech Stack:** stabilizer-qec 1.2 (engine under test), stim 1.16 + PyMatching 2.4 (cross-check).

## Global Constraints

- SD6 noise (stabilizer-qec `noise="sd6"`), correlated matching, rotated surface code.
- Every recorded dataset stores shots, seed, counts per flip pattern, versions and the command.
- Values outside the directly measured range are marked extrapolated wherever shown.

## Tasks

### Task 1: sampling/decoding helper + flip-pattern counts
`ftalgo/calib/qec.py`: `run(circuit, *, shots, seed, target=None, batch=65536, correlated=True)
-> Counts(shots, patterns: dict[int,int], num_obs)`; stops early once `target` shots with any flip
are seen. Test: a d=3 memory at p=0.5% gives a failure rate equal (within 5σ) to a direct
sample/decode of the same shots.

### Task 2: Pauli-fidelity algebra
`ftalgo/calib/channels.py`: `fidelities(counts) -> dict[s, f]`, `probabilities(fid) -> pattern
probs` (inverse Walsh–Hadamard), `per_round(fid_T1, fid_T2, T1, T2)`, `deconvolve(fid_M, fid_B)`,
`pauli_channel_2(xpart, zpart) -> 15-vector` (Stim order). Tests: round trips on random
channels; composing two channels by convolution equals multiplying fidelities.

### Task 3: reference-qubit (Choi) memory circuit
`ftalgo/calib/choi.py`: `choi_memory(d, rounds, p) -> text`. Tests: stabilizer-qec builds its DEM
(all detectors and both observables deterministic); noiseless circuit gives no flips; with
noise, P(obs0) and P(obs1) agree with X- and Z-basis memory experiments of the same length
within 5σ.

### Task 4: experiment grid + Stim/PyMatching cross-check
`tools/calibrate.py` (resumable, one JSON per experiment): memory X/Z at T = d and 3d, Choi
memory at T = d and 3d, CNOT (z and x inputs, merged = d) with 2d-round memory baselines,
Z⊗Z (basis z: outcome, Z1, Z2) and X⊗X (basis x) with pre = post = d and 2d-round baselines,
repeated_zz k = 2, 4 and line n = 3 for composition; d ∈ {3, 5, 7, 9, 11} (as memory allows),
p ∈ {0.1, 0.2, 0.3, 0.5}%; target 2000 failures or a shot cap. `tools/xcheck.py` reruns a subset
with Stim + PyMatching.

### Task 5: model and fits
`ftalgo/calib/model.py`: `LogicalModel.load(path)`, per component `rate(d, p)` with
ε = A·(p/p*)^((d+1)/2) fitted (weighted least squares on log ε over measured points with ≥ 20
failures), `idle(d,p) -> (px,py,pz)` per round, `cnot(d,p) -> 15-vector`, `zz(d,p)`, flags for
extrapolation. Tests on recorded data: fit residuals within uncertainty; extrapolation monotone.

### Task 6: composition checks
`ftalgo/calib/compose.py`: predict repeated_zz and line outcomes, and the CNOT experiment from
Z⊗Z/X⊗X components + idle; compare with measurement (ratio and σ) → `data/calibration/composition.json`.
