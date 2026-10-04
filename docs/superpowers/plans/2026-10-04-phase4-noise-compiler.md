# Phase 4: noise model, scheduling and experiments — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Turn an ideal logical circuit into a scheduled noisy ftsim program for a surface-code
machine at (d, p, factory), estimate algorithm success with fault-count stratification, and
attribute failures to operation classes.

**Architecture:** `ftalgo.arch` holds the operation model (durations, channels from the
calibrated `LogicalModel`, magic-state factories as cited inputs) and the physical footprint.
`ftalgo.schedule` walks the IR (gates, measurements, TABLE, REPEAT) with per-qubit clocks
(ASAP, reaction time before feed-forward, barriers at REPEAT), emitting each op, its channel and
idle channels for every gap, and accumulates expected faults per class. `ftalgo.estimate` runs
stratified sampling and returns success with error bars and the error budget.

## Operation model (documented in `arch.py`; † = cited input, not measured here)

| op | rounds | channel |
|---|---|---|
| X, Y, Z | 0 | none (Pauli frame) |
| H | 0 | none: transversal, boundary orientation tracked (fast block) |
| S, S_DAG | d | Z⊗Z-type merge channel on the patch (Y-state consumption) |
| T, T_DAG | d | merge channel on the patch + Z with ε_T † + wrong outcome → Z w.p. ½ |
| CX, CZ, CY | 2d | measured surgery CNOT channel (includes idling) |
| SWAP | 6d | three CNOT channels |
| CCX, CCZ | 2d | three merge channels + Z-type with ε_CCZ † + wrong outcome → twirled CZ on the other two; then d rounds idle |
| M, MX | 1 | one round of idle |
| R, RX | 1 | one round of idle |
| idle | per round | measured (p_X, p_Y, p_Z) |

Reaction time: 10 rounds between a measurement and a TABLE that reads it.

## Tasks

1. `arch.py`: `Factory` catalog (injected, 15-to-1, cultivation; ε_T(p), ε_CCZ(p), qubit·rounds per
   state, sources), `Architecture(model, d, p, factory)` with op specs; tests: channels are valid
   probability vectors, monotone in d.
2. `schedule.py`: `compile_noisy(circuit, arch) -> NoisyProgram(text, rounds, budget, counts,
   qubits)`; tests: idle gap noise equals composed idle; TABLE padding; REPEAT barrier and
   duration; expected-fault budget equals ftsim's `expected_faults`; a Clifford-only program
   under the noise matches a direct Pauli-frame calculation.
3. `estimate.py`: `stratified(program, success_fn, shots_per_k, kmax|tail)`; tests: equals plain
   sampling within error on a small noisy program; P(k) tail accounting.
4. `tools/run_shor.py`, `tools/run_qpe.py`: grids → `data/results/*.json` (seeded, resumable).
