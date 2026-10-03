# Algorithms on a simulated fault-tolerant quantum computer: design

Date: 2026-10-04. Status: approved for autonomous build ("Go ahead and do it. Plan/execute. Do big
things."), with open questions self-answered by the standing rule: take the fullest option unless
it is extremely slow or needs very large storage.

## Goal

Run real quantum algorithms, end to end, on a simulated surface-code quantum computer whose
logical error rates are **measured** by circuit-level simulation, and report what it takes for
them to give the right answer. The headline results are curves: probability that Shor's
algorithm factors N, and probability that phase estimation finds a molecule's ground-state energy
to chemical accuracy, as functions of code distance d and physical error rate p, with an error
budget saying which operations cause the failures. The same measured numbers then extrapolate to
cryptographically and chemically relevant sizes.

The project sits beside the QEC simulator (`stabilizer-qec`, PyPI 1.2.0) and keeps its identity:
every number is measured or computed, every component is cross-checked against an independent
tool, every assumption is stated where the number is shown, and results reproduce from a seed.

Non-goals: new decoders, new codes, hardware-specific compilation, claiming any speedup.

## The stack

```
physical  SD6 circuit noise at p           stabilizer-qec: sample + decode (Stim-format circuits)
   │      memory, lattice-surgery ZZ/XX, CNOT, reference-qubit (Choi) experiments
   ▼
logical   calibrated Pauli channels per operation and per round, fitted ε(d, p)
   │      magic states from published factory models (marked †)
   ▼
program   Clifford + T + Toffoli(CCZ) logical circuits: Shor (full modular exponentiation),
   │      phase estimation for H2 / HeH+ / H3+ / H4 Hamiltonians we compute ourselves
   ▼
run       ftsim (Rust): sparse and dense state vectors, Pauli-channel Monte Carlo, classical
          feed-forward, seeded and thread-count independent; scheduled with op durations so
          every live qubit idles honestly
```

### Why per-operation Pauli channels are the right abstraction

With Pauli circuit noise, a Clifford circuit, and a decoder that applies a Pauli correction,
the logical action of one surface-code operation is exactly a Pauli channel on its logical
qubits (each syndrome history yields a logical Pauli; averaging gives a stochastic mixture). The
one approximation is **composition**: treating the channels of consecutive operations as
independent. Correlations across operation boundaries come from the decoder's view of the
boundary rounds. This is checked directly: sequences simulated physically (repeated Z⊗Z,
multi-patch merges, the full lattice-surgery CNOT against its components) are compared with the
composed prediction.

Non-Clifford gates are implemented by consuming magic states (|T⟩ for rotations, |CCZ⟩ for
Toffolis). After twirling by each state's stabilizer group, an imperfect magic state is the ideal
one followed by a Z-type Pauli error with the state's error probability ε, so the gate is the
ideal gate followed by a Z-type channel; consumption itself is a lattice-surgery measurement
whose channel is calibrated. ε comes from published factory models (injection, 15-to-1,
cultivation, 8T-to-CCZ), cited and marked †; these are the only non-measured error rates.

## Components

### 1. `ftsim` — the logical-level simulator (Rust crate, PyO3 module, WASM build)

- **Program format**: a Stim-like text format for logical circuits. Gates: `I X Y Z H S S_DAG T
  T_DAG SQRT_X SQRT_X_DAG CX CY CZ SWAP CCX CCZ`; `M MX MY R RX`; noise `PAULI_CHANNEL_1(px,py,pz)`,
  `PAULI_CHANNEL_2(15 probabilities, Stim order)`, `Z_CHANNEL_3(7 probabilities)` (Z-type errors on
  three qubits, for CCZ); `REPEAT n { }`; feed-forward `TABLE(rec[-k] ... ) { case0 } { case1 } ...`
  selecting one block by the bits of earlier measurements (semiclassical QFT corrections and
  phase-estimation feedback); `TICK` as a scheduling marker; `MARK(label)` to attribute faults to
  operation classes.
- **Backends**: sparse (basis states as `u128` keys in a hash map, for permutation-heavy
  arithmetic where few amplitudes are non-zero) and dense (`Vec<Complex64>`, up to ~28 qubits).
  The same program runs on both and they agree to machine precision.
- **Monte Carlo**: per-shot RNG derived from (seed, shot index), so results are identical for any
  number of threads; rayon across shots. Two sampling modes: plain (every noise location fires
  independently) and **stratified by fault count** (the exact probability of exactly k faults
  from a Poisson-binomial recursion, and shots conditioned on k faults), which resolves success
  probabilities to tight error bars from the noiseless regime to the failing one.
- **Fault attribution**: each sampled fault carries its `MARK` class, so failed shots can be
  attributed to memory, CNOT, T, CCZ, measurement, etc.
- **Outputs**: measurement records per shot (bit-packed numpy arrays) and fault logs.
- **Checks**: unit tests per gate against matrices; dense vs sparse; against Qiskit Aer
  (noiseless statevector amplitudes on random circuits; noisy output distributions under the same
  Pauli channels, by exact density-matrix simulation and a chi-square test on our samples).

### 2. `ftalgo` — compiler and algorithms (Python)

- **Circuit builder**: a small IR with registers, emitting `ftsim` text; gate counts and depth.
- **Synthesis**: Rz(θ) into Clifford+T by Ross–Selinger gridsynth (`pygridsynth`), cached per
  (θ, ε); every sequence verified as a matrix to be within ε of the target.
- **Arithmetic** (Toffoli based, no rotations): Cuccaro ripple-carry adder, constant loading,
  comparator, controlled modular addition of a classical constant (VBE/Beckman structure with a
  flag uncomputed), controlled modular multiplication by a constant with uncomputation by its
  inverse and a controlled swap. Verified by exhaustive classical simulation of the reversible
  circuit for every input at small n (a bit-level simulator, independent of `ftsim`).
- **Shor**: order finding with one recycled control qubit (semiclassical QFT, Griffiths–Niu),
  Coppersmith's approximate QFT (corrections below π/2^m dropped, m chosen so the noiseless
  success stays within a stated margin of the exact QFT), classical continued fractions and gcd.
  The circuit uses only classically computed constants a^(2^k) mod N, never the order: no
  "compiled Shor". Instances: N = 15, 21, 35, 51, 77, 143, 221, 247 and as far as runtime allows
  (n up to ~10 bits with the sparse backend).
- **Chemistry**: our own STO-3G integrals for s-orbital molecules (overlap, kinetic, nuclear
  attraction and electron repulsion in closed form via the Boys function), restricted Hartree–Fock,
  second quantization, parity mapping with two-qubit reduction. Molecules: H2, HeH+, H3+, H4.
  Checked against Szabo & Ostlund's published H2 / HeH+ numbers and against PySCF (integrals,
  HF and FCI energies).
- **Phase estimation**: iterative QPE (one ancilla, feedback by `TABLE`), second-order Trotter
  controlled evolutions with synthesized rotations, Hartree–Fock input state. Parameters chosen so
  the noiseless circuit reaches chemical accuracy (1.6 mHa) with high probability; Trotter and
  finite-bit bias reported separately from logical noise.
- **Stretch**: Hamiltonian simulation of a Heisenberg chain (observable error vs d), if time allows.

### 3. Calibration (Python on `stabilizer-qec`)

Physical model: SD6 (standard circuit-level depolarizing) at p ∈ {0.05, 0.1, 0.2, 0.3, 0.5}%.

- **Idle memory**: per-round X and Z logical flip rates of a rotated patch, from memory
  experiments at several lengths (slope fit, removing preparation and readout), d = 3..11.
- **X/Z correlation (logical Y)**: a reference-qubit (Choi) memory experiment: our own SD6 rotated
  memory generator, a noiseless reference qubit entangled by a noiseless Pauli-product
  measurement, and noiseless joint readout of X_R X_L and Z_R Z_L, giving the full single-qubit
  Pauli channel. The generator is checked against `stabilizer-qec`'s memory circuit (same logical
  error rate within statistics).
- **Lattice surgery**: Z⊗Z and X⊗X measurement channels (outcome flip, errors on each patch) vs
  merged rounds; the CNOT channel measured directly (`surgery.cnot`, both input bases) and
  predicted from its components (composition check); `repeated_zz` and `line` for longer
  sequences.
- **Cross-check**: the same circuits sampled by Stim and decoded by PyMatching agree within
  statistics.
- **Fits**: ε(d, p) = A·(p/p*)^((d+1)/2) per channel component, with bootstrap uncertainties;
  values beyond the directly measured range are marked extrapolated wherever shown.

### 4. Noise model, schedule and experiments

- **Operation model** (durations in syndrome rounds, 1 round = 1 µs; each op cites its source):
  Paulis free (tracked in the frame); CNOT by surgery (measured channel and its 2d + O(1) rounds);
  T by |T⟩ consumption (Z⊗Z with a magic patch, ε_T, S fix-up half the time); Toffoli by |CCZ⟩
  consumption (ε_CCZ); H and S by patch-level procedures modeled as idle-equivalent rounds
  (stated, not measured); transversal measurement and preparation.
- **Schedule**: ASAP over qubit availability with ample routing (Litinski fast-block style);
  every live qubit receives the calibrated idle channel for every round it waits. Factory count
  sized so magic-state supply never stalls; its footprint counted in the qubit total.
- **Experiments**: for each algorithm instance, success probability vs d at each p and for each
  factory model; error budget; physical qubits and wall-clock time. Smallest d at which each
  instance reaches its noiseless success.
- **Scaling**: our compiled gate counts as functions of n (exact counts from the generator, not
  simulation) with the fitted ε(d, p) give d, qubits and time for RSA-2048 with our textbook
  arithmetic, and with Gidney's 2025 published counts as an alternative input to validate the
  estimator against his ~1M qubits / <1 week. Chemistry scaled the same way with published
  FeMoco counts (†).

### 5. Presentation

- **Site** (static, paper-figure style like the QEC site): the stack, calibration, composition
  check, Shor, phase estimation, error budgets, scaling, methods and cross-checks. An interactive
  demo runs Shor on N = 15 in the browser (`ftsim` compiled to WebAssembly) at the reader's p, d
  and factory, drawing the outcome histogram live.
- **Technical report** built from the committed data (pandoc, HTML and PDF), every number filled
  in by script.
- **README** with install, reproduction commands and result tables generated from data.

## Testing and verification

- Rust: unit tests (gates, channels, TABLE, REPEAT, seeding invariance, dense = sparse),
  property tests on random circuits.
- Python: pytest + Hypothesis: arithmetic exhaustive checks, synthesis distance bounds, Shor
  noiseless success against the exact formula, chemistry vs PySCF, QPE noiseless accuracy,
  calibration fits on recorded data, Aer cross-checks (skipped when Aer is absent).
- Every recorded dataset carries its seed, versions and command; `tools/reproduce.py` reruns it.

## Delivery

Local git repository `~/Desktop/vibes/ft-algorithms`, a commit per task on feature branches
merged to `master` when each phase is verified. Creating the GitHub repository and publishing the
site are left for Jasper (one command each), since that publishes under his name.

## Phases

1. `ftsim` engine + Python bindings + Aer cross-checks.
2. Compiler: synthesis, arithmetic, Shor; chemistry and QPE (noiseless correctness).
3. Calibration with `stabilizer-qec` (+ Stim/PyMatching cross-check, Choi experiment, composition).
4. Noise model, scheduling, experiments, error budgets.
5. Scaling estimates.
6. Site, WASM demo, report, README.
