# Phase 2: compiler and algorithms — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Logical circuits for Shor's algorithm (full modular exponentiation, Toffoli arithmetic,
semiclassical approximate QFT) and for iterative phase estimation of molecular Hamiltonians we
compute ourselves, all in Clifford + T + Toffoli, correct in the noiseless limit.

**Architecture:** A Python IR (`Circuit`) keeps REPEAT and TABLE structure and emits `ftsim`
text; it is also what Phase 4 schedules and decorates with noise. Rotations come from gridsynth
and are verified as matrices. Arithmetic is verified by an independent bit-level reversible
simulator. Chemistry is computed from closed-form s-Gaussian integrals and checked against PySCF.

**Tech Stack:** Python 3.13, numpy, scipy, pygridsynth 2.0 + mpmath, PySCF 2.14 (tests only),
the `ftalgo` engine from Phase 1.

## Global Constraints

- No "compiled Shor": the circuit depends only on N, a and classically computed a^(2^k) mod N;
  every multiplier is the generic circuit, even when its constant is 1.
- Every rotation is a verified Clifford+T sequence; every arithmetic block is verified on all
  inputs at small n; every chemistry number is checked against PySCF or a published value.
- Gate names are `ftsim` names. Qubit q of a register is bit q (little-endian).

## Files

```
python/ftalgo/circuit.py    IR: Circuit (alloc, gates, measure→record id, repeat, table, mark), to_text(), counts()
python/ftalgo/revsim.py     bit-level simulator of X/CX/CCX/SWAP/CSWAP circuits on Python ints
python/ftalgo/synth.py      rz(theta, eps) → Clifford+T gate list; verification; JSON cache in data/synth/
python/ftalgo/arith.py      Cuccaro add/sub, constant (un)loading, cc modular add, controlled modular multiply
python/ftalgo/phase.py      semiclassical (iterative) phase estimation with AQFT feedback by TABLE
python/ftalgo/shor.py       order-finding circuit, ideal distribution, classical post-processing, instances
python/ftalgo/chem/basis.py     STO-3G data; s-Gaussian integrals (S, T, V, ERI) via the Boys function
python/ftalgo/chem/scf.py       restricted Hartree–Fock, MO integrals, FCI by exact diagonalization
python/ftalgo/chem/qubit.py     Pauli algebra, Jordan–Wigner, Z2 tapering, molecules
python/ftalgo/qpe.py        Trotter steps, controlled evolution, energy from phase records
tests/test_circuit.py tests/test_synth.py tests/test_arith.py tests/test_shor.py tests/test_chem.py tests/test_qpe.py
```

## Interfaces

```python
# circuit.py
class Circuit:
    def __init__(self, num_qubits: int = 0)
    def alloc(self, n: int, name: str) -> list[int]
    def gate(self, name: str, *targets: int) -> None          # plus helpers x h s t cx cz ccx swap …
    def gates(self, seq: list[str], q: int) -> None            # a single-qubit sequence on q
    def measure(self, q: int, basis: str = "Z") -> int         # absolute record id
    def reset(self, q: int, basis: str = "Z") -> None
    def repeat(self, n: int, body: "Circuit") -> None          # body must not contain TABLE/measure-dependent refs
    def table(self, records: list[int], cases: list["Circuit"]) -> None
    def mark(self, label: str) -> None
    def extend(self, other: "Circuit") -> None
    def inverse(self) -> "Circuit"                             # gates only; reverses and daggers
    def to_text(self) -> str
    def counts(self) -> collections.Counter                    # REPEAT multiplied; TABLE = max over cases per gate kind
    num_qubits: int; num_records: int; ops: list

# revsim.py
def run_reversible(c: Circuit, state: int) -> int              # X CX CCX SWAP only (else ValueError)

# synth.py
def rz(theta: float, eps: float) -> list[str]                   # gates in time order; ≈ diag(1, e^{iθ}) up to phase
def phase_seq(theta, eps) == rz(theta, eps)
def matrix(seq: list[str]) -> np.ndarray
def distance_up_to_phase(U, V) -> float

# arith.py   (all registers little-endian lists of qubits)
def add(c, a, b, carry) / sub(...)                              # b ← b ± a  (mod 2^len(b)), len(a) == len(b)
def load_const(c, k, reg, ctrl=None)                            # reg ^= k (controlled by ctrl qubit if given)
def cc_mod_add(c, k, N, ctrls: list[int], b, tmp, carry, flag)  # b ← b + k·[all ctrls] mod N
def c_mod_mul(c, a, N, ctrl, x, b, tmp, carry, flag, anc)       # x ← a·x mod N if ctrl (b, tmp… returned clean)
class ShorLayout(n) -> qubit allocation used by shor.py

# phase.py
def semiclassical_qpe(c, ctrl, m, cutoff, emit_power, eps) -> list[int]   # records, least significant bit first

# shor.py
def order_finding(N, a, *, m=None, cutoff=None, eps=1e-6) -> (Circuit, info)
def ideal_distribution(N, a, m) -> np.ndarray                  # exact-QFT outcome probabilities
def outcome_value(records_row) -> int
def factors_from_outcome(y, m, N, a) -> tuple[int, int] | None
def choose_base(N) -> int
INSTANCES = [15, 21, 35, 51, 77, 143, 221, 247, 391, 899]

# chem
def molecule(name: str, R: float | None = None) -> Molecule     # "H2", "HeH+", "H3+", "H4"
def qubit_hamiltonian(mol) -> (PauliSum, hf_bits, info)        # tapered
# qpe.py
def qpe_circuit(ham, hf_bits, *, bits, tau, steps, eps) -> (Circuit, info)
def energy_from_records(rows, info) -> np.ndarray
```

## Tasks

### Task 1: Circuit IR + reversible simulator
- [ ] Tests: text round trip parses in `ftalgo.Program`; record offsets in TABLE are right
  (measure 3 qubits, table on the first and third); counts multiply through REPEAT; inverse of a
  random Clifford+T circuit composed with it is identity (via engine statevector); revsim on a
  random X/CX/CCX/SWAP circuit equals ftsim on every basis input of 5 qubits.
- [ ] Implement, pass, commit.

### Task 2: Rotation synthesis
- [ ] Tests: exact multiples of π/4 return ≤ 3 gates without calling gridsynth; for 50 random θ
  at ε ∈ {1e-3, 1e-6} the sequence matrix is within ε of diag(1, e^{iθ}) up to phase; sequences
  contain only H, S, T, X, Y, Z (S_DAG/T_DAG allowed); cache round trip; engine statevector of
  `H 0` + sequence equals the target state up to phase.
- [ ] Implement, pass, commit.

### Task 3: Arithmetic
- [ ] Tests (revsim, exhaustive): add/sub for every (a, b) at widths 1..5; load_const; cc_mod_add
  for N ∈ {3, 5, 7, 11, 13, 15, 21} (all b < N, all control values, all k < N), ancillas clean;
  c_mod_mul for N ∈ {15, 21, 35} with every x < N, both control values, every invertible a;
  gate counts reported (Toffoli = CCX).
- [ ] Implement, pass, commit.

### Task 4: Shor
- [ ] Tests: `ideal_distribution(15, 7, 8)` peaks at multiples of 64; noiseless ftsim sampling
  with exact QFT (cutoff = m) matches ideal_distribution by chi-square for N = 15, 21; AQFT
  cutoff default gives success probability within 0.01 of exact; post-processing finds factors
  from the ideal peaks; `choose_base` returns a valid base for every instance.
- [ ] Implement, pass, commit.

### Task 5: Chemistry
- [ ] Tests: STO-3G H2 at 1.4 bohr: HF −1.1167 and FCI −1.1373 hartree (Szabo & Ostlund) to
  1e-4; H2, HeH+, H3+, H4 integrals, HF and FCI energies equal PySCF to 1e-8; JW Hamiltonian's
  ground energy in the HF particle sector equals FCI; tapered Hamiltonian's ground energy equals
  FCI and the tapered HF basis state has energy E_HF.
- [ ] Implement, pass, commit.

### Task 6: Phase estimation
- [ ] Tests: on a 1-qubit Hamiltonian H = 0.3 Z with exact eigenstate input, noiseless iterative
  QPE returns the eigenvalue's bits deterministically when representable; Trotter step matrix
  error matches the second-order bound; H2 noiseless QPE with chosen (bits, tau, steps, eps)
  lands within 1.6 mHa of FCI in ≥ 90% of shots (sampled with ftsim).
- [ ] Implement, pass, commit.
