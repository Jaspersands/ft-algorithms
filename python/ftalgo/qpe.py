"""Iterative phase estimation of a molecular Hamiltonian.

U = exp(−i(H − E_ref)τ) is approximated by r second-order Trotter steps; each controlled Pauli
rotation exp(−iθP) is a basis change, a CNOT parity ladder and a controlled Z rotation built
from two synthesized Rz's and two CNOTs:

    controlled-exp(−iθ Z_t) = exp(−iθZ_t/2) · CX(c,t) · exp(+iθZ_t/2) · CX(c,t).

Only uncontrolled rotations are synthesized, so their global phases stay global. The identity
term becomes a phase gate on the control. The eigenphase φ = −(E − E_ref)τ/2π (mod 1) is read
by ``phase.semiclassical_qpe`` with m bits; E = E_ref − 2π φ/τ with φ taken in [−½, ½).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from . import synth
from .chem.qubit import PauliSum, _letter
from .circuit import Circuit
from .phase import estimate, semiclassical_qpe

CHEMICAL_ACCURACY = 1.6e-3  # hartree (1 kcal/mol)


def _terms(H: PauliSum, e_ref: float) -> tuple[float, list[tuple[tuple[int, int], float]]]:
    c0 = H.terms.get((0, 0), 0.0) - e_ref
    rest = sorted(((k, float(v)) for k, v in H.terms.items() if k != (0, 0)), key=lambda kv: (kv[0][1], kv[0][0]))
    return c0, rest


SUZUKI_P = 1 / (4 - 4 ** (1 / 3))


def trotter_sequence(H: PauliSum, dt: float, e_ref: float, order: int = 2) -> tuple[float, list[tuple[tuple[int, int], float]]]:
    """One product-formula step: (identity angle, [(Pauli, θ)]) meaning Π exp(−iθP) in time
    order. Order 2 is the symmetric Strang step S2(dt); order 4 is Suzuki's
    S2(p·dt)² S2((1 − 4p)·dt) S2(p·dt)² with p = 1/(4 − 4^(1/3)). Adjacent rotations about the
    same Pauli are merged."""
    c0, rest = _terms(H, e_ref)

    def s2(t):
        return [(k, v * t / 2) for k, v in rest] + [(k, v * t / 2) for k, v in reversed(rest)]

    if order == 2:
        seq = s2(dt)
    elif order == 4:
        p = SUZUKI_P
        seq = s2(p * dt) + s2(p * dt) + s2((1 - 4 * p) * dt) + s2(p * dt) + s2(p * dt)
    else:
        raise ValueError("order must be 2 or 4")
    merged: list = []
    for k, th in seq:
        if merged and merged[-1][0] == k:
            merged[-1] = (k, merged[-1][1] + th)
        else:
            merged.append((k, th))
    return c0 * dt, merged


def trotter_unitary(H: PauliSum, dt: float, e_ref: float, order: int = 2) -> np.ndarray:
    """The exact matrix of one product-formula step (identity phase included), for checks."""
    phi0, seq = trotter_sequence(H, dt, e_ref, order)
    U = np.eye(1 << H.n, dtype=complex) * np.exp(-1j * phi0)
    for k, th in seq:
        P = PauliSum(H.n, {k: 1.0}).matrix()
        U = (np.cos(th) * np.eye(1 << H.n) - 1j * np.sin(th) * P) @ U
    return U


def controlled_pauli_rotation(c: Circuit, ctrl: int, sys: list[int], key: tuple[int, int], theta: float, eps: float) -> None:
    """controlled-exp(−iθP) for a Pauli string P on the system register."""
    n = len(sys)
    support = [q for q in range(n) if _letter(*key, q) != 0]
    letters = {q: _letter(*key, q) for q in support}
    for q in support:  # V: rotate each factor to Z
        if letters[q] == 1:
            c.h(sys[q])
        elif letters[q] == 2:
            c.sdg(sys[q])
            c.h(sys[q])
    t = sys[support[-1]]
    for q in support[:-1]:
        c.cx(sys[q], t)
    c.cx(ctrl, t)
    c.gates(synth.rz(-theta, eps), t)  # exp(+iθZ/2)
    c.cx(ctrl, t)
    c.gates(synth.rz(theta, eps), t)  # exp(−iθZ/2)
    for q in reversed(support[:-1]):
        c.cx(sys[q], t)
    for q in support:  # V†
        if letters[q] == 1:
            c.h(sys[q])
        elif letters[q] == 2:
            c.h(sys[q])
            c.s(sys[q])


def controlled_trotter_step(H: PauliSum, ctrl: int, sys: list[int], dt: float, e_ref: float, eps: float, num_qubits: int, order: int = 2) -> Circuit:
    phi0, seq = trotter_sequence(H, dt, e_ref, order)
    step = Circuit(num_qubits)
    for k, th in seq:
        controlled_pauli_rotation(step, ctrl, sys, k, th, eps)
    if abs(phi0) > 1e-15:
        step.gates(synth.rz(-phi0, eps), ctrl)  # phase e^{−iφ0} on the control's |1⟩
    return step


@dataclass
class QPEInfo:
    m: int
    tau: float
    steps: int
    order: int
    e_ref: float
    eps: float
    records: list[int]
    rotations_per_run: int

    def energies(self, rows: np.ndarray) -> np.ndarray:
        y = estimate(rows, self.records).astype(float)
        phi = y / (1 << self.m)
        phi = np.where(phi >= 0.5, phi - 1.0, phi)
        return self.e_ref - 2 * math.pi * phi / self.tau

    def resolution(self) -> float:
        return 2 * math.pi / (self.tau * (1 << self.m))


def qpe_circuit(H: PauliSum, hf_bits: int, *, m: int, tau: float, steps: int, e_ref: float, order: int = 2, eps: float | None = None, cutoff: int | None = None) -> tuple[Circuit, QPEInfo]:
    """Iterative QPE with m bits of U = exp(−i(H − e_ref)τ) ≈ (product-formula step(τ/steps))^steps.
    eps defaults to 0.05 / (number of synthesized rotations in a run)."""
    _, seq = trotter_sequence(H, tau / steps, e_ref, order)
    per_step = 2 * len(seq) + 1
    n_rot = ((1 << m) - 1) * steps * per_step + m
    eps = 0.05 / n_rot if eps is None else eps
    c = Circuit()
    ctrl = c.alloc(1, "ctrl")[0]
    sys = c.alloc(H.n, "sys")
    c.mark("prepare")
    for q in range(H.n):
        if hf_bits >> q & 1:
            c.x(sys[q])
    step = controlled_trotter_step(H, ctrl, sys, tau / steps, e_ref, eps, c.num_qubits, order)

    def power(circ: Circuit, k: int) -> None:
        circ.mark("evolve")
        circ.repeat((1 << k) * steps, step)

    recs = semiclassical_qpe(c, ctrl, m, m if cutoff is None else cutoff, power, eps)
    return c, QPEInfo(m, tau, steps, order, e_ref, eps, recs, n_rot)


def ideal_distribution(H: PauliSum, hf_bits: int, info: QPEInfo) -> np.ndarray:
    """Outcome probabilities of noiseless QPE with exact feedback for the Trotterized U (exact
    rotations): Σ_j |⟨v_j|HF⟩|² F(φ_j − y/2^m), with F the Fejér kernel."""
    step = trotter_unitary(H, info.tau / info.steps, info.e_ref, info.order)
    U = np.linalg.matrix_power(step, info.steps)
    w, v = np.linalg.eig(U)
    phases = np.angle(w) / (2 * math.pi) % 1.0
    ov = np.abs(np.linalg.solve(v, np.eye(len(w))[:, hf_bits])) ** 2
    M = 1 << info.m
    y = np.arange(M)
    p = np.zeros(M)
    for ph, o in zip(phases, ov):
        d = ph - y / M
        s = np.sin(np.pi * d)
        with np.errstate(divide="ignore", invalid="ignore"):
            f = np.where(np.abs(s) < 1e-12, 1.0, np.sin(np.pi * M * d) ** 2 / (M**2 * s**2))
        p += o * f
    return p / p.sum()
