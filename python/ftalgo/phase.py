"""Semiclassical (iterative) phase estimation with one recycled control qubit.

Griffiths and Niu (1996): the inverse QFT of phase estimation becomes single-qubit rotations
conditioned on earlier measurement results. Round j (j = 0 … m − 1) applies the controlled
power U^(2^(m−1−j)), which puts the phase 2π·0.φ_{m−j}φ_{m−j+1}…φ_m on |1⟩ of the control; the
rotation by −π Σ_{i≥1} r_{j−i} / 2^i removes the already measured bits, leaving ±1, which H and a
Z measurement read as r_j = φ_{m−j}. The estimate is y = Σ_j r_j 2^j ≈ φ·2^m.

Coppersmith's approximate QFT keeps only the terms i < cutoff, so each correction depends on at
most cutoff − 1 earlier records: a TABLE with ≤ 2^(cutoff−1) cases, each a synthesized rotation.
"""

from __future__ import annotations

import math
from typing import Callable

from . import synth
from .circuit import Circuit


def correction_angle(bits: list[int]) -> float:
    """−π Σ_i bits[i−1] / 2^i, where bits[0] is the most recent record."""
    return -math.pi * sum(b / (1 << (i + 1)) for i, b in enumerate(bits))


def semiclassical_qpe(
    c: Circuit,
    ctrl: int,
    m: int,
    cutoff: int,
    emit_power: Callable[[Circuit, int], None],
    eps: float,
) -> list[int]:
    """Appends m rounds to c; emit_power(c, k) must append the controlled-U^(2^k) with control
    `ctrl`. Returns the record ids, least significant bit of the estimate first."""
    if cutoff < 1:
        raise ValueError("cutoff ≥ 1")
    recs: list[int] = []
    for j in range(m):
        k = m - 1 - j
        c.mark("prepare")
        c.reset(ctrl)
        c.h(ctrl)
        emit_power(c, k)
        prev = [recs[j - i] for i in range(1, min(j, cutoff - 1) + 1)]  # most recent first
        if prev:
            c.mark("rotation")
            cases = []
            for idx in range(1 << len(prev)):
                bits = [(idx >> t) & 1 for t in range(len(prev))]
                case = Circuit(c.num_qubits)
                case.gates(synth.rz(correction_angle(bits), eps), ctrl)
                cases.append(case)
            c.table(prev, cases)
        c.mark("measure")
        c.h(ctrl)
        recs.append(c.measure(ctrl))
    synth.save_cache()
    return recs


def estimate(rows, recs: list[int]):
    """y = Σ_j r_j 2^j for each shot (rows: shots × records bool array, recs: record ids)."""
    import numpy as np

    weights = np.array([1 << j for j in range(len(recs))], dtype=object if len(recs) > 62 else np.int64)
    return (rows[:, recs].astype(np.int64) * weights).sum(axis=1)
