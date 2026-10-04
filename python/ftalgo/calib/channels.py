"""Pauli channels through their fidelities.

A distribution P over flip patterns (or Pauli errors) of n bits has fidelities
f(s) = Σ_x P(x)(−1)^{s·x}. Independent channels compose by multiplying fidelities, so

- per-round fidelity from two lengths: f_round = (f(T2)/f(T1))^{1/(T2−T1)}
  (preparation and readout cancel), and
- an operation's own channel: f_op = f_measured / f_baseline.

``probabilities`` inverts (Walsh–Hadamard) and clips the tiny negative values statistical
noise can leave.
"""

from __future__ import annotations

import numpy as np


def _signs(n: int) -> np.ndarray:
    k = np.arange(1 << n)
    return np.array([[(-1) ** bin(s & x).count("1") for x in k] for s in k], dtype=float)


def fidelities(p: np.ndarray) -> np.ndarray:
    n = int(np.log2(len(p)))
    return _signs(n) @ p


def probabilities(f: np.ndarray, clip: bool = True) -> np.ndarray:
    n = int(np.log2(len(f)))
    p = _signs(n) @ f / (1 << n)
    if clip:
        p = np.clip(p, 0, None)
        p[0] = 0
        p[0] = 1 - p.sum()
    return p


def per_round(f1: np.ndarray, f2: np.ndarray, t1: int, t2: int) -> np.ndarray:
    r = np.clip(f2 / f1, 1e-300, None)
    return r ** (1.0 / (t2 - t1))


def deconvolve(f_measured: np.ndarray, f_baseline: np.ndarray) -> np.ndarray:
    return f_measured / f_baseline


def product_fidelities(*per_bit_flip: float) -> np.ndarray:
    """Fidelities of independent single-bit flips with probabilities q_j (bit j)."""
    n = len(per_bit_flip)
    f = np.ones(1 << n)
    for s in range(1 << n):
        for j, q in enumerate(per_bit_flip):
            if s >> j & 1:
                f[s] *= 1 - 2 * q
    return f


def remap(p: np.ndarray, matrix: list[list[int]]) -> np.ndarray:
    """Distribution of y = M·x over GF(2) for x ~ p (rows of M are output bits as masks of x)."""
    out = np.zeros(1 << len(matrix))
    for x, px in enumerate(p):
        y = 0
        for i, row in enumerate(matrix):
            y |= (bin(row & x).count("1") & 1) << i
        out[y] += px
    return out


def pauli_channel_1(px: float, py: float, pz: float) -> tuple[float, float, float]:
    return (px, py, pz)


def pauli_channel_2(xpart: np.ndarray, zpart: np.ndarray) -> list[float]:
    """The 15 PAULI_CHANNEL_2 probabilities (Stim order: index 4a + b − 1, a on the first qubit)
    from independent X-part and Z-part distributions over 2 bits (bit 0 = first qubit)."""
    out = []
    for k in range(1, 16):
        a, b = k // 4, k % 4  # 0 I, 1 X, 2 Y, 3 Z
        xb = (a in (1, 2)) | ((b in (1, 2)) << 1)
        zb = (a in (2, 3)) | ((b in (2, 3)) << 1)
        out.append(float(xpart[xb] * zpart[zb]))
    return out
