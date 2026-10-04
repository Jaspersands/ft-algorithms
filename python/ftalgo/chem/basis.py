"""STO-3G s-type Gaussian basis functions and their integrals in closed form.

For s Gaussians g(r) = (2α/π)^{3/4} e^{−α|r−A|²} every one- and two-electron integral has a
closed form involving only the Boys function F₀(t) = ½√(π/t) erf(√t) (Szabo & Ostlund,
Appendix A). That covers H and He, so H2, HeH+, H3+ and H4 need nothing else.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.special import erf

# STO-3G contraction coefficients of a 1s Slater function (shared by H and He).
STO3G_COEFFS = (0.15432897, 0.53532814, 0.44463454)
# Exponents for ζ = 1 (Szabo & Ostlund Table 3.8 / Hehre, Stewart and Pople 1969).
STO3G_ZETA1 = (2.227660, 0.405771, 0.109818)
# Basis Set Exchange STO-3G exponents (what PySCF's "sto-3g" uses).
STO3G_BSE = {
    "H": (3.42525091, 0.62391373, 0.16885540),
    "He": (6.36242139, 1.15892300, 0.31364979),
}
CHARGE = {"H": 1, "He": 2}


def sto3g_exponents(element: str, zeta: float | None = None) -> tuple[float, ...]:
    """BSE exponents, or the ζ = 1 set scaled by ζ² (Szabo & Ostlund use ζ_H = 1.24, ζ_He = 2.0925)."""
    if zeta is None:
        return STO3G_BSE[element]
    return tuple(a * zeta**2 for a in STO3G_ZETA1)


@dataclass(frozen=True)
class SFunction:
    """A contracted s function: Σ_k d_k N(α_k) e^{−α_k |r − center|²}."""

    center: tuple[float, float, float]
    alphas: tuple[float, ...]
    coeffs: tuple[float, ...]

    def primitives(self):
        for a, d in zip(self.alphas, self.coeffs):
            yield a, d * (2 * a / math.pi) ** 0.75


def boys0(t: float) -> float:
    if t < 1e-12:
        return 1.0 - t / 3.0
    return 0.5 * math.sqrt(math.pi / t) * float(erf(math.sqrt(t)))


def _d2(A, B) -> float:
    return sum((a - b) ** 2 for a, b in zip(A, B))


def _gp(a, A, b, B):
    p = a + b
    return p, tuple((a * x + b * y) / p for x, y in zip(A, B))


def overlap(f: SFunction, g: SFunction) -> float:
    s = 0.0
    r2 = _d2(f.center, g.center)
    for a, ca in f.primitives():
        for b, cb in g.primitives():
            p = a + b
            s += ca * cb * (math.pi / p) ** 1.5 * math.exp(-a * b / p * r2)
    return s


def kinetic(f: SFunction, g: SFunction) -> float:
    s = 0.0
    r2 = _d2(f.center, g.center)
    for a, ca in f.primitives():
        for b, cb in g.primitives():
            p = a + b
            mu = a * b / p
            s += ca * cb * mu * (3 - 2 * mu * r2) * (math.pi / p) ** 1.5 * math.exp(-mu * r2)
    return s


def nuclear(f: SFunction, g: SFunction, C, Z: float) -> float:
    s = 0.0
    r2 = _d2(f.center, g.center)
    for a, ca in f.primitives():
        for b, cb in g.primitives():
            p, P = _gp(a, f.center, b, g.center)
            s += ca * cb * (-2 * math.pi / p) * Z * math.exp(-a * b / p * r2) * boys0(p * _d2(P, C))
    return s


def eri(f: SFunction, g: SFunction, h: SFunction, k: SFunction) -> float:
    """(fg|hk) in chemists' notation."""
    s = 0.0
    rab = _d2(f.center, g.center)
    rcd = _d2(h.center, k.center)
    for a, ca in f.primitives():
        for b, cb in g.primitives():
            p, P = _gp(a, f.center, b, g.center)
            eab = math.exp(-a * b / p * rab)
            for c, cc in h.primitives():
                for d, cd in k.primitives():
                    q, Q = _gp(c, h.center, d, k.center)
                    ecd = math.exp(-c * d / q * rcd)
                    s += (
                        ca * cb * cc * cd
                        * 2 * math.pi**2.5 / (p * q * math.sqrt(p + q))
                        * eab * ecd * boys0(p * q / (p + q) * _d2(P, Q))
                    )
    return s


@dataclass
class AOIntegrals:
    S: np.ndarray
    T: np.ndarray
    V: np.ndarray
    eri: np.ndarray  # (ij|kl)
    e_nuc: float

    @property
    def H(self) -> np.ndarray:
        return self.T + self.V


def ao_integrals(atoms: list[tuple[str, tuple[float, float, float]]], basis: list[SFunction]) -> AOIntegrals:
    n = len(basis)
    S = np.array([[overlap(basis[i], basis[j]) for j in range(n)] for i in range(n)])
    T = np.array([[kinetic(basis[i], basis[j]) for j in range(n)] for i in range(n)])
    V = np.zeros((n, n))
    for el, C in atoms:
        V += np.array([[nuclear(basis[i], basis[j], C, CHARGE[el]) for j in range(n)] for i in range(n)])
    G = np.zeros((n, n, n, n))
    for i in range(n):
        for j in range(i + 1):
            for k in range(n):
                for l in range(k + 1):
                    if i * (i + 1) // 2 + j < k * (k + 1) // 2 + l:
                        continue
                    v = eri(basis[i], basis[j], basis[k], basis[l])
                    for (a, b, c, d) in ((i, j, k, l), (j, i, k, l), (i, j, l, k), (j, i, l, k)):
                        G[a, b, c, d] = v
                        G[c, d, a, b] = v
    e_nuc = 0.0
    for x in range(len(atoms)):
        for y in range(x):
            e_nuc += CHARGE[atoms[x][0]] * CHARGE[atoms[y][0]] / math.sqrt(_d2(atoms[x][1], atoms[y][1]))
    return AOIntegrals(S, T, V, G, e_nuc)
