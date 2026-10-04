"""Shor's algorithm: order finding with the full modular exponentiation.

The circuit is generic: it depends only on N, the base a and the classically computed constants
a^(2^k) mod N. Every multiplication is the full controlled modular multiplier of ``arith``,
even when its constant happens to be 1, so nothing in the circuit uses the order r (the
"compiled Shor" shortcut Smolin, Smith and Vargo warned about). The order is used only to
compute the ideal output distribution against which results are judged.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from fractions import Fraction

import numpy as np

from .arith import ShorLayout
from .circuit import Circuit
from .phase import estimate, semiclassical_qpe

INSTANCES = [15, 21, 35, 51, 77, 143, 221, 247, 391, 899]


def multiplicative_order(a: int, N: int) -> int:
    if math.gcd(a, N) != 1:
        raise ValueError("a and N share a factor")
    r, x = 1, a % N
    while x != 1:
        x = x * a % N
        r += 1
    return r


def is_good_base(a: int, N: int) -> bool:
    """Whether the order r of a is even with a^(r/2) ≠ −1 mod N (so the order yields a factor)."""
    if math.gcd(a, N) != 1:
        return False
    r = multiplicative_order(a, N)
    return r % 2 == 0 and pow(a, r // 2, N) != N - 1


def choose_base(N: int) -> int:
    """The base used for N: 2, the smallest base, when it is coprime to N. (Whether it is a good
    base is reported separately; it is not used to choose.)"""
    for a in range(2, N):
        if math.gcd(a, N) == 1:
            return a
    raise ValueError(f"no base for {N}")


def default_cutoff(m: int) -> int:
    """AQFT: keep rotations down to angle π/2^(cutoff−1); ⌈log2 m⌉ + 3 bits keeps the output
    distribution within ~1% of the exact QFT's (checked in tests)."""
    return min(m, math.ceil(math.log2(m)) + 3)


@dataclass
class ShorInfo:
    N: int
    a: int
    n: int
    m: int
    cutoff: int
    eps: float
    records: list[int]
    layout: ShorLayout
    counts: dict = field(default_factory=dict)

    def outcomes(self, rows: np.ndarray) -> np.ndarray:
        return estimate(rows, self.records)


def order_finding(N: int, a: int | None = None, *, m: int | None = None, cutoff: int | None = None, eps: float | None = None) -> tuple[Circuit, ShorInfo]:
    """The order-finding circuit for a mod N with an m-bit phase estimate (default 2n)."""
    if N < 3 or N % 2 == 0:
        raise ValueError("N must be odd and ≥ 3")
    a = choose_base(N) if a is None else a
    n = N.bit_length()
    m = 2 * n if m is None else m
    cutoff = default_cutoff(m) if cutoff is None else cutoff
    eps = 1e-3 / m if eps is None else eps
    c = Circuit()
    lay = ShorLayout.allocate(c, n)
    c.mark("prepare")
    c.x(lay.x[0])

    def power(circ: Circuit, k: int) -> None:
        circ.mark("multiply")
        lay.mod_mul(circ, pow(a, 1 << k, N), N)

    recs = semiclassical_qpe(c, lay.ctrl, m, cutoff, power, eps)
    info = ShorInfo(N, a, n, m, cutoff, eps, recs, lay)
    info.counts = dict(c.counts())
    return c, info


def ideal_distribution(N: int, a: int, m: int) -> np.ndarray:
    """P(y) for y = 0 … 2^m − 1 under exact phase estimation: the register starts in
    |1⟩ = r^(−1/2) Σ_s |u_s⟩, eigenphases s/r, so P(y) = (1/r) Σ_s F(s/r − y/2^m) with the Fejér
    kernel F(δ) = sin²(π 2^m δ) / (2^(2m) sin²(π δ)) (F = 1 at integer δ)."""
    r = multiplicative_order(a, N)
    M = 1 << m
    y = np.arange(M)
    p = np.zeros(M)
    for s in range(r):
        d = s / r - y / M
        num = np.sin(np.pi * M * d) ** 2
        den = (M**2) * np.sin(np.pi * d) ** 2
        with np.errstate(divide="ignore", invalid="ignore"):
            f = np.where(np.abs(np.sin(np.pi * d)) < 1e-12, 1.0, num / den)
        p += f / r
    return p


def convergent_denominators(y: int, M: int, limit: int) -> list[int]:
    """Denominators < limit of the continued-fraction convergents of y/M."""
    out = []
    if y == 0:
        return out
    x = Fraction(y, M)
    h0, h1, k0, k1 = 0, 1, 1, 0
    while True:
        q = x.numerator // x.denominator
        h0, h1 = h1, q * h1 + h0
        k0, k1 = k1, q * k1 + k0
        if k1 >= limit:
            break
        out.append(k1)
        frac = x - q
        if frac == 0:
            break
        x = 1 / frac
    return out


def factors_from_outcome(y: int, m: int, N: int, a: int) -> tuple[int, int] | None:
    """Classical post-processing of one run: candidate orders are convergent denominators of
    y/2^m and their small multiples below N; the first q with a^q ≡ 1 is taken as the order, and
    gcd(a^(q/2) ± 1, N) gives the factors when q is even and a^(q/2) ≠ −1."""
    for q in convergent_denominators(y, 1 << m, N):
        t = 1
        while q * t < N:
            r = q * t
            if pow(a, r, N) == 1:
                if r % 2 == 0:
                    h = pow(a, r // 2, N)
                    if h != N - 1:
                        f = math.gcd(h - 1, N)
                        if 1 < f < N:
                            return (min(f, N // f), max(f, N // f))
                return None
            t += 1
    return None


def success_table(N: int, a: int, m: int) -> np.ndarray:
    """success[y] = whether outcome y yields the factors."""
    return np.array([factors_from_outcome(y, m, N, a) is not None for y in range(1 << m)])


def ideal_success(N: int, a: int, m: int) -> float:
    return float((ideal_distribution(N, a, m) * success_table(N, a, m)).sum())


def peak_table(N: int, a: int, m: int) -> np.ndarray:
    """peak[y] = whether y is one of the r ideal peaks round(s·2^m/r), s = 0 … r − 1. The share of
    shots landing on a peak is the primary score: its random baseline is only r/2^m, whereas for
    small N almost any outcome yields the factors after post-processing (Smolin, Smith & Vargo
    2013), which makes factoring success alone a weak test. The order is used to score, never
    in the circuit."""
    r = multiplicative_order(a, N)
    M = 1 << m
    t = np.zeros(M, dtype=bool)
    for s in range(r):
        t[round(s * M / r) % M] = True
    return t


def ideal_peak_probability(N: int, a: int, m: int) -> float:
    return float((ideal_distribution(N, a, m) * peak_table(N, a, m)).sum())
