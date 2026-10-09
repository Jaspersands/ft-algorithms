"""Modern fault-tolerant reversible arithmetic for Shor's algorithm.

Improvements over textbook Cuccaro arithmetic:
1. Gidney (2018) measurement-based carry uncomputation (arXiv:1709.06648):
   Reduces adder Toffolis from 2w to w - 1 for w-bit addition mod 2^w by replacing
   unitary carry uncomputation with X-basis measurements, reset, and classically
   controlled CZ fixups. Subtraction is performed via b <- ~(~b + a) mod 2^w with zero
   additional Toffolis.
2. Gidney (2019) windowed arithmetic (arXiv:1905.09749):
   Groups multiplier bits into k-bit windows (default k=2), performing 1 modular addition
   per window using QROM table lookups rather than k individual modular additions.
   This cuts the number of modular additions by roughly a factor of k.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from fractions import Fraction

import numpy as np

from .circuit import Circuit
from .phase import estimate, semiclassical_qpe
from .shor import ShorInfo, choose_base, default_cutoff


def gidney_add(c: Circuit, a: list[int], b: list[int], carries: list[int]) -> None:
    """b <- a + b mod 2^w using Gidney 2018 measurement-based carry uncomputation.

    Requires w = len(a) = len(b), and len(carries) == w - 1.
    Uses w - 1 Toffoli gates total (compared to 2w for Cuccaro).
    carries are clean ancillas and are restored clean to |0>.
    """
    w = len(a)
    if len(b) != w:
        raise ValueError("gidney_add: registers differ in length")
    if len(carries) != max(w - 1, 0):
        raise ValueError(f"gidney_add: expected {max(w - 1, 0)} carries, got {len(carries)}")
    if w == 1:
        c.cx(a[0], b[0])
        return

    # Forward carry propagation
    c.ccx(a[0], b[0], carries[0])
    for k in range(1, w - 1):
        c.cx(carries[k - 1], a[k])
        c.cx(carries[k - 1], b[k])
        c.ccx(a[k], b[k], carries[k])
        c.cx(carries[k - 1], carries[k])

    # Top bit (k = w - 1): carry into b[w - 1]
    c.cx(carries[w - 2], b[w - 1])
    c.cx(a[w - 1], b[w - 1])

    # Backward uncomputation of carries and sum computation
    for k in range(w - 2, 0, -1):
        c.cx(carries[k - 1], carries[k])
        # Measurement-based temporary AND uncomputation
        c.h(carries[k])
        m = c.measure(carries[k])
        case0 = Circuit(c.num_qubits)
        case1 = Circuit(c.num_qubits)
        case1.cz(a[k], b[k])
        c.table([m], [case0, case1])
        c.reset(carries[k])
        c.cx(carries[k - 1], a[k])
        c.cx(a[k], b[k])

    # Bit 0 uncompute
    c.h(carries[0])
    m0 = c.measure(carries[0])
    case0 = Circuit(c.num_qubits)
    case1 = Circuit(c.num_qubits)
    case1.cz(a[0], b[0])
    c.table([m0], [case0, case1])
    c.reset(carries[0])
    c.cx(a[0], b[0])


def gidney_sub(c: Circuit, a: list[int], b: list[int], carries: list[int]) -> None:
    """b <- b - a mod 2^w via b <- ~(~b + a) mod 2^w using bitwise NOT."""
    for q in b:
        c.x(q)
    gidney_add(c, a, b, carries)
    for q in b:
        c.x(q)


def load_const(c: Circuit, k: int, reg: list[int], ctrl: int | None = None) -> None:
    """reg ^= k, or reg ^= k * ctrl."""
    if k < 0 or k >> len(reg):
        raise ValueError(f"constant {k} does not fit in {len(reg)} bits")
    for i, q in enumerate(reg):
        if k >> i & 1:
            if ctrl is None:
                c.x(q)
            else:
                c.cx(ctrl, q)


def load_window_table(
    c: Circuit,
    C: list[int],
    tmp: list[int],
    ctrl: int,
    x_bits: list[int],
    anc: list[int],
) -> None:
    """Loads classical constant C[v] into tmp conditioned on ctrl and x_bits.

    Supports len(x_bits) in (1, 2).
    """
    if len(x_bits) == 1:
        c1_val = C[1]
        c.ccx(ctrl, x_bits[0], anc[0])
        for j, q in enumerate(tmp):
            if (c1_val >> j) & 1:
                c.cx(anc[0], q)
    elif len(x_bits) == 2:
        g1, g2, g3 = anc[0], anc[1], anc[2]
        c.ccx(ctrl, x_bits[0], g1)
        c.ccx(ctrl, x_bits[1], g2)
        c.ccx(g1, x_bits[1], g3)
        c1_val, c2_val, c3_val = C[1], C[2], C[3]
        for j, q in enumerate(tmp):
            b1 = (c1_val >> j) & 1
            b2 = (c2_val >> j) & 1
            b3 = ((c3_val >> j) & 1) ^ b1 ^ b2
            if b1:
                c.cx(g1, q)
            if b2:
                c.cx(g2, q)
            if b3:
                c.cx(g3, q)
    else:
        raise ValueError(f"Window size {len(x_bits)} not supported (use 1 or 2)")


def unload_window_table(
    c: Circuit,
    C: list[int],
    tmp: list[int],
    ctrl: int,
    x_bits: list[int],
    anc: list[int],
) -> None:
    """Unloads table C from tmp and clears ancillas."""
    if len(x_bits) == 1:
        c1_val = C[1]
        for j, q in enumerate(tmp):
            if (c1_val >> j) & 1:
                c.cx(anc[0], q)
        c.ccx(ctrl, x_bits[0], anc[0])
    elif len(x_bits) == 2:
        g1, g2, g3 = anc[0], anc[1], anc[2]
        c1_val, c2_val, c3_val = C[1], C[2], C[3]
        for j, q in enumerate(tmp):
            b1 = (c1_val >> j) & 1
            b2 = (c2_val >> j) & 1
            b3 = ((c3_val >> j) & 1) ^ b1 ^ b2
            if b3:
                c.cx(g3, q)
            if b2:
                c.cx(g2, q)
            if b1:
                c.cx(g1, q)
        c.ccx(g1, x_bits[1], g3)
        c.ccx(ctrl, x_bits[1], g2)
        c.ccx(ctrl, x_bits[0], g1)
    else:
        raise ValueError(f"Window size {len(x_bits)} not supported (use 1 or 2)")


def window_mod_add(
    c: Circuit,
    C: list[int],
    N: int,
    ctrl: int,
    x_bits: list[int],
    b: list[int],
    tmp: list[int],
    carries: list[int],
    flag: int,
    anc: list[int],
    subtract: bool = False,
) -> None:
    """b <- b +/- C(x_bits) mod N using Beauregard reduction with Gidney adders."""
    top = b[-1]
    if not subtract:
        load_window_table(c, C, tmp, ctrl, x_bits, anc)
        gidney_add(c, tmp, b, carries)
        unload_window_table(c, C, tmp, ctrl, x_bits, anc)

        load_const(c, N, tmp)
        gidney_sub(c, tmp, b, carries)
        load_const(c, N, tmp)

        c.cx(top, flag)

        load_const(c, N, tmp, flag)
        gidney_add(c, tmp, b, carries)
        load_const(c, N, tmp, flag)

        load_window_table(c, C, tmp, ctrl, x_bits, anc)
        gidney_sub(c, tmp, b, carries)
        unload_window_table(c, C, tmp, ctrl, x_bits, anc)

        c.x(top)
        c.cx(top, flag)
        c.x(top)

        load_window_table(c, C, tmp, ctrl, x_bits, anc)
        gidney_add(c, tmp, b, carries)
        unload_window_table(c, C, tmp, ctrl, x_bits, anc)
    else:
        load_window_table(c, C, tmp, ctrl, x_bits, anc)
        gidney_sub(c, tmp, b, carries)
        unload_window_table(c, C, tmp, ctrl, x_bits, anc)

        c.x(top)
        c.cx(top, flag)
        c.x(top)

        load_window_table(c, C, tmp, ctrl, x_bits, anc)
        gidney_add(c, tmp, b, carries)
        unload_window_table(c, C, tmp, ctrl, x_bits, anc)

        load_const(c, N, tmp, flag)
        gidney_sub(c, tmp, b, carries)
        load_const(c, N, tmp, flag)

        c.cx(top, flag)

        load_const(c, N, tmp)
        gidney_add(c, tmp, b, carries)
        load_const(c, N, tmp)

        load_window_table(c, C, tmp, ctrl, x_bits, anc)
        gidney_sub(c, tmp, b, carries)
        unload_window_table(c, C, tmp, ctrl, x_bits, anc)


def windowed_c_mod_mul(
    c: Circuit,
    a: int,
    N: int,
    ctrl: int,
    x: list[int],
    b: list[int],
    tmp: list[int],
    carries: list[int],
    flag: int,
    anc: list[int],
    k: int = 2,
) -> None:
    """x <- a * x mod N if ctrl (x < N, gcd(a, N) = 1) using windowed arithmetic."""
    n = len(x)
    if len(b) != n + 1:
        raise ValueError("windowed_c_mod_mul: b must have len(x) + 1 bits")
    a %= N
    a_inv = pow(a, -1, N)

    for i in range(0, n, k):
        x_bits = x[i : min(i + k, n)]
        K = (a * (1 << i)) % N
        num_vals = 1 << len(x_bits)
        C = [(v * K) % N for v in range(num_vals)]
        window_mod_add(c, C, N, ctrl, x_bits, b, tmp, carries, flag, anc, subtract=False)

    for i in range(n):
        c.cswap(ctrl, x[i], b[i])

    for i in range(0, n, k):
        x_bits = x[i : min(i + k, n)]
        K_inv = (a_inv * (1 << i)) % N
        num_vals = 1 << len(x_bits)
        C_inv = [(v * K_inv) % N for v in range(num_vals)]
        window_mod_add(c, C_inv, N, ctrl, x_bits, b, tmp, carries, flag, anc, subtract=True)


@dataclass
class ModernShorLayout:
    """Registers for modern Shor order finding: 4n + 7 qubits for an n-bit modulus."""

    n: int
    ctrl: int
    x: list[int]
    b: list[int]
    tmp: list[int]
    carries: list[int]
    flag: int
    anc: list[int]

    @staticmethod
    def allocate(c: Circuit, n: int) -> "ModernShorLayout":
        ctrl = c.alloc(1, "ctrl")[0]
        x = c.alloc(n, "x")
        b = c.alloc(n + 1, "b")
        tmp = c.alloc(n + 1, "tmp")
        carries = c.alloc(n, "carries")  # w - 1 carries for w = n + 1
        flag = c.alloc(1, "flag")[0]
        anc = c.alloc(3, "anc")
        return ModernShorLayout(n, ctrl, x, b, tmp, carries, flag, anc)

    def mod_mul(self, c: Circuit, a: int, N: int, k: int = 2) -> None:
        windowed_c_mod_mul(
            c, a, N, self.ctrl, self.x, self.b, self.tmp, self.carries, self.flag, self.anc, k=k
        )

    @property
    def ancillas(self) -> list[int]:
        return [*self.b, *self.tmp, *self.carries, self.flag, *self.anc]


def modern_toffolis(n: int, k: int = 2) -> int:
    """Total Toffoli count of order_finding_modern for an n-bit modulus."""
    # Number of windows
    num_windows = (n + k - 1) // k
    # Toffolis per window in window_mod_add:
    # 5 additions of width n + 1 (each n Toffolis) = 5n
    # 3 load/unloads (each 3 Toffolis for k=2, 1 for k=1): 3 * 2 * 3 = 18 for k=2
    tof_per_window_2 = 5 * n + 18
    tof_per_window_1 = 5 * n + 6

    tof_mul = 0
    for i in range(0, n, k):
        w_size = min(k, n - i)
        w_cost = tof_per_window_2 if w_size == 2 else tof_per_window_1
        tof_mul += 2 * w_cost  # forward and backward
    tof_mul += n  # n CSWAPs

    m = 2 * n  # 2n multiplications in semiclassical QPE
    return m * tof_mul


def order_finding_modern(
    N: int,
    a: int | None = None,
    *,
    m: int | None = None,
    cutoff: int | None = None,
    eps: float | None = None,
    k: int = 2,
) -> tuple[Circuit, ShorInfo]:
    """The modern order-finding circuit for a mod N with an m-bit phase estimate."""
    if N < 3 or N % 2 == 0:
        raise ValueError("N must be odd and >= 3")
    a = choose_base(N) if a is None else a
    n = N.bit_length()
    m = 2 * n if m is None else m
    cutoff = default_cutoff(m) if cutoff is None else cutoff
    eps = 1e-3 / m if eps is None else eps
    c = Circuit()
    lay = ModernShorLayout.allocate(c, n)
    c.mark("prepare")
    c.x(lay.x[0])

    def power(circ: Circuit, pwr_idx: int) -> None:
        circ.mark("multiply")
        lay.mod_mul(circ, pow(a, 1 << pwr_idx, N), N, k=k)

    recs = semiclassical_qpe(c, lay.ctrl, m, cutoff, power, eps)
    info = ShorInfo(N, a, n, m, cutoff, eps, recs, lay)  # type: ignore[arg-type]
    info.counts = dict(c.counts())
    return c, info
