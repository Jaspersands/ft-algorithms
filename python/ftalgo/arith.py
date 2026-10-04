"""Reversible arithmetic from X, CX and Toffoli gates (no rotations), for Shor's algorithm.

Registers are little-endian lists of qubits.

- ``add``/``sub``: Cuccaro–Draper–Kutin–Moulton ripple-carry adder (2004), b ← b ± a mod 2^w,
  one clean carry ancilla, 2w Toffolis.
- ``cc_mod_add``: b ← b + k·[controls] mod N for a classical constant k, following the
  Vedral–Barenco–Ekert / Beauregard structure (add k, subtract N, copy the sign to a flag, add N
  back if negative, subtract k to clear the flag, add k), with k and N loaded into a scratch
  register by CX gates.
- ``c_mod_mul``: x ← a·x mod N if ctrl, by accumulating a·x into a zeroed register, swapping it
  with x under control, and uncomputing the old x by subtracting a⁻¹·(a·x).

Every block is checked on all inputs by ``tests/test_arith.py`` with the independent bit-level
simulator ``ftalgo.revsim``.
"""

from __future__ import annotations

from dataclasses import dataclass

from .circuit import Circuit


def _maj(c: Circuit, x: int, y: int, z: int) -> None:
    c.cx(z, y)
    c.cx(z, x)
    c.ccx(x, y, z)


def _uma(c: Circuit, x: int, y: int, z: int) -> None:
    c.ccx(x, y, z)
    c.cx(z, x)
    c.cx(x, y)


def add(c: Circuit, a: list[int], b: list[int], carry: int) -> None:
    """b ← a + b mod 2^w (w = len(a) = len(b)); a and the clean ancilla `carry` restored."""
    w = len(a)
    if len(b) != w:
        raise ValueError("add: registers differ in length")
    if w == 1:
        c.cx(a[0], b[0])
        return
    _maj(c, carry, b[0], a[0])
    for i in range(1, w):
        _maj(c, a[i - 1], b[i], a[i])
    for i in range(w - 1, 0, -1):
        _uma(c, a[i - 1], b[i], a[i])
    _uma(c, carry, b[0], a[0])


def sub(c: Circuit, a: list[int], b: list[int], carry: int) -> None:
    """b ← b − a mod 2^w: the adder run backwards."""
    t = Circuit(c.num_qubits)
    add(t, a, b, carry)
    c.extend(t.inverse())


def load_const(c: Circuit, k: int, reg: list[int], ctrl: int | None = None) -> None:
    """reg ^= k, or reg ^= k·ctrl."""
    if k < 0 or k >> len(reg):
        raise ValueError(f"constant {k} does not fit in {len(reg)} bits")
    for i, q in enumerate(reg):
        if k >> i & 1:
            if ctrl is None:
                c.x(q)
            else:
                c.cx(ctrl, q)


def _add_const(c: Circuit, k: int, b: list[int], tmp: list[int], carry: int, ctrl: int | None, subtract: bool = False) -> None:
    load_const(c, k, tmp, ctrl)
    (sub if subtract else add)(c, tmp, b, carry)
    load_const(c, k, tmp, ctrl)


def cc_mod_add(
    c: Circuit,
    k: int,
    N: int,
    ctrls: list[int],
    b: list[int],
    tmp: list[int],
    carry: int,
    flag: int,
    anc: int,
) -> None:
    """b ← b + k·[all ctrls] mod N, for 0 ≤ k < N and b < N held in n + 1 bits (top bit 0).

    tmp (n + 1 bits), carry, flag and anc are clean ancillas and are returned clean. Up to two
    controls; with two, their AND is computed into anc."""
    n1 = len(b)
    if len(tmp) != n1 or N >= 1 << (n1 - 1) or not 0 <= k < N:
        raise ValueError("cc_mod_add: need k < N < 2^(len(b) − 1) and len(tmp) == len(b)")
    if len(ctrls) == 2:
        c.ccx(ctrls[0], ctrls[1], anc)
        g = anc
    elif len(ctrls) == 1:
        g = ctrls[0]
    elif len(ctrls) == 0:
        g = None
    else:
        raise ValueError("cc_mod_add: at most two controls")
    top = b[-1]
    _add_const(c, k, b, tmp, carry, g)                    # b + k·g
    _add_const(c, N, b, tmp, carry, None, subtract=True)  # b + k·g − N
    c.cx(top, flag)                                        # flag = [b + k·g < N]
    _add_const(c, N, b, tmp, carry, flag)                 # (b + k·g) mod N
    _add_const(c, k, b, tmp, carry, g, subtract=True)     # … − k·g: negative iff no wrap
    c.x(top)
    c.cx(top, flag)                                        # clear the flag
    c.x(top)
    _add_const(c, k, b, tmp, carry, g)                    # (b + k·g) mod N
    if len(ctrls) == 2:
        c.ccx(ctrls[0], ctrls[1], anc)


def cc_mod_sub(c: Circuit, k: int, N: int, ctrls, b, tmp, carry, flag, anc) -> None:
    """b ← b − k·[all ctrls] mod N: the inverse of cc_mod_add."""
    t = Circuit(c.num_qubits)
    cc_mod_add(t, k, N, ctrls, b, tmp, carry, flag, anc)
    c.extend(t.inverse())


def c_mod_mul(c: Circuit, a: int, N: int, ctrl: int, x: list[int], b: list[int], tmp, carry, flag, anc) -> None:
    """x ← a·x mod N if ctrl (x < N, gcd(a, N) = 1); b (n + 1 bits) and the ancillas clean."""
    n = len(x)
    if len(b) != n + 1:
        raise ValueError("c_mod_mul: b must have len(x) + 1 bits")
    a %= N
    a_inv = pow(a, -1, N)
    for i in range(n):
        cc_mod_add(c, a * (1 << i) % N, N, [ctrl, x[i]], b, tmp, carry, flag, anc)
    for i in range(n):
        c.cswap(ctrl, x[i], b[i])
    for i in range(n):
        cc_mod_sub(c, a_inv * (1 << i) % N, N, [ctrl, x[i]], b, tmp, carry, flag, anc)


@dataclass
class ShorLayout:
    """Qubits of the order-finding circuit: 3n + 6 for an n-bit modulus."""

    n: int
    ctrl: int
    x: list[int]
    b: list[int]
    tmp: list[int]
    carry: int
    flag: int
    anc: int

    @staticmethod
    def allocate(c: Circuit, n: int) -> "ShorLayout":
        ctrl = c.alloc(1, "ctrl")[0]
        x = c.alloc(n, "x")
        b = c.alloc(n + 1, "b")
        tmp = c.alloc(n + 1, "tmp")
        carry = c.alloc(1, "carry")[0]
        flag = c.alloc(1, "flag")[0]
        anc = c.alloc(1, "anc")[0]
        return ShorLayout(n, ctrl, x, b, tmp, carry, flag, anc)

    def mod_mul(self, c: Circuit, a: int, N: int) -> None:
        c_mod_mul(c, a, N, self.ctrl, self.x, self.b, self.tmp, self.carry, self.flag, self.anc)

    @property
    def ancillas(self) -> list[int]:
        return [*self.b, *self.tmp, self.carry, self.flag, self.anc]
