"""Exhaustive checks of the arithmetic blocks with the bit-level reversible simulator."""

from math import gcd

import pytest

from ftalgo import arith
from ftalgo.circuit import Circuit
from ftalgo.revsim import get, put, run_reversible


@pytest.mark.parametrize("w", [1, 2, 3, 4, 5])
def test_add_and_sub_every_input(w):
    c = Circuit()
    a = c.alloc(w, "a")
    b = c.alloc(w, "b")
    carry = c.alloc(1, "carry")[0]
    add_c, sub_c = Circuit(c.num_qubits), Circuit(c.num_qubits)
    arith.add(add_c, a, b, carry)
    arith.sub(sub_c, a, b, carry)
    if w > 1:
        assert add_c.toffoli_count() == 2 * w
    for x in range(1 << w):
        for y in range(1 << w):
            s = put(put(0, a, x), b, y)
            out = run_reversible(add_c, s)
            assert get(out, a) == x and get(out, b) == (x + y) % (1 << w) and not out >> carry & 1
            out = run_reversible(sub_c, s)
            assert get(out, a) == x and get(out, b) == (y - x) % (1 << w) and not out >> carry & 1


def test_load_const():
    c = Circuit()
    r = c.alloc(4, "r")
    ctrl = c.alloc(1, "c")[0]
    arith.load_const(c, 0b1011, r, ctrl)
    assert get(run_reversible(c, 0), r) == 0
    assert get(run_reversible(c, 1 << ctrl), r) == 0b1011
    with pytest.raises(ValueError):
        arith.load_const(c, 16, r)


def mod_add_setup(N):
    n = N.bit_length()
    c = Circuit()
    ctl = c.alloc(2, "ctl")
    b = c.alloc(n + 1, "b")
    tmp = c.alloc(n + 1, "tmp")
    carry, flag, anc = c.alloc(3, "anc")
    return c, ctl, b, tmp, carry, flag, anc


@pytest.mark.parametrize("N", [3, 5, 7, 11, 13, 15, 21])
def test_cc_mod_add_every_input(N):
    c0, ctl, b, tmp, carry, flag, anc = mod_add_setup(N)
    clean = [*tmp, carry, flag, anc]
    for k in range(N):
        add_c = Circuit(c0.num_qubits)
        arith.cc_mod_add(add_c, k, N, ctl, b, tmp, carry, flag, anc)
        sub_c = Circuit(c0.num_qubits)
        arith.cc_mod_sub(sub_c, k, N, ctl, b, tmp, carry, flag, anc)
        for y in range(N):
            for cv in range(4):
                s = put(put(0, ctl, cv), b, y)
                on = cv == 3
                out = run_reversible(add_c, s)
                assert get(out, b) == ((y + k) % N if on else y), (N, k, y, cv)
                assert get(out, ctl) == cv and all(not out >> q & 1 for q in clean)
                out = run_reversible(sub_c, s)
                assert get(out, b) == ((y - k) % N if on else y)
                assert all(not out >> q & 1 for q in clean)


@pytest.mark.parametrize("N", [15, 21, 35])
def test_c_mod_mul_every_input(N):
    n = N.bit_length()
    c0 = Circuit()
    lay = arith.ShorLayout.allocate(c0, n)
    for a in range(2, N):
        if gcd(a, N) != 1:
            continue
        mul = Circuit(c0.num_qubits)
        lay.mod_mul(mul, a, N)
        assert mul.toffoli_count() == 2 * n * (10 * n + 12) + n
        for x in range(N):
            for ctrl in (0, 1):
                s = put(0, lay.x, x) | (ctrl << lay.ctrl)
                out = run_reversible(mul, s)
                assert get(out, lay.x) == (a * x % N if ctrl else x), (N, a, x, ctrl)
                assert out >> lay.ctrl & 1 == ctrl
                assert all(not out >> q & 1 for q in lay.ancillas), (N, a, x, ctrl)
