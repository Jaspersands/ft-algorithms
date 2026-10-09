"""Exhaustive bit-level and quantum checks of modern arithmetic modules."""

from math import gcd

import pytest

from ftalgo import arith_modern as am
from ftalgo.circuit import Circuit
from ftalgo.engine import Program
from ftalgo.revsim import get, put, run_reversible


@pytest.mark.parametrize("w", [1, 2, 3, 4, 5])
def test_gidney_add_and_sub_every_input(w):
    c = Circuit()
    a = c.alloc(w, "a")
    b = c.alloc(w, "b")
    carries = c.alloc(max(w - 1, 0), "carries")
    add_c, sub_c = Circuit(c.num_qubits), Circuit(c.num_qubits)
    am.gidney_add(add_c, a, b, carries)
    am.gidney_sub(sub_c, a, b, carries)
    if w > 1:
        assert add_c.toffoli_count() == w - 1
        assert sub_c.toffoli_count() == w - 1

    for x in range(1 << w):
        for y in range(1 << w):
            s = put(put(0, a, x), b, y)
            out = run_reversible(add_c, s, ignore_phases=True)
            assert get(out, a) == x
            assert get(out, b) == (x + y) % (1 << w)
            assert all(not (out >> q & 1) for q in carries)

            out_sub = run_reversible(sub_c, s, ignore_phases=True)
            assert get(out_sub, a) == x
            assert get(out_sub, b) == (y - x) % (1 << w)
            assert all(not (out_sub >> q & 1) for q in carries)


def test_window_table_load_unload():
    for C in [[0, 5, 10, 15], [0, 3, 6, 9], [0, 7, 2, 9]]:
        for ctrl in (0, 1):
            for xv in range(4):
                c = Circuit()
                ctl = c.alloc(1, "c")[0]
                x = c.alloc(2, "x")
                tmp = c.alloc(5, "tmp")
                anc = c.alloc(3, "anc")
                am.load_window_table(c, C, tmp, ctl, x, anc)
                s = put(put(0, [ctl], ctrl), x, xv)
                out = run_reversible(c, s, ignore_phases=True)
                expected = C[xv] if ctrl else 0
                assert get(out, tmp) == expected

                c_un = Circuit(c.num_qubits)
                am.load_window_table(c_un, C, tmp, ctl, x, anc)
                am.unload_window_table(c_un, C, tmp, ctl, x, anc)
                out2 = run_reversible(c_un, s, ignore_phases=True)
                assert get(out2, tmp) == 0
                assert all(not (out2 >> q & 1) for q in anc)


@pytest.mark.parametrize("N", [3, 5, 7, 15])
def test_window_mod_add_every_input(N):
    n = N.bit_length()
    n1 = n + 1
    C = [0, 1 % N, 2 % N, 3 % N]
    for ctrl in (0, 1):
        for xv in range(4):
            c = Circuit()
            ctl = c.alloc(1, "ctl")[0]
            x = c.alloc(2, "x")
            b = c.alloc(n1, "b")
            tmp = c.alloc(n1, "tmp")
            carries = c.alloc(n1 - 1, "carries")
            flag = c.alloc(1, "flag")[0]
            anc = c.alloc(3, "anc")
            am.window_mod_add(c, C, N, ctl, x, b, tmp, carries, flag, anc, subtract=False)
            clean = [*tmp, *carries, flag, *anc]

            c_sub = Circuit(c.num_qubits)
            am.window_mod_add(c_sub, C, N, ctl, x, b, tmp, carries, flag, anc, subtract=True)

            for y in range(N):
                s = put(put(put(0, [ctl], ctrl), x, xv), b, y)
                term = C[xv] if ctrl else 0
                out = run_reversible(c, s, ignore_phases=True)
                assert get(out, b) == (y + term) % N
                assert all(not (out >> q & 1) for q in clean)

                out_sub = run_reversible(c_sub, s, ignore_phases=True)
                assert get(out_sub, b) == (y - term) % N
                assert all(not (out_sub >> q & 1) for q in clean)


@pytest.mark.parametrize("N", [15, 21])
def test_windowed_c_mod_mul_every_input(N):
    n = N.bit_length()
    c0 = Circuit()
    lay = am.ModernShorLayout.allocate(c0, n)
    for a in range(2, N):
        if gcd(a, N) != 1:
            continue
        mul = Circuit(c0.num_qubits)
        lay.mod_mul(mul, a, N, k=2)

        # Verify against analytical formula for Toffolis
        num_w = (n + 1) // 2
        exp_toffolis = 2 * (num_w * (5 * n + 18) if n % 2 == 0 else (num_w - 1) * (5 * n + 18) + (5 * n + 6)) + n
        assert mul.toffoli_count() == exp_toffolis

        for x in range(N):
            for ctrl in (0, 1):
                s = put(0, lay.x, x) | (ctrl << lay.ctrl)
                out = run_reversible(mul, s, ignore_phases=True)
                assert get(out, lay.x) == (a * x % N if ctrl else x), (N, a, x, ctrl)
                assert (out >> lay.ctrl & 1) == ctrl
                assert all(not (out >> q & 1) for q in lay.ancillas), (N, a, x, ctrl)


def test_modern_shor_quantum_simulation():
    import numpy as np

    # Test N = 15 with a = 7 on ftsim
    c, info = am.order_finding_modern(15, 7, m=4, cutoff=2)
    p = Program(c.to_text())
    assert p.num_qubits == c.num_qubits
    amps, recs = p.final_state()
    # 4 QPE rounds
    assert len(info.records) == 4
    outcomes = info.outcomes(np.array([recs]))
    # For a = 7 mod 15 (r = 4), outcome must be a multiple of 2^(m)/r = 16/4 = 4
    assert outcomes[0] % 4 == 0
