import numpy as np
import pytest

import ftalgo
from ftalgo.arch import FACTORIES, Architecture, compose_p1
from ftalgo.circuit import Circuit
from ftalgo.estimate import stratified
from ftalgo.schedule import compile_noisy


class FakeModel:
    """A LogicalModel stand-in with simple, known channels."""

    def idle(self, d, p):
        return (1e-3, 1e-4, 2e-3)

    def cnot(self, d, p):
        v = [0.0] * 15
        v[3] = 0.01  # X on the control
        v[14] = 0.005  # ZZ
        return v

    def pauli_measurement(self, kind, d, p):
        j = np.zeros(8)
        j[1] = 0.002  # wrong outcome
        j[2] = 0.003  # X on patch 1
        j[0] = 1 - j.sum()
        return {"joint": j, "o": 0.004}

    def extrapolated(self, d, p):
        return False


def arch(d=5):
    return Architecture(FakeModel(), d, 0.001, FACTORIES["injected"])


def test_compose_p1():
    q = (0.01, 0.002, 0.03)
    two = compose_p1(q, 2)
    # Direct convolution of the Pauli group: I X Y Z multiply as the Klein group.
    p = {"I": 1 - sum(q), "X": q[0], "Y": q[1], "Z": q[2]}
    mul = {("X", "Y"): "Z", ("Y", "X"): "Z", ("X", "Z"): "Y", ("Z", "X"): "Y", ("Y", "Z"): "X", ("Z", "Y"): "X"}
    out = {k: 0.0 for k in "IXYZ"}
    for a, pa in p.items():
        for b, pb in p.items():
            c = "I" if a == b else (b if a == "I" else a if b == "I" else mul[(a, b)])
            out[c] += pa * pb
    assert two == pytest.approx((out["X"], out["Y"], out["Z"]))


def test_idle_gaps_and_durations():
    c = Circuit(2)
    c.cx(0, 1)  # both start at 0, end at 2d = 10
    c.t(0)  # 10 → 15
    c.t(0)  # 15 → 20
    c.cx(0, 1)  # qubit 1 waited 10 → idle(10)
    np_ = compile_noisy(c, arch())
    assert np_.rounds == 30
    texts = np_.text
    idle10 = compose_p1((1e-3, 1e-4, 2e-3), 10)
    assert f"PAULI_CHANNEL_1({idle10[0]!r},{idle10[1]!r},{idle10[2]!r}) 1" in texts
    prog = ftalgo.Program(texts)
    assert prog.expected_faults == pytest.approx(sum(np_.budget.values()), rel=1e-9)


def test_table_padding_and_reaction():
    c = Circuit(2)
    c.h(0)
    r = c.measure(0)
    a, b = Circuit(2), Circuit(2)
    a.t(1)
    b.t(1)
    b.t(1)
    b.s(1)
    c.table([r], [a, b])
    c.measure(1)
    np_ = compile_noisy(c, arch())
    # Measure at 0 → 1; reaction 10 → table starts at 11; longest case 3·5 = 15 → 26; M → 27.
    assert np_.rounds == 27
    prog = ftalgo.Program(np_.text)
    assert prog.num_measurements == 2


def test_repeat_barrier():
    body = Circuit(3)
    body.cx(0, 1)
    c = Circuit(3)
    c.t(2)  # qubit 2 live from 0 to 5
    c.repeat(4, body)
    np_ = compile_noisy(c, arch())
    # Barrier at 5 (qubits 0, 1 start there: not yet live), 4 × 10 rounds.
    assert np_.rounds == 5 + 40
    prog = ftalgo.Program(np_.text)
    assert prog.expected_faults == pytest.approx(sum(np_.budget.values()), rel=1e-9)
    assert np_.counts["CX"] == 4


def test_stratified_matches_plain():
    c = Circuit(3)
    for _ in range(6):
        c.h(0)
        c.cx(0, 1)
        c.ccx(0, 1, 2)
        c.t(2)
    for q in range(3):
        c.measure(q)
    np_ = compile_noisy(c, arch())
    prog = ftalgo.Program(np_.text)
    ideal = ftalgo.Program(c.to_text())
    ref_rows = ideal.sample(1, seed=0, faults="none")[0]

    def score(rows):
        return (rows == ref_rows).all(axis=1)

    # The circuit is deterministic only if H's randomness cancels; check against plain sampling.
    plain = score(prog.sample(400_000, seed=1)).mean()
    est = stratified(prog, score, seed=2, shots=200_000)
    assert abs(est.value - plain) < 4 * (est.sigma + np.sqrt(plain * (1 - plain) / 400_000)) + est.tail
