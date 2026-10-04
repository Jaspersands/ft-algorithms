import numpy as np
import pytest

import ftalgo
from ftalgo.circuit import Circuit
from ftalgo.revsim import get, put, run_reversible

ONE = ["X", "Y", "Z", "H", "S", "S_DAG", "T", "T_DAG", "SQRT_X", "SQRT_X_DAG"]
TWO = ["CX", "CY", "CZ", "SWAP"]


def random_circuit(rng, n, count, gates1=ONE, gates2=TWO, gates3=("CCX", "CCZ")):
    c = Circuit(n)
    for _ in range(count):
        r = rng.random()
        if r < 0.5:
            c.gate(rng.choice(gates1), int(rng.integers(n)))
        elif r < 0.8:
            a, b = rng.choice(n, 2, replace=False)
            c.gate(rng.choice(gates2), int(a), int(b))
        else:
            a, b, t = rng.choice(n, 3, replace=False)
            c.gate(rng.choice(gates3), int(a), int(b), int(t))
    return c


def test_text_round_trip_and_table_offsets():
    c = Circuit()
    q = c.alloc(4, "q")
    c.x(q[0])
    c.x(q[2])
    r0 = c.measure(q[0])
    c.measure(q[1])
    r2 = c.measure(q[2])
    cases = []
    for i in range(4):
        k = Circuit(4)
        if i == 3:  # both records 1 → flip qubit 3
            k.x(q[3])
        cases.append(k)
    c.table([r0, r2], cases)
    c.measure(q[3])
    text = c.to_text()
    assert "TABLE rec[-3] rec[-1]" in text
    p = ftalgo.Program(text)
    rec = p.sample(5, seed=0)
    assert rec.tolist() == [[True, False, True, True]] * 5


def test_counts_through_repeat_and_table():
    body = Circuit(2)
    body.t(0)
    body.ccx  # attribute exists
    body.cx(0, 1)
    c = Circuit(3)
    c.repeat(10, body)
    c.measure(0)
    a, b = Circuit(3), Circuit(3)
    a.t(2)
    b.t(2)
    b.t(1)
    b.ccx(0, 1, 2)
    c.table([0], [a, b])
    cnt = c.counts()
    assert cnt["T"] == 12 and cnt["CX"] == 10 and cnt["CCX"] == 1 and cnt["M"] == 1
    assert c.t_count() == 12 and c.toffoli_count() == 1


@pytest.mark.parametrize("seed", range(20))
def test_inverse_is_identity(seed):
    rng = np.random.default_rng(seed)
    c = random_circuit(rng, 4, 60)
    body = random_circuit(rng, 4, 5)
    c.repeat(3, body)
    full = Circuit(4)
    full.h(0)
    full.h(2)
    prep = full.to_text()
    full.extend(c)
    full.extend(c.inverse())
    v = ftalgo.Program(full.to_text()).statevector()
    w = ftalgo.Program(prep + "I 3\n").statevector()
    assert np.allclose(v, w, atol=1e-10)


@pytest.mark.parametrize("seed", range(10))
def test_revsim_matches_engine(seed):
    rng = np.random.default_rng(100 + seed)
    c = random_circuit(rng, 5, 40, gates1=["X"], gates2=["CX", "SWAP"], gates3=["CCX"])
    for x in range(32):
        prep = Circuit(5)
        for q in range(5):
            if x >> q & 1:
                prep.x(q)
        prep.extend(c)
        for q in range(5):
            prep.measure(q)
        rec = ftalgo.Program(prep.to_text()).sample(1, seed=0)[0]
        got = sum(int(b) << i for i, b in enumerate(rec))
        assert got == run_reversible(c, x)


def test_register_helpers_and_errors():
    s = put(0, [3, 1, 4], 0b101)
    assert s == (1 << 3) | (1 << 4)
    assert get(s, [3, 1, 4]) == 0b101
    c = Circuit()
    with pytest.raises(ValueError):
        c.gate("CX", 1, 1)
    with pytest.raises(ValueError):
        c.gate("FOO", 1)
    c.h(0)
    with pytest.raises(ValueError):
        run_reversible(c, 0)
    with pytest.raises(ValueError):
        c.table([0], [Circuit(), Circuit()])
