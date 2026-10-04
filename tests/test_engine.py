import numpy as np
import pytest

import ftalgo


def test_bell_pairs():
    p = ftalgo.Program("H 0\nCX 0 1\nM 0 1")
    r = p.sample(1000, seed=1)
    assert r.shape == (1000, 2)
    assert (r[:, 0] == r[:, 1]).all()
    assert 400 < r[:, 0].sum() < 600


def test_determinism_across_threads_and_backends():
    text = "H 0 1 2\nREPEAT 30 {\nCCX 0 1 2\nT 0\nH 1\nDEPOLARIZE1(0.05) 0 1 2\nM 2\nTABLE rec[-1] {S 0} {H 0}\n}\nM 0 1"
    p = ftalgo.Program(text)
    a = p.sample(500, seed=3, threads=1)
    for kw in ({"threads": 4}, {"threads": 0}, {"backend": "dense"}, {"backend": "sparse"}):
        assert np.array_equal(a, p.sample(500, seed=3, **kw)), kw
    b = np.concatenate([p.sample(200, seed=3), p.sample(300, seed=3, first_shot=200)])
    assert np.array_equal(a, b)


def test_fault_logs_and_counts():
    p = ftalgo.Program("MARK(cnot)\nDEPOLARIZE2(0.01) 0 1\nMARK(idle)\nX_ERROR(0.02) 0 1 2\nM 0 1 2")
    assert p.num_sites == 4
    assert p.classes == ["unmarked", "cnot", "idle"]
    dist = p.fault_count_distribution(4)
    qs = [0.01, 0.02, 0.02, 0.02]
    assert dist[0] == pytest.approx(np.prod([1 - q for q in qs]))
    assert dist.sum() == pytest.approx(1.0)
    rec, f = p.sample(2000, seed=5, faults=2, log_faults=True)
    assert len(f.shot) == 4000
    assert np.all(np.bincount(f.shot.astype(np.int64), minlength=2000) == 2)
    assert set(np.unique(f.cls)) <= {1, 2}


def test_errors_raise_value_error():
    with pytest.raises(ValueError, match="line 2"):
        ftalgo.Program("H 0\nCX 0 0")
    p = ftalgo.Program("H 0\nM 0")
    with pytest.raises(ValueError):
        p.sample(10, seed=0, faults=1)
    with pytest.raises(ValueError):
        p.sample(10, seed=0, faults="sometimes")
    with pytest.raises(ValueError):
        p.sample(10, seed=0, backend="gpu")


def test_statevector():
    v = ftalgo.Program("H 0\nT 0").statevector()
    assert np.allclose(v, np.array([1, np.exp(1j * np.pi / 4)]) / np.sqrt(2))
    amps, _ = ftalgo.Program("X 100").final_state()
    assert amps == {1 << 100: 1}
