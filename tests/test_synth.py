import math

import numpy as np
import pytest

import ftalgo
from ftalgo import synth
from ftalgo.circuit import Circuit


@pytest.fixture(autouse=True)
def tmp_cache(tmp_path, monkeypatch):
    monkeypatch.setattr(synth, "_CACHE_PATH", str(tmp_path / "cache.json"))
    monkeypatch.setattr(synth, "_cache", None)
    yield


def test_exact_multiples_of_pi_over_4():
    for k in range(-16, 17):
        seq = synth.rz(k * math.pi / 4, 1e-9)
        assert len(seq) <= 2
        assert synth.distance_up_to_phase(synth.matrix(seq), synth.target(k * math.pi / 4)) < 1e-12


@pytest.mark.parametrize("eps", [1e-3, 1e-6])
def test_random_angles_within_eps(eps):
    rng = np.random.default_rng(int(-math.log10(eps)))
    for theta in rng.uniform(-math.pi, math.pi, size=25):
        seq = synth.rz(theta, eps)
        assert set(seq) <= {"H", "S", "T", "X", "Y", "Z", "S_DAG", "T_DAG"}
        d = synth.distance_up_to_phase(synth.matrix(seq), synth.target(theta))
        assert d <= eps
        # T-count grows like 3·log2(1/ε) (Ross–Selinger): a loose sanity bound.
        assert synth.t_count(seq) <= 4 * math.log2(1 / eps) + 20


def test_cache_round_trip():
    a = synth.rz(0.123, 1e-4)
    synth.save_cache()
    synth._cache = None
    b = synth.rz(0.123, 1e-4)
    assert a == b


def test_engine_agrees():
    theta = 1.234
    seq = synth.rz(theta, 1e-8)
    c = Circuit(1)
    c.h(0)
    c.gates(seq, 0)
    v = ftalgo.Program(c.to_text()).statevector()
    want = np.array([1, np.exp(1j * theta)]) / math.sqrt(2)
    overlap = abs(np.vdot(want, v))
    assert overlap > 1 - 1e-12
