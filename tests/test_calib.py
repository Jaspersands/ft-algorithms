import numpy as np
import pytest

sq = pytest.importorskip("stabilizer_qec")

from ftalgo.calib import channels  # noqa: E402
from ftalgo.calib.choi import choi_memory  # noqa: E402
from ftalgo.calib.qec import Counts, run  # noqa: E402


def test_run_matches_direct_decoding():
    c = sq.memory_circuit(distance=3, rounds=3, p=0.005)
    res = run(c, shots=1 << 17, seed=11)
    dets, obs = c.compile_detector_sampler(seed=11).sample(1 << 17, separate_observables=True)
    pred = sq.Matching(c.detector_error_model(decompose_errors=True), enable_correlations=True).decode_batch(dets)
    direct = (pred != obs).any(axis=1).sum()
    assert res.failures == direct  # same seed, same stream, same decoder
    j = Counts.from_json(res.to_json())
    assert j.patterns == res.patterns and j.shots == res.shots


def test_early_stop():
    c = sq.memory_circuit(distance=3, rounds=3, p=0.01)
    res = run(c, shots=1 << 22, seed=1, target=100, batch=4096)
    assert res.failures >= 100 and res.shots < 1 << 22 and res.shots % 4096 == 0


def test_fidelity_round_trip_and_composition():
    rng = np.random.default_rng(0)
    for n in (1, 2, 3):
        a = rng.dirichlet(np.ones(1 << n) * 0.3)
        b = rng.dirichlet(np.ones(1 << n) * 0.3)
        assert np.allclose(channels.probabilities(channels.fidelities(a), clip=False), a)
        conv = np.zeros(1 << n)
        for x in range(1 << n):
            for y in range(1 << n):
                conv[x ^ y] += a[x] * b[y]
        assert np.allclose(channels.fidelities(conv), channels.fidelities(a) * channels.fidelities(b))
        assert np.allclose(channels.probabilities(channels.deconvolve(channels.fidelities(conv), channels.fidelities(b)), clip=False), a)


def test_per_round_from_two_lengths():
    rnd = np.array([0.97, 0.01, 0.015, 0.005])
    prep = np.array([0.9, 0.05, 0.03, 0.02])
    f = lambda t: channels.fidelities(prep) * channels.fidelities(rnd) ** t  # noqa: E731
    assert np.allclose(channels.probabilities(channels.per_round(f(3), f(9), 3, 9)), rnd)


def test_pauli_channel_2_marginals():
    x = np.array([0.9, 0.04, 0.05, 0.01])
    z = np.array([0.8, 0.1, 0.07, 0.03])
    v = channels.pauli_channel_2(x, z)
    assert len(v) == 15
    assert sum(v) == pytest.approx(1 - x[0] * z[0])
    # X on the first qubit alone: a = X (1), b = I (0) → index 4·1 + 0 − 1 = 3.
    assert v[3] == pytest.approx(x[1] * z[0])


def test_choi_memory_is_valid_and_consistent():
    d, T, p = 3, 6, 0.004
    text = choi_memory(d, T, p)
    c = sq.Circuit(text)
    assert c.num_observables == 2
    c.detector_error_model(decompose_errors=True)  # raises if anything is non-deterministic
    noiseless = sq.Circuit(choi_memory(d, T, 0.0))
    dets, obs = noiseless.compile_detector_sampler(seed=0).sample(256, separate_observables=True)
    assert not dets.any() and not obs.any()
    shots = 1 << 19
    res = run(c, shots=shots, seed=5)
    zmem = run(sq.memory_circuit(distance=d, rounds=T, p=p, basis="z"), shots=shots, seed=6)
    xmem = run(sq.memory_circuit(distance=d, rounds=T, p=p, basis="x"), shots=shots, seed=7)
    # obs0 (X_R X_L) sees Z/Y like the X memory; obs1 sees X/Y like the Z memory. The Choi
    # experiment ends in a perfect round instead of a noisy readout, so it fails a little less.
    for got, mem in ((res.marginal(0), xmem.marginal(0)), (res.marginal(1), zmem.marginal(0))):
        assert 0.6 * mem < got <= mem * 1.05 + 5 * np.sqrt(mem / shots)
