import numpy as np
import pytest

from ftalgo import shor
from ftalgo.calib.model import LogicalModel
from ftalgo.estimate import _pool
from ftalgo.experiment import merge_runs, run_config


def test_pool_matches_concatenation():
    rng = np.random.default_rng(1)
    groups = [rng.random((n, 2)) for n in (5, 17, 40)]
    mean, var = _pool(np.array([len(g) for g in groups], float), np.array([g.mean(axis=0) for g in groups]),
                      np.array([g.std(axis=0, ddof=1) for g in groups]))
    allv = np.concatenate(groups)
    assert np.allclose(mean, allv.mean(axis=0)) and np.allclose(var, allv.var(axis=0, ddof=1))


@pytest.fixture(scope="module")
def n15():
    circ, info = shor.order_finding(15)
    peaks = shor.peak_table(15, info.a, info.m)
    return circ, (lambda rows: np.stack([peaks[info.outcomes(rows)], peaks[info.outcomes(rows)]], axis=1)), LogicalModel.load()


@pytest.mark.parametrize("d", [7, 11])  # plain sampling at d = 7, stratified at d = 11
def test_chunks_merge_like_one_run(n15, d):
    circ, score, model = n15
    one = run_config(circ, score, model, d=d, p=0.001, factory="injected", seed=7, shots=800, min_shots=40)
    assert merge_runs([one])["scores"] == pytest.approx(one["scores"])
    assert merge_runs([one])["score_sigmas"] == pytest.approx(one["score_sigmas"])
    parts = [run_config(circ, score, model, d=d, p=0.001, factory="injected", seed=100 + c, shots=200, min_shots=10) for c in range(4)]
    m = merge_runs(parts)
    assert m["estimate"]["method"] == one["estimate"]["method"]
    assert m["chunks"] == 4 and m["estimate"]["shots"] == sum(q["estimate"]["shots"] for q in parts)
    s = np.hypot(m["sigma"], one["sigma"])
    assert abs(m["success"] - one["success"]) < 5 * s + 1e-9
    assert m["sigma"] < max(q["sigma"] for q in parts)
