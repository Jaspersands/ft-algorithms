import math

import pytest

from ftalgo import scale
from ftalgo.arch import ROUND_SECONDS
from ftalgo.calib.model import LogicalModel
from ftalgo.sensitivity import GidneyProjection, IdleBand, SteeperTail, fit_covariance


@pytest.fixture(scope="module")
def model():
    return LogicalModel.load()


def test_gidney_workload_is_one_shot():
    w = scale.GIDNEY_2025
    assert w.rounds_fixed * ROUND_SECONDS == pytest.approx(12.07 * 3600)
    assert w.toffolis * scale.GIDNEY_2025_SHOTS == pytest.approx(6.5e9)
    # 1e-15 per logical-qubit round over 6.9e13 qubit-rounds: 93.3% of shots succeed.
    assert math.exp(-scale.GIDNEY_2025_PER_ROUND * 6.9e13) == pytest.approx(0.933, abs=1e-3)


def test_gidney_projection_reproduces_his_distance():
    g = GidneyProjection()
    # His Section 3.2: d = 25 is where his projection first reaches 1e-15 per patch-round.
    assert scale.distance_for_idle(g, scale.GIDNEY_2025_PER_ROUND) == scale.GIDNEY_2025_D
    # The estimator, given his projection and his budget, returns his distance.
    e = scale.estimate(scale.GIDNEY_2025, g, target=scale.GIDNEY_2025_TARGET)
    assert e.d == scale.GIDNEY_2025_D


def test_fit_covariance_reproduces_model_fit(model):
    for name in ("idle_xy", "idle_zy"):
        fc = fit_covariance(model, name, 1e-3)
        f = model.per_p[name][1e-3]
        assert fc["a"] == pytest.approx(f.a) and fc["b"] == pytest.approx(f.b)
        assert fc["lambda_lo"] < fc["lambda"] < fc["lambda_hi"]


def test_idle_band_brackets_the_fit_and_widens_with_extrapolation(model):
    lo, hi = IdleBand(model, -1), IdleBand(model, 1)
    width = []
    for d in (7, 25):
        c = scale.idle_total(model, d, 1e-3)
        assert scale.idle_total(lo, d, 1e-3) < c < scale.idle_total(hi, d, 1e-3)
        width.append(scale.idle_total(hi, d, 1e-3) / scale.idle_total(lo, d, 1e-3))
    assert width[1] > width[0]


def test_steeper_tail_keeps_measured_points(model):
    s = SteeperTail(model, 1.15)
    dmax = model.per_p["idle_xy"][1e-3].dmax
    assert s.rate("idle_xy", dmax, 1e-3)[0] == pytest.approx(model.rate("idle_xy", dmax, 1e-3)[0])
    assert s.rate("idle_xy", dmax + 10, 1e-3)[0] < model.rate("idle_xy", dmax + 10, 1e-3)[0]
    assert s.rate("cnot_x2", 25, 1e-3) == model.rate("cnot_x2", 25, 1e-3)


@pytest.mark.parametrize("N", [15, 21, 35, 77, 143])
def test_textbook_counts_used_at_scale_match_compiled_circuits(N):
    from ftalgo import shor

    c, info = shor.order_finding(N)
    assert c.toffoli_count() == scale.textbook_toffolis(info.n)
    assert c.num_qubits == 3 * info.n + 6
