import pytest
import math

from ftalgo.biased import (
    BiasedNoise,
    XZZXPatch,
    logical_rates_analytic,
    logical_rates_calibrated,
    optimize_patch,
    evaluate_bias_scaling,
)
from ftalgo.calib.model import LogicalModel


def test_biased_noise_probabilities():
    noise = BiasedNoise(p_total=1e-3, eta=100.0)
    assert math.isclose(noise.px + noise.py + noise.pz, 1e-3, rel_tol=1e-9)
    assert noise.px == noise.py
    assert math.isclose(noise.pz / (2 * noise.px), 100.0, rel_tol=1e-9)
    assert noise.px < 5e-6
    assert noise.pz > 9.9e-4


def test_xzzx_patch_dimensions():
    patch = XZZXPatch(dx=25, dz=7)
    assert patch.tile_qubits == 2 * (25 + 1) * (7 + 1)  # 2 * 26 * 8 = 416
    assert patch.data_qubits == (25 + 1) * (7 + 1)  # 208

    with pytest.raises(ValueError):
        XZZXPatch(dx=24, dz=7)
    with pytest.raises(ValueError):
        XZZXPatch(dx=25, dz=2)


def test_logical_rates():
    noise = BiasedNoise(p_total=1e-3, eta=100.0)
    patch = XZZXPatch(dx=25, dz=7)
    rates = logical_rates_analytic(patch, noise)
    assert rates.pz_logical > 0
    assert rates.px_logical > 0
    assert rates.ptot_logical == rates.pz_logical + rates.px_logical
    # With dz=7 and px=5e-6 vs p*=0.0109, px_log should be tiny
    assert rates.px_logical < 1e-12

    model = LogicalModel.load()
    calib_rates = logical_rates_calibrated(patch, noise, model)
    assert calib_rates.ptot_logical > 0
    assert math.isclose(calib_rates.pz_logical, rates.pz_logical, rel_tol=0.2)


def test_optimize_patch():
    target = 1e-10
    patch_sym, _ = optimize_patch(target, p_total=1e-3, eta=0.5)
    patch_biased, _ = optimize_patch(target, p_total=1e-3, eta=100.0)

    # dx should be equal or close (since pz ≈ p), but dz should be much smaller
    assert patch_biased.dx >= patch_biased.dz
    assert patch_biased.dz < patch_sym.dz
    assert patch_biased.tile_qubits < patch_sym.tile_qubits


def test_evaluate_bias_scaling():
    model = LogicalModel.load()
    rows = evaluate_bias_scaling(
        total_space_time_rounds=1e6,
        num_logical_qubits=100,
        target_faults=0.01,
        p_total=1e-3,
        etas=[1.0, 10.0, 50.0, 100.0, 500.0],
        model=model,
    )
    assert len(rows) == 5
    # Reductions should grow with eta
    for r in rows:
        assert r.target_error_achieved <= 0.01
    assert rows[-1].qubit_reduction_ratio > rows[0].qubit_reduction_ratio
    # For eta=500, reduction ratio should be > 2.5x
    assert rows[-1].qubit_reduction_ratio >= 2.5
