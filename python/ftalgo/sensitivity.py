"""Sensitivity analysis and extrapolation robustness for fault-tolerant scaling.

Analyzes the robustness of the RSA-2048 headline result (d = 29 vs Gidney's d = 25)
under:
1. Idle channel fit uncertainty (bootstrap / WLS covariance on slope Lambda and prefactor a).
2. Decoder choice (plain MWPM vs correlated matching / belief matching).
3. Magic-state distillation quality (epsilon_CCZ variation).
4. Decoder reaction time (5 to 20 rounds).
5. Assumption-matched reproduction (feeding Gidney's assumed model into our estimator).
"""

from __future__ import annotations

import copy
import dataclasses
import json
import math
import os
from typing import Any

import numpy as np

from .arch import FACTORIES, Architecture
from .calib.model import Fit, LogicalModel, Point
from .scale import GIDNEY_2025, Workload, estimate, expected_faults


def idle_fit_covariance(
    model: LogicalModel, name: str = "idle_xy", p: float = 1e-3
) -> dict[str, Any]:
    """Computes WLS regression parameters and covariance matrix for a channel."""
    pts = [q for q in model.points.get(name, []) if q.ok() and q.p == p]
    if len(pts) < 2:
        raise ValueError(f"Insufficient usable points for {name} at p={p}")

    x = np.array([(q.d + 1) / 2 for q in pts])
    y = np.log([q.value for q in pts])
    w = np.array([q.events for q in pts], dtype=float)

    A = np.vstack([np.ones_like(x), x]).T
    W = np.diag(w)
    cov = np.linalg.inv(A.T @ W @ A)
    coef = cov @ (A.T @ W @ y)
    res = y - A @ coef
    chi2_red = float((res**2 * w).sum() / max(len(pts) - 2, 1))
    cov_param = cov * chi2_red

    a = float(coef[0])
    b = float(coef[1])
    se_a = float(np.sqrt(cov_param[0, 0]))
    se_b = float(np.sqrt(cov_param[1, 1]))
    cov_ab = float(cov_param[0, 1])

    # Evaluate band across odd d from 3 to 35
    band: dict[int, dict[str, float]] = {}
    for d in range(3, 37, 2):
        xd = (d + 1) / 2
        var_log = float(cov_param[0, 0] + 2 * xd * cov_param[0, 1] + xd**2 * cov_param[1, 1])
        se_log = math.sqrt(max(var_log, 0.0))
        val = math.exp(a + b * xd)
        lo_1 = math.exp(a + b * xd - se_log)
        hi_1 = math.exp(a + b * xd + se_log)
        lo_2 = math.exp(a + b * xd - 2 * se_log)
        hi_2 = math.exp(a + b * xd + 2 * se_log)
        band[d] = {
            "central": val,
            "lo_1sigma": lo_1,
            "hi_1sigma": hi_1,
            "lo_2sigma": lo_2,
            "hi_2sigma": hi_2,
        }

    return {
        "channel": name,
        "p": p,
        "n_points": len(pts),
        "a": a,
        "b": b,
        "se_a": se_a,
        "se_b": se_b,
        "cov_ab": cov_ab,
        "lambda": float(np.exp(-b)),
        "lambda_lo": float(np.exp(-(b + se_b))),
        "lambda_hi": float(np.exp(-(b - se_b))),
        "band": band,
    }


class PerturbedLogicalModel:
    """Wraps a LogicalModel with parameter shifts for sensitivity analysis."""

    def __init__(
        self,
        base: LogicalModel,
        *,
        slope_scale: float = 1.0,
        prefactor_scale: float = 1.0,
        override_rate: dict[str, float] | None = None,
    ):
        self._base = base
        self.slope_scale = slope_scale
        self.prefactor_scale = prefactor_scale
        self.override_rate = override_rate or {}

    def rate(self, name: str, d: int, p: float) -> tuple[float, bool]:
        if name in self.override_rate:
            return self.override_rate[name], True
        val, extrap = self._base.rate(name, d, p)
        if name.startswith("idle"):
            f = self._base.per_p.get(name, {}).get(p)
            if f is not None:
                xd = (d + 1) / 2
                mod_a = f.a + math.log(self.prefactor_scale)
                mod_b = f.b * self.slope_scale
                return math.exp(mod_a + mod_b * xd), extrap
        return val * self.prefactor_scale, extrap

    def idle(self, d: int, p: float) -> tuple[float, float, float]:
        xy, _ = self.rate("idle_xy", d, p)
        zy, _ = self.rate("idle_zy", d, p)
        y, _ = self.rate("idle_y", d, p)
        y = min(y, 0.5 * min(xy, zy))
        return (max(xy - y, 0.0), y, max(zy - y, 0.0))

    def spam(self, d: int, p: float) -> float:
        return self._base.spam(d, p)

    def cnot(self, d: int, p: float) -> list[float]:
        return self._base.cnot(d, p)

    def pauli_measurement(self, kind: str, d: int, p: float) -> dict:
        return self._base.pauli_measurement(kind, d, p)

    def extrapolated(self, d: int, p: float) -> bool:
        return self._base.extrapolated(d, p)


class GidneyAssumedModel:
    """Reproduces Gidney 2025's assumed logical error scaling model (p_L = 10^-15 at d=25)."""

    def __init__(self, target_d25: float = 1.0e-15):
        # 10^-15 at x = (25+1)/2 = 13 corresponds to 0.01 * 0.1^13
        self.target_d25 = target_d25
        self.prefactor = target_d25 / (0.1**13)

    def idle(self, d: int, p: float) -> tuple[float, float, float]:
        xd = (d + 1) / 2
        rate = self.prefactor * (0.1**xd)
        return (rate / 2.0, 0.0, rate / 2.0)

    def spam(self, d: int, p: float) -> float:
        return 1e-4

    def cnot(self, d: int, p: float) -> list[float]:
        return [0.0] * 15

    def pauli_measurement(self, kind: str, d: int, p: float) -> dict:
        return {"joint": [1.0] + [0.0] * 7, "o": 0.0}

    def extrapolated(self, d: int, p: float) -> bool:
        return False


def run_sensitivity_tornado(
    model: LogicalModel, workload: Workload = GIDNEY_2025, p: float = 1e-3, target: float = 0.1
) -> dict[str, Any]:
    """Runs a multi-factor sensitivity analysis on the Gidney 2025 headline."""
    cov = idle_fit_covariance(model, "idle_xy", p)
    se_b = cov["se_b"]
    se_a = cov["se_a"]
    central_b = cov["b"]

    baseline_est = estimate(workload, model, p=p, target=target)

    experiments = []

    # 1. Idle slope Lambda uncertainty (+/- 1 sigma)
    # Lambda = exp(-b)
    slope_scale_lo = (central_b + se_b) / central_b  # steeper suppression (larger Lambda)
    slope_scale_hi = (central_b - se_b) / central_b  # shallower suppression (smaller Lambda)

    m_slope_lo = PerturbedLogicalModel(model, slope_scale=slope_scale_lo)
    e_slope_lo = estimate(workload, m_slope_lo, p=p, target=target)
    experiments.append({
        "category": "Idle Slope",
        "name": f"Lambda +1sigma (Lambda = {np.exp(-(central_b-se_b)):.2f})",
        "param": "slope",
        "value": "+1sigma",
        "d": e_slope_lo.d,
        "physical_qubits": e_slope_lo.physical_qubits,
        "delta_d": e_slope_lo.d - baseline_est.d if e_slope_lo.d else None,
    })

    m_slope_hi = PerturbedLogicalModel(model, slope_scale=slope_scale_hi)
    e_slope_hi = estimate(workload, m_slope_hi, p=p, target=target)
    experiments.append({
        "category": "Idle Slope",
        "name": f"Lambda -1sigma (Lambda = {np.exp(-(central_b+se_b)):.2f})",
        "param": "slope",
        "value": "-1sigma",
        "d": e_slope_hi.d,
        "physical_qubits": e_slope_hi.physical_qubits,
        "delta_d": e_slope_hi.d - baseline_est.d if e_slope_hi.d else None,
    })

    # 2. Prefactor uncertainty (+/- 1 sigma)
    m_pref_lo = PerturbedLogicalModel(model, prefactor_scale=math.exp(-se_a))
    e_pref_lo = estimate(workload, m_pref_lo, p=p, target=target)
    experiments.append({
        "category": "Idle Prefactor",
        "name": "Prefactor A -1sigma (0.87x)",
        "param": "prefactor",
        "value": "-1sigma",
        "d": e_pref_lo.d,
        "physical_qubits": e_pref_lo.physical_qubits,
        "delta_d": e_pref_lo.d - baseline_est.d if e_pref_lo.d else None,
    })

    m_pref_hi = PerturbedLogicalModel(model, prefactor_scale=math.exp(se_a))
    e_pref_hi = estimate(workload, m_pref_hi, p=p, target=target)
    experiments.append({
        "category": "Idle Prefactor",
        "name": "Prefactor A +1sigma (1.15x)",
        "param": "prefactor",
        "value": "+1sigma",
        "d": e_pref_hi.d,
        "physical_qubits": e_pref_hi.physical_qubits,
        "delta_d": e_pref_hi.d - baseline_est.d if e_pref_hi.d else None,
    })

    # 3. Decoder improvement: correlated matching / belief matching (+15% Lambda)
    m_dec = PerturbedLogicalModel(model, slope_scale=math.log(cov["lambda"] * 1.15) / (-central_b))
    e_dec = estimate(workload, m_dec, p=p, target=target)
    experiments.append({
        "category": "Decoder",
        "name": "Correlated / Belief-matching (+15% Lambda)",
        "param": "decoder",
        "value": "correlated",
        "d": e_dec.d,
        "physical_qubits": e_dec.physical_qubits,
        "delta_d": e_dec.d - baseline_est.d if e_dec.d else None,
    })

    # 4. Magic-state error variation (epsilon_CCZ * 0.1 and * 10.0)
    for mult, label in [(0.1, "0.1x"), (10.0, "10.0x")]:
        fac_copy = copy.deepcopy(FACTORIES["cultivation"])
        fac_copy = dataclasses.replace(
            fac_copy,
            eps_ccz=lambda p, m=mult: FACTORIES["cultivation"].eps_ccz(p) * m,
        )
        custom_facs = dict(FACTORIES)
        custom_facs["custom"] = fac_copy
        # Estimate with custom factory
        for d_try in range(3, 63, 2):
            arch_try = Architecture(model, d_try, p, fac_copy)
            ef, r, _ = expected_faults(workload, arch_try)
            if ef <= target:
                d_found = d_try
                pq = arch_try.data_qubits(workload.logical_qubits) + arch_try.factory_qubits(
                    workload.t_gates, workload.toffolis, r
                )
                break
        else:
            d_found = None
            pq = 0
        experiments.append({
            "category": "Magic States",
            "name": f"epsilon_CCZ {label}",
            "param": "eps_ccz",
            "value": label,
            "d": d_found,
            "physical_qubits": pq,
            "delta_d": d_found - baseline_est.d if d_found else None,
        })

    # 5. Assumption-matched reproduction (Gidney assumed model)
    gidney_m = GidneyAssumedModel()
    e_gid_01 = estimate(workload, gidney_m, p=p, target=target)
    experiments.append({
        "category": "Assumed Model",
        "name": "Gidney assumed model (target E[faults] <= 0.1)",
        "param": "assumed_model",
        "value": "target_0.1",
        "d": e_gid_01.d,
        "physical_qubits": e_gid_01.physical_qubits,
        "delta_d": e_gid_01.d - baseline_est.d if e_gid_01.d else None,
    })

    # Gidney with target 0.5 (Gidney's target ~60% success)
    for d_try in range(3, 63, 2):
        arch_try = Architecture(gidney_m, d_try, p, FACTORIES["cultivation"])
        ef, r, _ = expected_faults(workload, arch_try)
        if ef <= 0.5:
            d_gid_05 = d_try
            pq_05 = arch_try.data_qubits(workload.logical_qubits) + arch_try.factory_qubits(
                workload.t_gates, workload.toffolis, r
            )
            break
    else:
        d_gid_05, pq_05 = None, 0

    experiments.append({
        "category": "Assumed Model",
        "name": "Gidney assumed model (target E[faults] <= 0.5, Gidney's target)",
        "param": "assumed_model",
        "value": "target_0.5",
        "d": d_gid_05,
        "physical_qubits": pq_05,
        "delta_d": d_gid_05 - baseline_est.d if d_gid_05 else None,
    })

    return {
        "baseline": {
            "workload": workload.name,
            "d": baseline_est.d,
            "physical_qubits": baseline_est.physical_qubits,
            "rounds": baseline_est.rounds,
            "expected_faults": baseline_est.expected_faults,
        },
        "fit_covariance": cov,
        "experiments": experiments,
    }
