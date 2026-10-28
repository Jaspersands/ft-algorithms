"""How firmly the measured model sets the code distance of Gidney 2025's RSA-2048 algorithm.

Gidney (arXiv:2505.15917, Section 3.2) picks d = 25 because his projection of surface-code memory
errs at most 1e-15 per patch-round there, which over one 12-hour shot leaves 93.3% of shots free of
logical errors. This module asks the same two questions of the measured model: the smallest d
whose idle patch meets 1e-15 per round (``scale.distance_for_idle``), and the smallest d for which
one shot of his algorithm, run through ``scale.estimate``, meets his budget. It then checks how the
second answer moves with

1. the idle fits' uncertainty (weighted least squares with covariance; both idle components moved
   to the edge of their k-sigma band together);
2. a steeper extrapolation beyond the measured distances, as a better decoder would give;
3. the CCZ states' error rate;
4. the failure budget;

and that the same estimator, given his projection in place of the measured model, returns his d.
"""

from __future__ import annotations

import dataclasses
import math
from typing import Any

import numpy as np

from .arch import FACTORIES
from .calib.model import LogicalModel
from .scale import (GIDNEY_2025, GIDNEY_2025_D, GIDNEY_2025_PER_ROUND, GIDNEY_2025_TARGET, Workload,
                    distance_for_idle, estimate, idle_total)

IDLE = ("idle_xy", "idle_zy")


def fit_covariance(model: LogicalModel, name: str, p: float = 1e-3) -> dict[str, Any]:
    """The per-p fit log c = a + b·(d+1)/2 of one component with its parameter covariance (weighted
    least squares, weights = events, scaled by the reduced χ²). Reproduces the model's own fit."""
    pts = [q for q in model.points.get(name, []) if q.ok() and q.p == p]
    if len(pts) < 3:
        raise ValueError(f"need at least three usable points for {name} at p = {p}")
    x = np.array([(q.d + 1) / 2 for q in pts])
    y = np.log([q.value for q in pts])
    w = np.array([q.events for q in pts], dtype=float)
    A = np.vstack([np.ones_like(x), x]).T
    inv = np.linalg.inv(A.T @ (w[:, None] * A))
    coef = inv @ (A.T @ (w * y))
    res = y - A @ coef
    cov = inv * float((res**2 * w).sum() / (len(pts) - 2))
    f = model.per_p[name][p]
    assert abs(coef[0] - f.a) < 1e-9 and abs(coef[1] - f.b) < 1e-9, "covariance fit differs from the model's"
    b, se_b = float(coef[1]), float(math.sqrt(cov[1, 1]))
    return {
        "name": name, "p": p, "a": float(coef[0]), "b": b, "cov": cov.tolist(), "dmax": f.dmax, "npts": len(pts),
        "lambda": math.exp(-b), "lambda_lo": math.exp(-(b + se_b)), "lambda_hi": math.exp(-(b - se_b)),
    }


def log_sigma(fc: dict[str, Any], d: int) -> float:
    """Standard error of the fitted log rate at distance d."""
    x = (d + 1) / 2
    c = fc["cov"]
    return math.sqrt(max(c[0][0] + 2 * x * c[0][1] + x * x * c[1][1], 0.0))


class _Derived(LogicalModel):
    """The measured model with its idle rates replaced by ``_idle_rate``."""

    def __init__(self, base: LogicalModel):
        super().__init__(base.points, base.per_p, base.glob)

    def rate(self, name: str, d: int, p: float) -> tuple[float, bool]:
        v, ex = super().rate(name, d, p)
        if name in IDLE:
            return self._idle_rate(name, d, p, v), ex
        return v, ex

    def _idle_rate(self, name: str, d: int, p: float, v: float) -> float:
        raise NotImplementedError


class IdleBand(_Derived):
    """Both idle components moved k standard errors of their fitted log rate (covariance included);
    the logical-Y component, about 2% of the total, stays central."""

    def __init__(self, base: LogicalModel, k: float, p: float = 1e-3):
        super().__init__(base)
        self.k = k
        self.fc = {n: fit_covariance(base, n, p) for n in IDLE}

    def _idle_rate(self, name, d, p, v):
        fc = self.fc[name]
        return v * math.exp(self.k * log_sigma(fc, d)) if p == fc["p"] else v


class SteeperTail(_Derived):
    """Beyond the largest measured distance the idle rates fall by factor·Λ per two units of
    distance instead of Λ; every measured point is kept."""

    def __init__(self, base: LogicalModel, factor: float):
        super().__init__(base)
        self.factor = factor

    def _idle_rate(self, name, d, p, v):
        f = self.per_p.get(name, {}).get(p)
        if f is None or d <= f.dmax:
            return v
        return f(f.dmax) * (math.exp(-f.b) * self.factor) ** (-(d - f.dmax) / 2)


class GidneyProjection:
    """Gidney 2025's projection of surface-code memory at p = 1e-3 (the legend of his Figure 6:
    p_L ≈ r·n·3.5^-d / 40 for r rounds of n patches; it reaches 1e-15 per patch-round first at d = 25,
    as his Section 3.2 says), split evenly between X and Z. Nothing but idling errs: his estimate
    budgets only storage error."""

    GAIN = 3.5  # per unit of distance, so Λ = 3.5² = 12.25 per two units

    @staticmethod
    def rate_per_round(d: int) -> float:
        return GidneyProjection.GAIN ** -d / 40

    def idle(self, d: int, p: float) -> tuple[float, float, float]:
        r = self.rate_per_round(d)
        return (r / 2, 0.0, r / 2)

    def spam(self, d: int, p: float) -> float:
        return 0.0

    def cnot(self, d: int, p: float) -> list[float]:
        return [0.0] * 15

    def pauli_measurement(self, kind: str, d: int, p: float) -> dict:
        return {"joint": np.array([1.0] + [0.0] * 7), "o": 0.0}

    def extrapolated(self, d: int, p: float) -> bool:
        return False


class FittedMemory(GidneyProjection):
    """Idle error from per-round fits of a memory's two logicals (data/calibration/noise_models.json):
    a Z-basis memory's flips are X or Y errors, an X-basis memory's Z or Y. Y is counted in both,
    which over-counts by the Y share, about 2%. Like GidneyProjection, nothing but idling errs."""

    def __init__(self, fits: dict):
        self.fits = fits

    def basis_rate(self, basis: str, d: int) -> float:
        f = self.fits[basis]
        return math.exp(f["a"] + f["b"] * (d + 1) / 2)

    def idle(self, d: int, p: float) -> tuple[float, float, float]:
        return (self.basis_rate("z", d), 0.0, self.basis_rate("x", d))


def noise_models(workload: Workload = GIDNEY_2025, path: str | None = None) -> dict[str, Any] | None:
    """The memory models of tools/noise_models.py (SD6 as calibrated, Stim's uniform depolarizing,
    SI1000) against Gidney's projection: each model's distance for 1e-15 per patch-round, counting
    the worse basis alone or either basis, and for one shot of his algorithm at his budget."""
    import json
    import os

    from .calib.model import DATA

    path = path or os.path.join(DATA, "noise_models.json")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        nm = json.load(f)
    out = {}
    for name, m in nm["models"].items():
        if set(m["fits"]) != {"z", "x"}:
            continue
        fm = FittedMemory(m["fits"])
        d_basis = next((d for d in range(3, 62, 2) if max(fm.basis_rate("z", d), fm.basis_rate("x", d)) <= GIDNEY_2025_PER_ROUND), None)
        measured = [{"d": r["d"], "z": r["z"]["q"], "x": r["x"]["q"], "gidney": GidneyProjection.rate_per_round(r["d"])} for r in m["rows"]]
        out[name] = {
            "lambda": {b: m["fits"][b]["lambda"] for b in ("z", "x")},
            "dmax": min(m["fits"][b]["dmax"] for b in ("z", "x")),
            "measured": measured,
            "fit": {str(d): {"z": fm.basis_rate("z", d), "x": fm.basis_rate("x", d)} for d in range(3, 37, 2)},
            "d_basis": d_basis,
            "d_any": distance_for_idle(fm, GIDNEY_2025_PER_ROUND),
            "d_shot": _est(fm, workload, nm["p"], GIDNEY_2025_TARGET)["d"],
        }
    return out


def _est(model, w: Workload, p: float, target: float, factory=FACTORIES["cultivation"]) -> dict[str, Any]:
    e = estimate(w, model, p=p, factory=factory, target=target)
    return {"d": e.d, "expected_faults": e.expected_faults, "breakdown": e.breakdown, "extrapolated": e.extrapolated}


def run(model: LogicalModel, workload: Workload = GIDNEY_2025, p: float = 1e-3) -> dict[str, Any]:
    target = GIDNEY_2025_TARGET
    fits = {n: fit_covariance(model, n, p) for n in IDLE}

    rows = []
    for d in range(3, 37, 2):
        bands = {k: idle_total(IdleBand(model, k, p), d, p) for k in (-2, -1, 1, 2)}
        rows.append({"d": d, "measured": idle_total(model, d, p), "lo1": bands[-1], "hi1": bands[1], "lo2": bands[-2],
                     "hi2": bands[2], "gidney": GidneyProjection.rate_per_round(d),
                     "extrapolated": any(d > fits[n]["dmax"] for n in IDLE)})
    per_round = {
        "target": GIDNEY_2025_PER_ROUND,
        "rows": rows,
        "d_measured": distance_for_idle(model, GIDNEY_2025_PER_ROUND, p),
        "d_band": {str(k): distance_for_idle(IdleBand(model, k, p), GIDNEY_2025_PER_ROUND, p) for k in (-2, -1, 1, 2)},
        "d_gidney": distance_for_idle(GidneyProjection(), GIDNEY_2025_PER_ROUND, p),
    }

    base = _est(model, workload, p, target)
    scenarios = []

    def add(category, name, r):
        scenarios.append({"category": category, "name": name, "d": r["d"], "expected_faults": r["expected_faults"],
                          "delta_d": None if r["d"] is None or base["d"] is None else r["d"] - base["d"]})

    for k in (-2, -1, 1, 2):
        add("Idle fit", f"rates at the {'+' if k > 0 else '−'}{abs(k)}σ edge of the fit", _est(IdleBand(model, k, p), workload, p, target))
    add("Extrapolation", "Λ × 1.15 beyond the measured distances", _est(SteeperTail(model, 1.15), workload, p, target))
    cult = FACTORIES["cultivation"]
    for mult in (0.1, 10.0):
        fac = dataclasses.replace(cult, eps_ccz=lambda q, m=mult: cult.eps_ccz(q) * m)
        add("Magic states", f"ε_CCZ × {mult:g}", _est(model, workload, p, target, fac))
    for t, label in ((0.1, "90% of shots succeed (E[faults] ≤ 0.1)"), (0.01, "99% of shots succeed (E[faults] ≤ 0.01)")):
        add("Budget", label, _est(model, workload, p, t))

    return {
        "p": p,
        "workload": dataclasses.asdict(workload),
        "published_d": GIDNEY_2025_D,
        "gidney_lambda": GidneyProjection.GAIN ** 2,
        "target": target,
        "fits": fits,
        "per_round": per_round,
        "baseline": base,
        "gidney_projection": _est(GidneyProjection(), workload, p, target),
        "scenarios": scenarios,
        "noise_models": noise_models(workload),
    }
