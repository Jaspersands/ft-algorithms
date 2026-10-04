"""Success probabilities with error bars, stratified by the number of faults.

With K the number of static faults, P(success) = Σ_k P(K = k)·S_k, where P(K = k) is exact
(``Program.fault_count_distribution``) and S_k is estimated from shots with exactly k faults
(dynamic faults, inside TABLE cases, fire at their own rates in every stratum). Strata beyond
kmax carry the tail mass T; they are assumed to behave like the last sampled stratum, and the
reported interval widens by T on both sides. When the expected number of faults is large,
stratification has no advantage and plain sampling is used instead.

The k = 1 stratum, with faults logged, gives the harm of a single fault in each operation class:
the error budget is the expected number of faults per class times that harm.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable

import numpy as np

from .engine import Program

PLAIN_ABOVE = 12.0  # expected faults above which plain sampling is used


@dataclass
class Estimate:
    """value and sigma are floats for a single score, lists for several (score columns)."""

    value: float
    sigma: float
    tail: float
    method: str
    strata: list = field(default_factory=list)
    single_fault: dict = field(default_factory=dict)  # class → (shots, mean score)
    shots: int = 0

    def interval(self, z: float = 2.0, i: int = 0) -> tuple[float, float]:
        v = self.value if np.isscalar(self.value) else self.value[i]
        s = self.sigma if np.isscalar(self.sigma) else self.sigma[i]
        return (max(v - z * s - self.tail, 0.0), min(v + z * s + self.tail, 1.0))

    def to_json(self) -> dict:
        return {
            "value": self.value,
            "sigma": self.sigma,
            "tail": self.tail,
            "method": self.method,
            "strata": self.strata,
            "single_fault": self.single_fault,
            "shots": self.shots,
        }


def stratified(
    prog: Program,
    score: Callable[[np.ndarray], np.ndarray],
    *,
    seed: int,
    shots: int = 4000,
    min_shots: int = 200,
    tail: float = 1e-4,
    kmax_limit: int = 60,
    threads: int = 0,
    attribute: bool = True,
) -> Estimate:
    """score(rows) → per-shot scores in [0, 1]: a vector (one score) or a (shots × S) array."""

    def sc(rows):
        s = np.asarray(score(rows), dtype=float)
        return s[:, None] if s.ndim == 1 else s

    def out(v):
        v = [float(x) for x in np.atleast_1d(v)]
        return v[0] if single else v

    ef = prog.expected_faults
    raw = np.asarray(score(prog.sample(1, seed=seed, faults="none", threads=1)))
    single = raw.ndim == 1
    S = 1 if single else raw.shape[1]
    if ef > PLAIN_ABOVE:
        rows = prog.sample(shots, seed=seed, threads=threads)
        s = sc(rows)
        return Estimate(out(s.mean(axis=0)), out(s.std(axis=0, ddof=1) / math.sqrt(len(s))), 0.0, "plain", shots=shots)
    K = 4
    while True:
        pk = prog.fault_count_distribution(K)
        if pk[-1] < tail or K >= kmax_limit:
            break
        K = min(2 * K, kmax_limit)
    # Drop strata that carry negligible weight from the top.
    while K > 0 and pk[K] < tail / 10:
        pk[-1] += pk[K]
        K -= 1
    pk = np.concatenate([pk[: K + 1], [pk[-1]]])
    value = np.zeros(S)
    var = np.zeros(S)
    strata = []
    per_class: dict = {}
    total = 0
    last = None
    classes = prog.classes
    for k in range(K + 1):
        w = float(pk[k])
        if w <= 0 or k > prog.num_sites:
            continue
        n = int(min(max(round(shots * w), min_shots), shots))
        if k == 1 and attribute:
            rows, faults = prog.sample(n, seed=seed * 1000 + k, faults=k, threads=threads, log_faults=True)
        else:
            rows = prog.sample(n, seed=seed * 1000 + k, faults=k, threads=threads)
        s = sc(rows)
        sk = s.mean(axis=0)
        vk = s.var(axis=0, ddof=1) / n if n > 1 else np.full(S, 0.25)
        value += w * sk
        var += w * w * vk
        total += n
        last = sk
        strata.append({"k": k, "weight": w, "shots": n, "score": out(sk)})
        if k == 1 and attribute:
            static = faults.site != faults.DYNAMIC
            shot_cls = np.full(n, -1)
            shot_cls[faults.shot[static].astype(np.int64)] = faults.cls[static]
            for c in np.unique(shot_cls):
                if c < 0:
                    continue
                sel = shot_cls == c
                per_class[classes[c]] = (int(sel.sum()), out(s[sel].mean(axis=0)))
    t = float(pk[-1])
    if last is not None:
        value += t * last
    return Estimate(out(value), out(np.sqrt(var)), t, "stratified", strata, per_class, total)


def error_budget(budget: dict, est: Estimate, i: int = 0) -> dict:
    """Per class: expected faults, the score (column i) of shots with one fault there, and the
    expected score lost (expected faults × (S₀ − S₁(class)))."""
    pick = lambda v: v if np.isscalar(v) else v[i]  # noqa: E731
    s0 = next((pick(s["score"]) for s in est.strata if s["k"] == 0), None)
    out = {}
    for cls, ef in budget.items():
        sf = est.single_fault.get(cls)
        if s0 is None or sf is None:
            out[cls] = {"expected_faults": ef, "single_fault_score": None, "loss": None}
        else:
            out[cls] = {"expected_faults": ef, "single_fault_score": pick(sf[1]), "single_fault_shots": sf[0], "loss": ef * (s0 - pick(sf[1]))}
    return out
