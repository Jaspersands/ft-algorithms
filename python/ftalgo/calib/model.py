"""Logical channels from the calibration data, fitted in d and p.

Components (each a probability, measured at points (d, p) with a bootstrap uncertainty):

- ``idle_xy``, ``idle_zy``: per-round flip rates of a patch's Z and X logicals (X∪Y and Z∪Y
  errors), from the slope of memory experiments over T = d, 3d, 4d rounds; ``idle_y`` from the
  reference-qubit experiments; ``spam_z``/``spam_x`` from the memories' intercept (preparation
  plus readout).
- ``cnot_x{1,2,3}`` / ``cnot_z{1,2,3}``: the surgery CNOT's own channel, X part and Z part
  (pattern bit 0 = control, bit 1 = target), after dividing out a 4d-round memory baseline.
- ``zz_m``, ``zz_x{1,2,3}``, ``zz_z``: Z⊗Z measurement: wrong outcome, X errors on the two
  patches (joint), Z error per patch (from the X1X2 parity, patches assumed alike), after a
  3d-round baseline. ``xx_*`` the mirror image.

Each component c(d, p) is fitted per measured p as log c = a_p + b_p·(d + 1)/2 and globally as
c = A·(p/p*)^((d+1)/2); ``rate`` uses the per-p fit at measured p (global otherwise) and reports
whether (d, p) lies outside the directly measured range.
"""

from __future__ import annotations

import glob
import json
import math
import os
from dataclasses import dataclass, field

import numpy as np

from . import channels as ch
from .qec import Counts

DATA = os.path.join(os.path.dirname(__file__), "..", "..", "..", "data", "calibration")
MIN_EVENTS = 10
BOOT = 300


def load_raw(path: str | None = None) -> dict[str, dict]:
    path = path or os.path.join(DATA, "raw")
    out = {}
    for f in glob.glob(os.path.join(path, "*.json")):
        with open(f) as fh:
            r = json.load(fh)
        out[r["name"]] = r
    return out


def _probs(rec: dict) -> tuple[np.ndarray, int]:
    c = Counts.from_json(rec["counts"])
    return c.probabilities(), c.shots


def _boot(rng, probs: np.ndarray, shots: int) -> np.ndarray:
    return rng.multinomial(shots, probs / probs.sum()) / shots


@dataclass
class Point:
    d: int
    p: float
    value: float
    sigma: float
    events: float  # effective number of events behind the value

    def ok(self) -> bool:
        return self.events >= MIN_EVENTS and self.value > 0


def _slope(fids: list[tuple[int, np.ndarray]]) -> tuple[np.ndarray, np.ndarray]:
    """Least-squares fit of log f(T) = a + b·T per fidelity entry; returns (per-round f, intercept f)."""
    Ts = np.array([t for t, _ in fids], dtype=float)
    F = np.array([f for _, f in fids])
    F = np.clip(F, 1e-300, None)
    A = np.vstack([np.ones_like(Ts), Ts]).T
    coef, *_ = np.linalg.lstsq(A, np.log(F), rcond=None)
    return np.exp(coef[1]), np.exp(coef[0])


class Components:
    """Derives component values (with bootstrap σ) at every measured (d, p)."""

    def __init__(self, raw: dict[str, dict], seed: int = 0):
        self.raw = raw
        self.rng = np.random.default_rng(seed)
        self.points: dict[str, list[Point]] = {}

    def _add(self, name: str, d: int, p: float, samples: list[float], central: float, events: float) -> None:
        s = float(np.std(samples)) if samples else float("nan")
        self.points.setdefault(name, []).append(Point(d, p, float(central), s, float(events)))

    def _grid(self):
        dps = sorted({(r["d"], r["p"]) for r in self.raw.values()})
        return dps

    def build(self) -> "Components":
        for d, p in self._grid():
            self._memory(d, p)
            self._choi(d, p)
            self._cnot(d, p)
            self._zz(d, p, "zz")
            self._zz(d, p, "xx")
        return self

    # -- memory --------------------------------------------------------------------------------
    def _mem_recs(self, b, d, p):
        out = []
        for T in (d, 3 * d, 4 * d):
            r = self.raw.get(f"mem-{b}-d{d}-p{p}-T{T}")
            if r:
                out.append((T, r))
        return out

    def memory_round(self, b: str, d: int, p: float, boot: bool = False) -> tuple[float, float] | None:
        """(per-round flip q, spam flip q) for basis b; with boot=True from a resample."""
        recs = self._mem_recs(b, d, p)
        if len(recs) < 2:
            return None
        fids = []
        for T, r in recs:
            pr, n = _probs(r)
            if boot:
                pr = _boot(self.rng, pr, n)
            fids.append((T, ch.fidelities(pr)))
        fr, f0 = _slope([(T, f[1:]) for T, f in fids])
        return (1 - fr[0]) / 2, (1 - f0[0]) / 2

    def _memory(self, d, p):
        for b, name in (("z", "idle_xy"), ("x", "idle_zy")):
            c = self.memory_round(b, d, p)
            if c is None:
                continue
            samples = [self.memory_round(b, d, p, boot=True) for _ in range(BOOT)]
            events = sum(Counts.from_json(r["counts"]).failures for _, r in self._mem_recs(b, d, p))
            self._add(name, d, p, [s[0] for s in samples], c[0], events)
            self._add(f"spam_{b}", d, p, [s[1] for s in samples], c[1], events)

    # -- reference-qubit memory ----------------------------------------------------------------
    def choi_round(self, d, p, boot=False):
        a, b = self.raw.get(f"choi-d{d}-p{p}-T{d}"), self.raw.get(f"choi-d{d}-p{p}-T{3 * d}")
        if not (a and b):
            return None
        fs = []
        for r in (a, b):
            pr, n = _probs(r)
            if boot:
                pr = _boot(self.rng, pr, n)
            fs.append(ch.fidelities(pr))
        pr = ch.probabilities(ch.per_round(fs[0], fs[1], d, 3 * d), clip=False)
        return pr  # [I, Z, X, Y] (bit 0 = obs0 = Z∪Y, bit 1 = obs1 = X∪Y)

    def _choi(self, d, p):
        c = self.choi_round(d, p)
        if c is None:
            return
        samples = [self.choi_round(d, p, boot=True) for _ in range(BOOT)]
        a, b = self.raw[f"choi-d{d}-p{p}-T{d}"], self.raw[f"choi-d{d}-p{p}-T{3 * d}"]
        ev = sum(Counts.from_json(r["counts"]).patterns.get(3, 0) for r in (a, b))
        self._add("idle_y", d, p, [s[3] for s in samples], c[3], ev)

    # -- baselines ------------------------------------------------------------------------------
    def half(self, b: str, d: int, p: float, boot=False) -> float | None:
        """Flip probability of half of a 2d-round memory (preparation + d rounds, or d rounds +
        readout; taken equal by time symmetry): f_half = √f_mem(2d)."""
        r = self.raw.get(f"mem-{b}-d{d}-p{p}-T{2 * d}")
        if r is not None:
            pr, n = _probs(r)
            if boot:
                pr = _boot(self.rng, pr, n)
            f2 = 1 - 2 * pr[1]
        else:
            a, c = self.raw.get(f"mem-{b}-d{d}-p{p}-T{d}"), self.raw.get(f"mem-{b}-d{d}-p{p}-T{3 * d}")
            if not (a and c):
                return None
            fs = []
            for rr in (a, c):
                pr, n = _probs(rr)
                if boot:
                    pr = _boot(self.rng, pr, n)
                fs.append(1 - 2 * pr[1])
            f2 = math.sqrt(max(fs[0] * fs[1], 1e-300))
        return (1 - math.sqrt(max(f2, 0.0))) / 2

    # -- CNOT ----------------------------------------------------------------------------------
    def cnot_part(self, inp: str, d: int, p: float, boot=False):
        """The CNOT's channel over its 2d merge rounds (idling included), X part for inputs "z"
        and Z part for "x"; bit 0 = control, bit 1 = target. Errors of the d rounds before the
        gate are pushed through it (X: control → both; Z: target → both) before dividing out."""
        r = self.raw.get(f"cnot-{inp}-d{d}-p{p}")
        qh = self.half(inp, d, p, boot)
        if r is None or qh is None:
            return None
        pr, n = _probs(r)
        if boot:
            pr = _boot(self.rng, pr, n)
        indep = ch.probabilities(ch.product_fidelities(qh, qh), clip=False)
        if inp == "z":
            q = pr  # obs0 = Z_C ↔ X on C, obs1 = Z_T (+frame) ↔ X on T
            pre = ch.remap(indep, [0b01, 0b11])
        else:
            q = ch.remap(pr, [0b11, 0b01])  # (o0 = X_T, o1 = X_C X_T) → (z_C, z_T)
            pre = ch.remap(indep, [0b11, 0b10])
        f_op = ch.fidelities(q) / (ch.fidelities(pre) * ch.product_fidelities(qh, qh))
        return ch.probabilities(f_op, clip=False)

    def _cnot(self, d, p):
        for inp, tag in (("z", "x"), ("x", "z")):
            c = self.cnot_part(inp, d, p)
            if c is None:
                continue
            samples = [self.cnot_part(inp, d, p, boot=True) for _ in range(BOOT)]
            cnt = Counts.from_json(self.raw[f"cnot-{inp}-d{d}-p{p}"]["counts"])
            pats = cnt.patterns if inp == "z" else {((k & 1) ^ (k >> 1 & 1)) | ((k & 1) << 1): v for k, v in cnt.patterns.items()}
            for k in (1, 2, 3):
                self._add(f"cnot_{tag}{k}", d, p, [s[k] for s in samples], c[k], pats.get(k, 0))

    # -- Z⊗Z / X⊗X measurement -----------------------------------------------------------------
    def zz_parts(self, kind: str, d: int, p: float, boot=False):
        """Z⊗Z ("zz") or X⊗X ("xx") over its d merge rounds: the joint law of (wrong outcome,
        measured-type error on patch 1, on patch 2) as an ideal measurement followed by those
        flips, and the other-type error per patch (from the preserved parity, patches alike)."""
        same = "z" if kind == "zz" else "x"
        other = "x" if kind == "zz" else "z"
        r1 = self.raw.get(f"{kind}-{same}-d{d}-p{p}")
        r2 = self.raw.get(f"{kind}-{other}-d{d}-p{p}")
        qs, qo = self.half(same, d, p, boot), self.half(other, d, p, boot)
        if not (r1 and r2) or qs is None or qo is None:
            return None
        p1, n1 = _probs(r1)
        p2, n2 = _probs(r2)
        if boot:
            p1, p2 = _boot(self.rng, p1, n1), _boot(self.rng, p2, n2)
        indep = ch.probabilities(ch.product_fidelities(qs, qs), clip=False)
        pre = ch.remap(indep, [0b11, 0b01, 0b10])  # an error before the merge changes the outcome too
        f_op = ch.fidelities(p1) / (ch.fidelities(pre) * ch.product_fidelities(0.0, qs, qs))
        joint = ch.probabilities(f_op, clip=False)
        fo = 1 - 2 * qo
        fpar = (1 - 2 * p2[1]) / fo**4
        q_other = (1 - math.sqrt(max(fpar, 0.0))) / 2
        return joint, q_other

    def _zz(self, d, p, kind):
        c = self.zz_parts(kind, d, p)
        if c is None:
            return
        samples = [self.zz_parts(kind, d, p, boot=True) for _ in range(BOOT)]
        same = "z" if kind == "zz" else "x"
        other = "x" if kind == "zz" else "z"
        n1 = Counts.from_json(self.raw[f"{kind}-{same}-d{d}-p{p}"]["counts"])
        n2 = Counts.from_json(self.raw[f"{kind}-{other}-d{d}-p{p}"]["counts"])
        for k in range(1, 8):
            self._add(f"{kind}_j{k}", d, p, [s[0][k] for s in samples], c[0][k], n1.patterns.get(k, 0))
        self._add(f"{kind}_o", d, p, [s[1] for s in samples], c[1], n2.failures)


@dataclass
class Fit:
    """log c = a + b·(d+1)/2 at one p (or globally with b = log(p/p*))."""

    a: float
    b: float
    dmax: int  # largest d with a usable point
    npts: int

    def __call__(self, d: int) -> float:
        return math.exp(self.a + self.b * (d + 1) / 2)


@dataclass
class GlobalFit:
    logA: float
    logpstar: float

    def __call__(self, d: int, p: float) -> float:
        return math.exp(self.logA + (d + 1) / 2 * (math.log(p) - self.logpstar))


def _wls(x: np.ndarray, y: np.ndarray, w: np.ndarray) -> tuple[float, float]:
    A = np.vstack([np.ones_like(x), x]).T * np.sqrt(w)[:, None]
    coef, *_ = np.linalg.lstsq(A, y * np.sqrt(w), rcond=None)
    return float(coef[0]), float(coef[1])


@dataclass
class LogicalModel:
    points: dict[str, list[Point]]
    per_p: dict[str, dict[float, Fit]] = field(default_factory=dict)
    glob: dict[str, GlobalFit] = field(default_factory=dict)

    @staticmethod
    def from_raw(raw: dict[str, dict] | None = None) -> "LogicalModel":
        comps = Components(raw if raw is not None else load_raw()).build()
        m = LogicalModel(comps.points)
        m.fit()
        return m

    def fit(self) -> None:
        for name, pts in self.points.items():
            good = [q for q in pts if q.ok()]
            self.per_p[name] = {}
            for p in sorted({q.p for q in good}):
                g = [q for q in good if q.p == p]
                if len(g) >= 2:
                    x = np.array([(q.d + 1) / 2 for q in g])
                    y = np.log([q.value for q in g])
                    w = np.array([q.events for q in g], dtype=float)
                    a, b = _wls(x, y, w)
                    self.per_p[name][p] = Fit(a, min(b, 0.0), max(q.d for q in g), len(g))
            if len({q.p for q in good}) >= 2 and len(good) >= 3:
                x = np.array([(q.d + 1) / 2 for q in good])
                y = np.log([q.value for q in good]) - x * np.log([q.p for q in good])
                w = np.array([q.events for q in good], dtype=float)
                a, b = _wls(x, y, w)
                self.glob[name] = GlobalFit(a, -b)

    def rate(self, name: str, d: int, p: float) -> tuple[float, bool]:
        """(value, extrapolated?)."""
        f = self.per_p.get(name, {}).get(p)
        if f is not None:
            return f(d), d > f.dmax
        g = self.glob.get(name)
        if g is None:
            return 0.0, True
        return g(d, p), True

    # -- channels for the noise compiler --------------------------------------------------------
    def idle(self, d: int, p: float) -> tuple[float, float, float]:
        """Per-round (p_X, p_Y, p_Z) of one patch."""
        xy, _ = self.rate("idle_xy", d, p)
        zy, _ = self.rate("idle_zy", d, p)
        y, _ = self.rate("idle_y", d, p)
        y = min(y, 0.5 * min(xy, zy))
        return (max(xy - y, 0.0), y, max(zy - y, 0.0))

    def spam(self, d: int, p: float) -> float:
        """Flip probability of a transversal preparation-plus-readout (mean of the two bases)."""
        return 0.5 * (self.rate("spam_z", d, p)[0] + self.rate("spam_x", d, p)[0])

    def cnot(self, d: int, p: float) -> list[float]:
        """The CNOT's own (beyond idling) channel as PAULI_CHANNEL_2 probabilities (control first)."""
        xp = np.array([0.0] + [self.rate(f"cnot_x{k}", d, p)[0] for k in (1, 2, 3)])
        zp = np.array([0.0] + [self.rate(f"cnot_z{k}", d, p)[0] for k in (1, 2, 3)])
        xp[0], zp[0] = 1 - xp[1:].sum(), 1 - zp[1:].sum()
        return ch.pauli_channel_2(xp, zp)

    def pauli_measurement(self, kind: str, d: int, p: float) -> dict:
        """Z⊗Z (kind "zz") or X⊗X ("xx") over its d rounds: "joint" = probabilities of the 8
        patterns (bit 0 wrong outcome, bit 1/2 measured-type error on patch 1/2) applied after an
        ideal measurement; "o" = the other-type error per patch."""
        j = np.array([0.0] + [self.rate(f"{kind}_j{k}", d, p)[0] for k in range(1, 8)])
        j[0] = 1 - j[1:].sum()
        return {"joint": j, "o": self.rate(f"{kind}_o", d, p)[0]}

    def extrapolated(self, d: int, p: float) -> bool:
        return any(self.rate(n, d, p)[1] for n in ("idle_xy", "idle_zy", "cnot_x3", "zz_j1"))

    def summary(self) -> dict:
        out = {}
        for name, pts in sorted(self.points.items()):
            out[name] = {
                "points": [{"d": q.d, "p": q.p, "value": q.value, "sigma": q.sigma, "events": q.events} for q in sorted(pts, key=lambda q: (q.p, q.d))],
                "per_p": {str(p): {"a": f.a, "b": f.b, "lambda": math.exp(-f.b), "dmax": f.dmax} for p, f in self.per_p.get(name, {}).items()},
                "global": ({"A": math.exp(self.glob[name].logA), "pstar": math.exp(self.glob[name].logpstar)} if name in self.glob else None),
            }
        return out
