"""Composition checks: whole lattice-surgery experiments predicted from their calibrated parts.

Each experiment is rebuilt as a *logical* program for ftsim: ideal Pauli-product measurements
(through an ancilla), with the calibrated channels in between: the half-memory channels before
and after (preparation + d rounds, d rounds + readout), the Z⊗Z / X⊗X joint channels of each
merge, and idling of patches not in a merge. Its predicted observable flips are compared with
the circuit-level measurement of the same experiment. Agreement says the per-operation channels
compose; the size of any disagreement is the composition error the algorithm results inherit.
"""

from __future__ import annotations

import math

import numpy as np

from ..engine import Program
from .model import Components
from .qec import Counts


def _f(x) -> str:
    return repr(float(max(x, 0.0)))


def _p1(q):
    px, py, pz = q
    return f"PAULI_CHANNEL_1({_f(px)},{_f(py)},{_f(pz)})"


def _z3(j):
    return "Z_CHANNEL_3(" + ",".join(_f(v) for v in j[1:]) + ")"


def _idle_pow(q, n):
    """Per-round (px, py, pz) composed n times, via fidelities."""
    px, py, pz = q
    fx = 1 - 2 * (py + pz)  # fidelity of X-type observables... per Pauli: f_X = 1 − 2(p_Y + p_Z)
    fy = 1 - 2 * (px + pz)
    fz = 1 - 2 * (px + py)
    fx, fy, fz = fx**n, fy**n, fz**n
    return ((1 + fx - fy - fz) / 4, (1 - fx + fy - fz) / 4, (1 - fx - fy + fz) / 4)


class Parts:
    """The measured components at one (d, p), straight from the data (no fit)."""

    def __init__(self, comps: Components, d: int, p: float):
        self.d, self.p = d, p
        self.hz = comps.half("z", d, p)  # X errors seen by Z-type logicals
        self.hx = comps.half("x", d, p)
        mz, mx = comps.memory_round("z", d, p), comps.memory_round("x", d, p)
        self.idle = (mz[0], 0.0, mx[0])  # per round: X∪Y as X, Z∪Y as Z (Y is small; see Choi)
        self.zz = comps.zz_parts("zz", d, p)
        self.xx = comps.zz_parts("xx", d, p)

    def ok(self) -> bool:
        return None not in (self.hz, self.hx, self.zz, self.xx)


def cnot_program(parts: Parts, inputs: str) -> str:
    """The surgery CNOT experiment as a logical program (C=0, A=1, T=2, ancilla 3); records:
    m1 (Z_C Z_A), m2 (X_A X_T), m3 (Z_A), then C and T read in the input basis."""
    d = parts.d
    idle_d = _idle_pow(parts.idle, d)
    zz_j, zz_o = parts.zz
    xx_j, xx_o = parts.xx
    L = []
    if inputs == "z":
        L += ["R 0 2", f"X_ERROR({_f(parts.hz)}) 0 2"]
    else:
        L += ["RX 0 2", f"Z_ERROR({_f(parts.hx)}) 0 2"]
    L += ["RX 1", f"Z_ERROR({_f(parts.hx)}) 1"]  # the ancilla patch: |+⟩ and d rounds
    # Z_C Z_A, d rounds; T idles.
    L += ["R 3", "CX 0 3 1 3", "H 3 0 1", f"{_z3(zz_j)} 3 0 1", "H 3 0 1", "M 3", f"Z_ERROR({_f(zz_o)}) 0 1", f"{_p1(idle_d)} 2"]
    # X_A X_T, d rounds; C idles.
    L += ["R 3", "H 3", "CX 3 1 3 2", f"{_z3(xx_j)} 3 1 2", "H 3", "M 3", f"X_ERROR({_f(xx_o)}) 1 2", f"{_p1(idle_d)} 0"]
    L += [f"{_p1(parts.idle)} 1", "M 1"]  # A read out in Z (one round)
    if inputs == "z":
        L += [f"X_ERROR({_f(parts.hz)}) 0 2", "M 0 2"]
    else:
        L += [f"Z_ERROR({_f(parts.hx)}) 0 2", "MX 0 2"]
    return "\n".join(L) + "\n"


def cnot_flips(rows: np.ndarray, inputs: str) -> np.ndarray:
    m1, m2, m3, c, t = (rows[:, i] for i in range(5))
    if inputs == "z":  # obs0 = Z_C, obs1 = Z_T ⊕ m1 ⊕ m3
        return np.stack([c, t ^ m1 ^ m3], axis=1)
    return np.stack([t, c ^ t ^ m2], axis=1)  # obs0 = X_T, obs1 = X_C X_T ⊕ m2


def repeated_program(parts: Parts, k: int) -> str:
    """k Z⊗Z measurements in a row on |00⟩, one round apart; records: outcomes, then Z1, Z2."""
    zz_j, zz_o = parts.zz
    L = ["R 0 1", f"X_ERROR({_f(parts.hz)}) 0 1"]
    for _ in range(k):
        L += ["R 3", "CX 0 3 1 3", "H 3 0 1", f"{_z3(zz_j)} 3 0 1", "H 3 0 1", "M 3", f"Z_ERROR({_f(zz_o)}) 0 1", f"{_p1(parts.idle)} 0 1"]
    L += [f"X_ERROR({_f(parts.hz)}) 0 1", "M 0 1"]
    return "\n".join(L) + "\n"


def line_program(parts: Parts) -> str:
    """Three patches merged at once: seams (0,1) and (1,2); records: two outcomes, then Z of each.
    The middle patch's errors are taken from the first seam's channel only (bit 2 of the second
    seam's pattern is dropped), so its idling is not counted twice."""
    zz_j, zz_o = parts.zz
    j2 = np.zeros(8)
    for k, v in enumerate(zz_j):
        j2[k & 0b101] += v  # keep outcome and the far patch
    L = ["R 0 1 2", f"X_ERROR({_f(parts.hz)}) 0 1 2"]
    L += ["R 3 4", "CX 0 3 1 3 1 4 2 4", "H 3 0 1", f"{_z3(zz_j)} 3 0 1", "H 3 0 1", "H 4 2", f"{_z3(j2)} 4 1 2", "H 4 2", "M 3 4"]
    L += [f"X_ERROR({_f(parts.hz)}) 0 1 2", "M 0 1 2"]
    return "\n".join(L) + "\n"


def predict(text: str, shots: int, seed: int) -> np.ndarray:
    return Program(text).sample(shots, seed=seed)


def compare(raw: dict, d: int, p: float, shots: int = 1 << 20, seed: int = 0) -> dict:
    comps = Components(raw)
    parts = Parts(comps, d, p)
    if not parts.ok():
        return {}
    out = {}

    def stat(flips: np.ndarray, rec: dict, labels: list[str]):
        meas = Counts.from_json(rec["counts"])
        res = {}
        for j, lab in enumerate(labels):
            pred = flips[:, j].mean()
            got = meas.marginal(j)
            sig = math.sqrt(got * (1 - got) / meas.shots + pred * (1 - pred) / flips.shape[0])
            res[lab] = {"measured": got, "predicted": float(pred), "ratio": float(pred / got) if got else None, "z": float((pred - got) / sig) if sig else None}
        anyp = flips.any(axis=1).mean()
        res["any"] = {"measured": meas.failures / meas.shots, "predicted": float(anyp)}
        return res

    for inp in ("z", "x"):
        rec = raw.get(f"cnot-{inp}-d{d}-p{p}")
        if rec:
            rows = predict(cnot_program(parts, inp), shots, seed)
            out[f"cnot-{inp}"] = stat(cnot_flips(rows, inp), rec, ["obs0", "obs1"])
    rec = raw.get(f"rep-k3-d{d}-p{p}")
    if rec:
        rows = predict(repeated_program(parts, 3), shots, seed + 1)
        out["repeated_zz_k3"] = stat(rows, rec, ["m1", "m2", "m3", "Z1", "Z2"])
    rec = raw.get(f"line-n3-d{d}-p{p}")
    if rec:
        rows = predict(line_program(parts), shots, seed + 2)
        out["line_n3"] = stat(rows, rec, ["seam01", "seam12", "Z1", "Z2", "Z3"])
    return out
