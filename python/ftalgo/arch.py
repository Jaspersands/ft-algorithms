"""The simulated machine: a rotated surface code at distance d under SD6 noise of strength p,
operated by lattice surgery, with magic states from a chosen factory.

Every logical channel comes from the calibration (``ftalgo.calib.model.LogicalModel``) except the
magic states' own error rates and footprints, which are cited inputs (marked † where shown).

Operation model (rounds of syndrome extraction; 1 round = 1 µs):

- X, Y, Z: Pauli frame, free.
- H: transversal, with the patch's boundary orientation tracked (Litinski's fast block gives
  access to both boundaries), free.
- S, S_DAG: one Z⊗Z-type merge with a |Y⟩ resource, d rounds: the merge channel on the patch; a
  wrong outcome leaves a wrong Z correction.
- T, T_DAG: one Z⊗Z merge with a |T⟩ patch, d rounds (auto-corrected π/8 rotation): the merge
  channel, Z with the magic state's error ε_T †, and a wrong outcome leaves a wrong S correction,
  which Pauli twirling turns into Z with probability ½.
- CX, CZ, CY: the surgery CNOT, 2d rounds, with its measured channel (idling included).
- SWAP: three CNOTs.
- CCX, CCZ: consumption of a |CCZ⟩ state: three Z⊗Z merges in parallel (d rounds) then a
  Clifford correction step (d rounds of idling); Z-type errors with the state's error ε_CCZ †
  spread evenly over the seven non-trivial patterns; a wrong outcome on one merge leaves a wrong
  CZ correction on the other two qubits, twirled into Z₁, Z₂ or Z₁Z₂ with probability ¼ each.
- M, MX, R, RX: transversal, one round of idling.
- Idling: the measured per-round (p_X, p_Y, p_Z) for every round a live patch waits.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable

import numpy as np

from .calib import channels as ch
from .calib.model import LogicalModel

ROUND_SECONDS = 1e-6
REACTION_ROUNDS = 10


@dataclass(frozen=True)
class Factory:
    name: str
    label: str
    eps_t: Callable[[float], float]
    eps_ccz: Callable[[float], float]
    t_qubit_rounds: Callable[[int], float]
    ccz_qubit_rounds: Callable[[int], float]
    source: str


def _tile(d: int) -> int:
    """Physical qubits of one patch with its share of routing space: 2(d + 1)²."""
    return 2 * (d + 1) ** 2


FACTORIES = {
    "injected": Factory(
        "injected",
        "magic-state injection",
        lambda p: p,
        lambda p: 28 * p**2,
        lambda d: _tile(d) * d,
        lambda d: 8 * _tile(d) * d + 15 * _tile(d) * 2 * d,
        "ε_T = p (Li 2015 finds ≈ 0.4p for injection under this noise; we round up); "
        "|CCZ⟩ by 8T-to-CCZ distillation, ε = 28 ε_T² (Jones 2013; Gidney & Fowler 2019).",
    ),
    "15to1": Factory(
        "15to1",
        "15-to-1 distillation",
        lambda p: 4.5e-8 * (p / 1e-3) ** 3,
        lambda p: 5.2e-11 * (p / 1e-3) ** 6,
        lambda d: 4620 * 42.6,
        lambda d: 47000 * 60.0,
        "Litinski 2019 (Quantum 3, 205), Table 1 at p = 1e-3: (15-to-1)17,7,7, 4.5e-8, 4,620 qubits, "
        "42.6 cycles; (15-to-1)⁶ × (8-to-CCZ), 5.2e-11 per CCZ, 47,000 qubits, 60 cycles. Other p "
        "scaled as p³ (one level) and p⁶ (two levels).",
    ),
    "cultivation": Factory(
        "cultivation",
        "magic-state cultivation",
        lambda p: 2e-9 * (p / 1e-3) ** 5.64,
        lambda p: 28 * (2e-9 * (p / 1e-3) ** 5.64) ** 2,
        lambda d: 2.4e6 / 8,
        lambda d: 2.4e6,
        "Gidney, Shutty & Jones 2024: ε_T = 2e-9 at p = 1e-3, 4e-11 at p = 5e-4 (power law between, "
        "exponent 5.64); |CCZ⟩ by 8T-to-CCZ, 28 ε_T² (Gidney 2025); footprint from Gidney 2025's "
        "factories: six of 3×4 hot patches (1,352 qubits each) give one CCZ per 25 rounds, "
        "2.4e6 qubit·rounds per CCZ.",
    ),
}


def _one_minus_pow(a: float, n: float) -> float:
    """1 − (1 − 2a)^n without cancellation when a is tiny."""
    return -math.expm1(n * math.log1p(-2 * a))


def compose_p1(q: tuple[float, float, float], n: int | float) -> tuple[float, float, float]:
    """A single-qubit Pauli channel applied n times. Through the Pauli fidelities
    f_X = 1 − 2(p_Y + p_Z) etc., written with g = 1 − f^n so that probabilities far below the
    double-precision epsilon survive (per-round rates reach 1e-20 at large d)."""
    px, py, pz = q
    ga, gb, gc = _one_minus_pow(py + pz, n), _one_minus_pow(px + pz, n), _one_minus_pow(px + py, n)
    return (max((-ga + gb + gc) / 4, 0.0), max((ga - gb + gc) / 4, 0.0), max((ga + gb - gc) / 4, 0.0))


def xor_flips(*qs: float) -> float:
    """Probability that an odd number of independent flips (probabilities qs) occur."""
    return -math.expm1(sum(math.log1p(-2 * q) for q in qs)) / 2


@dataclass
class Architecture:
    model: LogicalModel
    d: int
    p: float
    factory: Factory
    reaction: int = REACTION_ROUNDS
    _cache: dict = field(default_factory=dict, repr=False)

    def __post_init__(self):
        m, d, p = self.model, self.d, self.p
        self.idle_round = m.idle(d, p)
        self.cnot = [max(v, 0.0) for v in m.cnot(d, p)]
        zz = m.pauli_measurement("zz", d, p)
        j = np.clip(zz["joint"], 0, None)
        j = j / j.sum()
        self.merge_m = float(sum(j[k] for k in range(8) if k & 1))   # wrong outcome
        self.merge_x = float(sum(j[k] for k in range(8) if k & 2))   # X on the data patch
        self.merge_z = max(float(zz["o"]), 0.0)                      # Z on the data patch
        self.eps_t = self.factory.eps_t(p)
        self.eps_ccz = self.factory.eps_ccz(p)
        self.extrapolated = m.extrapolated(d, p)

    def idle(self, rounds: int) -> tuple[float, float, float]:
        key = ("idle", rounds)
        if key not in self._cache:
            self._cache[key] = compose_p1(self.idle_round, rounds)
        return self._cache[key]

    def single(self, z_extra: float) -> tuple[float, float, float]:
        """The data patch's channel in a merge, with extra Z flips z_extra."""
        qx = self.merge_x
        qz = xor_flips(self.merge_z, z_extra)
        return (qx * (1 - qz), qx * qz, qz * (1 - qx))

    def t_channel(self):
        return self.single(xor_flips(self.eps_t, 0.5 * self.merge_m))

    def s_channel(self):
        return self.single(self.merge_m)

    def ccz_channel(self) -> tuple[list[float], tuple]:
        """(Z_CHANNEL_3 probabilities, per-qubit X probability)."""
        key = "ccz"
        if key in self._cache:
            return self._cache[key]
        f = np.ones(8)
        # ε_CCZ spread over the 7 non-trivial Z patterns.
        dist = np.zeros(8)
        dist[0] = 1 - self.eps_ccz
        dist[1:] = self.eps_ccz / 7
        f *= ch.fidelities(dist)
        for i in range(3):
            single = np.zeros(8)
            single[0] = 1 - self.merge_z
            single[1 << i] = self.merge_z
            f *= ch.fidelities(single)
            # A wrong outcome on merge i: twirled CZ on the other two.
            j, k = [q for q in range(3) if q != i]
            twirl = np.zeros(8)
            twirl[0] = 1 - 0.75 * self.merge_m
            for pat in (1 << j, 1 << k, (1 << j) | (1 << k)):
                twirl[pat] += 0.25 * self.merge_m
            f *= ch.fidelities(twirl)
        probs = ch.probabilities(f)
        self._cache[key] = (list(probs[1:]), self.merge_x)
        return self._cache[key]

    # -- physical resources ----------------------------------------------------------------------
    def data_qubits(self, n_logical: int) -> int:
        """Litinski's fast block: 2n + ⌈√(8n)⌉ + 1 tiles of 2(d + 1)² qubits."""
        tiles = 2 * n_logical + math.ceil(math.sqrt(8 * n_logical)) + 1
        return tiles * _tile(self.d)

    def factory_qubits(self, t_count: float, ccz_count: float, rounds: float) -> float:
        """Average factory footprint needed to supply the run's magic states in its duration."""
        if rounds <= 0:
            return 0.0
        return (t_count * self.factory.t_qubit_rounds(self.d) + ccz_count * self.factory.ccz_qubit_rounds(self.d)) / rounds
