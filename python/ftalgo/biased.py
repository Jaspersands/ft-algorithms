"""XZZX surface codes and asymmetric rectangular patches under biased noise.

In platforms such as dual-rail superconducting qubits, Kerr-cat qubits, and acoustic-wave resonators,
dephasing noise dwarfs bit-flip noise:
    η = p_Z / p_X >> 1  (often η in [10, 1000]).

In the XZZX surface code (Bonilla Ataides et al., Nature Comms 12, 3812, 2021), stabilizer generators
X-Z-Z-X decouple pure dephasing noise into independent 1D repetition codes along diagonals, raising
the threshold to 50% under pure dephasing.

Under finite bias η:
    p_X = p_Y = p / (2(1 + η))
    p_Z = η * p / (1 + η)

Because bit flips occur with probability p_X << p, the code distance d_Z needed to protect against
bit flips is substantially smaller than the distance d_X needed to protect against phase flips.
An asymmetric rectangular patch of dimensions d_X × d_Z requires:
    Physical qubits per tile: 2(d_X + 1)(d_Z + 1)
compared to a standard symmetric patch:
    2(d_sym + 1)²

This module models:
1. Physical error rates under dephasing bias η.
2. Analytic & calibrated logical error rates P_L(d_X, d_Z, p, η).
3. Optimal distance selection (d_X, d_Z) under a target logical error budget.
4. Physical qubit savings and comparison across cryptographic and chemistry workloads.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

import numpy as np

from .calib.model import LogicalModel


@dataclass(frozen=True)
class BiasedNoise:
    p_total: float
    eta: float  # bias = p_Z / (2 * p_X) when p_X == p_Y

    @property
    def px(self) -> float:
        """Bit-flip probability p_X."""
        return self.p_total / (2.0 * (1.0 + self.eta))

    @property
    def py(self) -> float:
        """Y error probability p_Y = p_X."""
        return self.p_total / (2.0 * (1.0 + self.eta))

    @property
    def pz(self) -> float:
        """Phase-flip probability p_Z."""
        return self.eta * self.p_total / (1.0 + self.eta)


@dataclass(frozen=True)
class XZZXPatch:
    dx: int  # distance against Z errors (along X boundary)
    dz: int  # distance against X errors (along Z boundary)

    def __post_init__(self):
        if self.dx < 3 or self.dx % 2 == 0:
            raise ValueError(f"dx must be an odd integer >= 3, got {self.dx}")
        if self.dz < 3 or self.dz % 2 == 0:
            raise ValueError(f"dz must be an odd integer >= 3, got {self.dz}")

    @property
    def tile_qubits(self) -> int:
        """Physical qubits of one patch including its share of routing: 2(d_X + 1)(d_Z + 1)."""
        return 2 * (self.dx + 1) * (self.dz + 1)

    @property
    def data_qubits(self) -> int:
        """Data and ancilla qubits on the patch without routing: (d_X + 1)(d_Z + 1)."""
        return (self.dx + 1) * (self.dz + 1)


@dataclass(frozen=True)
class LogicalBiasedRates:
    pz_logical: float
    px_logical: float
    ptot_logical: float


def logical_rates_analytic(
    patch: XZZXPatch,
    noise: BiasedNoise,
    p_star: float = 0.008525875,
    a_prefactor: float = 0.06165385,
) -> LogicalBiasedRates:
    """Computes logical error rates per syndrome round for an asymmetric XZZX patch.

    P_{L,Z} = A * (p_Z / p_star)^((d_X + 1) / 2)
    P_{L,X} = A * (p_X / p_star)^((d_Z + 1) / 2)
    P_total = P_{L,Z} + P_{L,X}
    """
    pz_log = a_prefactor * ((noise.pz / p_star) ** ((patch.dx + 1) / 2))
    px_log = a_prefactor * ((noise.px / p_star) ** ((patch.dz + 1) / 2))
    return LogicalBiasedRates(
        pz_logical=pz_log,
        px_logical=px_log,
        ptot_logical=pz_log + px_log,
    )


def logical_rates_calibrated(
    patch: XZZXPatch,
    noise: BiasedNoise,
    model: LogicalModel,
) -> LogicalBiasedRates:
    """Uses the fitted calibration model parameters to evaluate logical rates."""
    gfit = model.glob.get("idle_xy")
    if gfit is not None:
        log_a = gfit.logA
        log_pstar = gfit.logpstar
        p_star = math.exp(log_pstar)
        a_prefactor = math.exp(log_a)
    else:
        p_star = 0.008525875
        a_prefactor = 0.06165385

    return logical_rates_analytic(patch, noise, p_star=p_star, a_prefactor=a_prefactor)


def optimize_patch(
    target_p_per_round: float,
    p_total: float = 1e-3,
    eta: float = 100.0,
    model: LogicalModel | None = None,
    max_d: int = 61,
) -> tuple[XZZXPatch, LogicalBiasedRates]:
    """Finds the odd pair (d_X, d_Z) minimizing patch footprint while achieving P_L <= target."""
    noise = BiasedNoise(p_total, eta)
    if model is not None and "idle_xy" in model.glob:
        p_star = math.exp(model.glob["idle_xy"].logpstar)
        a_prefactor = math.exp(model.glob["idle_xy"].logA)
    else:
        p_star = 0.008525875
        a_prefactor = 0.06165385

    # Target allocation: half budget to Z, half to X
    target_each = target_p_per_round / 2.0

    # Solve for dx: a * (pz / p_star)^((dx + 1) / 2) <= target_each
    ratio_z = noise.pz / p_star
    if ratio_z >= 1.0:
        raise ValueError(f"p_Z={noise.pz} is above or equal to threshold p*={p_star}")
    req_hx = math.log(target_each / a_prefactor) / math.log(ratio_z)
    req_dx = 2.0 * req_hx - 1.0
    dx = max(3, int(math.ceil(req_dx)))
    if dx % 2 == 0:
        dx += 1

    # Solve for dz: a * (px / p_star)^((dz + 1) / 2) <= target_each
    ratio_x = noise.px / p_star
    req_hz = math.log(target_each / a_prefactor) / math.log(ratio_x)
    req_dz = 2.0 * req_hz - 1.0
    dz = max(3, int(math.ceil(req_dz)))
    if dz % 2 == 0:
        dz += 1

    dx = min(dx, max_d)
    dz = min(dz, max_d)
    patch = XZZXPatch(dx, dz)
    rates = logical_rates_analytic(patch, noise, p_star=p_star, a_prefactor=a_prefactor)
    return patch, rates


@dataclass
class BiasedComparisonRow:
    eta: float
    dx: int
    dz: int
    tile_qubits: int
    symmetric_d: int
    symmetric_tile_qubits: int
    qubit_reduction_ratio: float
    total_physical_qubits: int
    symmetric_physical_qubits: int
    target_error_achieved: float


def evaluate_bias_scaling(
    total_space_time_rounds: float,
    num_logical_qubits: int,
    target_faults: float = 0.01,
    p_total: float = 1e-3,
    etas: Sequence[float] = (1.0, 10.0, 50.0, 100.0, 300.0, 500.0, 1000.0),
    model: LogicalModel | None = None,
) -> list[BiasedComparisonRow]:
    """Evaluates physical qubit savings across bias levels for a given space-time volume.

    A symmetric square patch must choose distance d_sym = max(d_X, d_Z) = d_X in order to
    suppress the dominant phase errors, forcing an isotropic footprint of 2(d_X + 1)².
    An asymmetric rectangular patch tailors d_Z to the suppressed bit-flip rate, reducing
    the footprint to 2(d_X + 1)(d_Z + 1).
    """
    target_p_per_round = target_faults / (num_logical_qubits * total_space_time_rounds)

    # Litinski fast-block formula: tiles = 2n + ceil(sqrt(8n)) + 1
    tiles = 2 * num_logical_qubits + math.ceil(math.sqrt(8 * num_logical_qubits)) + 1

    rows = []
    for eta in etas:
        patch, rates = optimize_patch(target_p_per_round, p_total, eta=eta, model=model)
        # Symmetric square patch required for this noise environment must have d = max(dx, dz)
        sym_d = max(patch.dx, patch.dz)
        sym_tile = 2 * (sym_d + 1) ** 2
        sym_total_phys = tiles * sym_tile

        rect_tile = patch.tile_qubits
        rect_total_phys = tiles * rect_tile
        reduction = sym_tile / rect_tile
        rows.append(
            BiasedComparisonRow(
                eta=float(eta),
                dx=patch.dx,
                dz=patch.dz,
                tile_qubits=rect_tile,
                symmetric_d=sym_d,
                symmetric_tile_qubits=sym_tile,
                qubit_reduction_ratio=reduction,
                total_physical_qubits=rect_total_phys,
                symmetric_physical_qubits=sym_total_phys,
                target_error_achieved=rates.ptot_logical * num_logical_qubits * total_space_time_rounds,
            )
        )
    return rows
