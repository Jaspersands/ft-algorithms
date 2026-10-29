"""From measured small instances to cryptographic and chemical scale.

A run's expected number of logical faults is

    E(d) = q_idle(d)·I + Σ_ops count_op · q_op(d)

with q the calibrated channel totals (``ftalgo.arch``) and I the idle qubit-rounds. For the
textbook Shor circuit the operation counts are exact functions of n (the generator's own
formula) and the schedule's shape is measured on compiled instances: rounds per Toffoli in
units of d, and the number of patches busy at a time (``fit_structure``). For published
algorithms the counts and the run time come from the paper (cited), and only the error model
is ours. The smallest odd d with E(d) ≤ target gives the physical qubits and the run time.

Published workloads are compared at their authors' own failure budgets, and the comparison is of
distance: their physical-qubit totals use layouts (yoked storage, hand-made floor plans) that this
model does not include.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .arch import FACTORIES, Architecture, Factory, ROUND_SECONDS
from .calib.model import LogicalModel


def textbook_toffolis(n: int) -> int:
    """Toffolis of ``shor.order_finding`` for an n-bit modulus: 2n controlled multiplications of
    2n modular additions (10n + 12 Toffolis each) and n controlled swaps."""
    return 2 * n * (2 * n * (10 * n + 12) + n)


@dataclass
class Structure:
    rounds_per_toffoli: float  # in units of d
    cx_per_toffoli: float
    busy: float                # patches busy at a time on average
    t_per_n: float             # synthesized-rotation T gates per bit of n
    s_per_n: float


def fit_structure(model: LogicalModel, Ns=(77, 143, 391, 899)) -> Structure:
    from . import shor
    from .schedule import compile_noisy

    rows = []
    for N in Ns:
        c, info = shor.order_finding(N)
        r = []
        for d in (11, 21):
            nz = compile_noisy(c, Architecture(model, d, 0.001, FACTORIES["cultivation"]))
            r.append((nz.rounds, nz.idle_qubit_rounds))
        R1 = (r[1][0] - r[0][0]) / 10
        cnt = c.counts()
        tof = c.toffoli_count()
        Q = c.num_qubits
        rows.append((R1 / tof, cnt["CX"] / tof, Q - r[1][1] / r[1][0], (cnt["T"] + cnt["T_DAG"]) / info.n, (cnt["S"] + cnt["S_DAG"]) / info.n))
    mean = [sum(x[i] for x in rows) / len(rows) for i in range(5)]
    return Structure(*mean)


def fit_modern_structure(model: LogicalModel, Ns=(15, 21, 35)) -> Structure:
    from . import arith_modern as am
    from .schedule import compile_noisy

    rows = []
    for N in Ns:
        c, info = am.order_finding_modern(N)
        r = []
        for d in (11, 21):
            nz = compile_noisy(c, Architecture(model, d, 0.001, FACTORIES["cultivation"]))
            r.append((nz.rounds, nz.idle_qubit_rounds))
        R1 = (r[1][0] - r[0][0]) / 10
        cnt = c.counts()
        tof = c.toffoli_count()
        Q = c.num_qubits
        rows.append((R1 / tof, cnt["CX"] / tof, Q - r[1][1] / r[1][0], (cnt["T"] + cnt["T_DAG"]) / info.n, (cnt["S"] + cnt["S_DAG"]) / info.n))
    mean = [sum(x[i] for x in rows) / len(rows) for i in range(5)]
    return Structure(*mean)


@dataclass
class Workload:
    name: str
    logical_qubits: int
    toffolis: float
    cnots: float
    t_gates: float
    s_gates: float
    rounds_per_d: float | None   # rounds = rounds_per_d · d (our schedule) …
    rounds_fixed: float | None   # … or a published run time in rounds
    source: str


def textbook_shor(n: int, st: Structure) -> Workload:
    tof = textbook_toffolis(n)
    return Workload(f"Shor, textbook arithmetic, n = {n}", 3 * n + 6, tof, st.cx_per_toffoli * tof, st.t_per_n * n, st.s_per_n * n,
                    st.rounds_per_toffoli * tof, None, "this work: exact counts from the generator, schedule shape fitted on compiled instances")


def modern_shor(n: int, st: Structure, k: int = 2) -> Workload:
    from . import arith_modern as am
    tof = am.modern_toffolis(n, k=k)
    return Workload(
        f"Shor, modern windowed arithmetic (k={k}), n = {n}",
        4 * n + 7,
        tof,
        st.cx_per_toffoli * tof,
        st.t_per_n * n,
        st.s_per_n * n,
        st.rounds_per_toffoli * tof,
        None,
        "this work: modern windowed arithmetic (k=2) with Gidney carry uncomputation",
    )


# Gidney 2025 budgets per shot: each shot is an independent attempt, and a shot lost to a logical
# error is simply repeated. Its expected 9.2 shots per factoring are not one long computation.
GIDNEY_2025_SHOTS = 9.2
GIDNEY_2025 = Workload(
    "RSA-2048, Gidney 2025 (one shot)", 1409, 6.5e9 / GIDNEY_2025_SHOTS, 0.0, 0.0, 0.0, None, 12.07 * 3600 / ROUND_SECONDS,
    "Gidney 2025 (arXiv:2505.15917), Section 3.2 and Table 5: one shot runs 12.07 h with at most 1,409 logical "
    "qubits active; 6.5e9 Toffolis per factoring over 9.2 expected shots",
)
# His budget: 1e-15 per logical-qubit round over 6.9e13 qubit-rounds, so 93.3% of shots see no
# logical error; reached at d = 25 by his projection of surface-code memory (his Figure 6).
GIDNEY_2025_PER_ROUND = 1e-15
GIDNEY_2025_TARGET = -math.log(0.933)
GIDNEY_2025_D = 25
FEMOCO_THC = Workload(
    "FeMoco (THC), Lee et al. 2021", 2142, 5.3e9, 0.0, 0.0, 0.0, None, 5.3e9 / 25e3 / ROUND_SECONDS,
    "Lee et al. 2021 (PRX Quantum 2, 030305), Table III: 2,142 logical qubits, 5.3e9 Toffolis; Toffolis at 25 kHz",
)
# Lee et al. choose d = 31 for a 1% total error budget (Section V).
FEMOCO_THC_TARGET = 0.01
FEMOCO_THC_D = 31

# Each published workload with its authors' failure budget (expected faults), distance, physical
# qubits and run time, as they report them.
PUBLISHED = (
    (GIDNEY_2025, GIDNEY_2025_TARGET, GIDNEY_2025_D, 897_864, "12.07 h per shot"),
    (FEMOCO_THC, FEMOCO_THC_TARGET, FEMOCO_THC_D, 4e6, "< 4 days"),
)


def idle_total(model, d: int, p: float) -> float:
    """Probability per round that an idle patch suffers any logical error (X, Y or Z)."""
    return sum(model.idle(d, p))


def distance_for_idle(model, target: float, p: float = 1e-3, dmax: int = 61) -> int | None:
    """Smallest odd d whose idle patch errs at most ``target`` per round: the criterion Gidney 2025
    uses to pick its distance, independent of how a workload is scheduled."""
    for d in range(3, dmax + 1, 2):
        if idle_total(model, d, p) <= target:
            return d
    return None


@dataclass
class Estimate:
    workload: str
    d: int | None
    expected_faults: float
    rounds: float
    seconds: float
    physical_qubits: float
    data_qubits: float
    factory_qubits: float
    extrapolated: bool
    breakdown: dict


def expected_faults(w: Workload, arch: Architecture, busy: float = 3.2) -> tuple[float, float, dict]:
    d = arch.d
    rounds = w.rounds_fixed if w.rounds_fixed is not None else w.rounds_per_d * d
    q_idle = sum(arch.idle(1))
    z3, qx = arch.ccz_channel()
    q_ccz = sum(z3) + 3 * qx + 3 * sum(arch.idle(d))
    br = {
        "idle": q_idle * max(w.logical_qubits - busy, 0) * rounds,
        "CCZ": w.toffolis * q_ccz,
        "CNOT": w.cnots * sum(arch.cnot),
        "T": w.t_gates * sum(arch.t_channel()),
        "S": w.s_gates * sum(arch.s_channel()),
    }
    return sum(br.values()), rounds, br


def estimate(w: Workload, model: LogicalModel, p: float = 1e-3, factory: str | Factory = "cultivation", target: float = 0.1, dmax: int = 61) -> Estimate:
    fac_model = FACTORIES[factory] if isinstance(factory, str) else factory
    for d in range(3, dmax + 1, 2):
        arch = Architecture(model, d, p, fac_model)
        ef, rounds, br = expected_faults(w, arch)
        if ef <= target:
            data = arch.data_qubits(w.logical_qubits)
            fac = arch.factory_qubits(w.t_gates, w.toffolis, rounds)
            return Estimate(w.name, d, ef, rounds, rounds * ROUND_SECONDS, data + fac, data, fac, bool(arch.extrapolated), br)
    return Estimate(w.name, None, float("inf"), 0, 0, 0, 0, 0, True, {})
