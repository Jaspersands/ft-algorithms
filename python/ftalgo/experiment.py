"""One algorithm run on one machine configuration: compile with noise, estimate success, record
the error budget and the physical resources. Results are JSON-ready dicts."""

from __future__ import annotations

import time
from typing import Callable

import numpy as np

from .arch import FACTORIES, Architecture
from .calib.model import LogicalModel
from .circuit import Circuit
from .engine import Program
from .estimate import error_budget, stratified
from .schedule import compile_noisy


def run_config(
    circ: Circuit,
    score: Callable[[np.ndarray], np.ndarray],
    model: LogicalModel,
    *,
    d: int,
    p: float,
    factory: str,
    seed: int,
    shots: int = 4000,
    threads: int = 0,
) -> dict:
    t0 = time.time()
    arch = Architecture(model, d, p, FACTORIES[factory])
    noisy = compile_noisy(circ, arch)
    prog = Program(noisy.text)
    est = stratified(prog, score, seed=seed, shots=shots, threads=threads)
    vals = est.value if isinstance(est.value, list) else [est.value]
    sigs = est.sigma if isinstance(est.sigma, list) else [est.sigma]
    return {
        "d": d,
        "p": p,
        "factory": factory,
        "success": vals[0],
        "sigma": sigs[0],
        "scores": vals,
        "score_sigmas": sigs,
        "interval": list(est.interval()),
        "estimate": est.to_json(),
        "expected_faults": prog.expected_faults,
        "budget": dict(noisy.budget),
        "error_budget": error_budget(dict(noisy.budget), est),
        "rounds": noisy.rounds,
        "seconds": noisy.seconds,
        "logical_qubits": noisy.num_qubits,
        "physical_qubits": noisy.physical_qubits(),
        "t_count": noisy.t_count,
        "ccz_count": noisy.ccz_count,
        "eps_t": arch.eps_t,
        "eps_ccz": arch.eps_ccz,
        "extrapolated": bool(arch.extrapolated),
        "wall_seconds": round(time.time() - t0, 1),
    }


def d_range(circ: Circuit, model: LogicalModel, p: float, factory: str, lo_faults: float = 40.0, hi_faults: float = 0.003, dmax: int = 31) -> list[int]:
    """Odd distances from where the run expects ≲ lo_faults faults to where it expects ≲ hi_faults."""
    out = []
    for d in range(3, dmax + 1, 2):
        ef = Program(compile_noisy(circ, Architecture(model, d, p, FACTORIES[factory])).text).expected_faults
        if ef <= lo_faults:
            out.append(d)
        if ef <= hi_faults:
            break
    return out
