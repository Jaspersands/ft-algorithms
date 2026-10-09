"""Simulates Lithium Hydride (LiH) ab initio chemistry:
1. Potential energy dissociation curve (HF, CASCI, full FCI) across R in [1.0, 3.0] Angstrom.
2. 10-bit QPE ground-state energy estimation on the simulated surface-code machine across distances.

Output: data/results/lih.json
"""

from __future__ import annotations

import json
import math
import pathlib
import sys
import zlib

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "python"))

import numpy as np

from ftalgo import qpe
from ftalgo.calib.model import LogicalModel
from ftalgo.chem.lih import lih_dissociation_curve, lih_problem
from ftalgo.engine import Program
from ftalgo.experiment import run_config

OUT = ROOT / "data" / "results" / "lih.json"


def main():
    print("Computing LiH potential energy dissociation curve...")
    r_sweep = np.linspace(1.0, 3.0, 15)
    curve_pts = lih_dissociation_curve(r_sweep)
    curve_data = [
        {
            "r_angstrom": round(p.r_angstrom, 4),
            "r_bohr": round(p.r_bohr, 4),
            "e_hf": p.e_hf,
            "e_casci": p.e_casci,
            "e_fci": p.e_fci,
        }
        for p in curve_pts
    ]

    print("Setting up LiH QPE at equilibrium (R = 1.595 A)...")
    prob = lih_problem(1.595, ncas=2)
    c, info = qpe.qpe_circuit(
        prob.hamiltonian,
        prob.hf_bits,
        m=10,
        tau=2 * math.pi,
        steps=4,
        e_ref=prob.e_hf,
        order=4,
    )

    def score(rows):
        e = info.energies(rows)
        return np.stack(
            [
                np.abs(e - prob.e_fci) < qpe.CHEMICAL_ACCURACY,
                np.abs(e - prob.e_fci) < 4 * qpe.CHEMICAL_ACCURACY,
            ],
            axis=1,
        )

    # Noiseless sampling baseline
    p_noiseless = Program(c.to_text())
    noiseless_shots = p_noiseless.sample(500, seed=1, faults="none")
    noiseless_score = score(noiseless_shots).mean(axis=0)

    model = LogicalModel.load()
    runs = []
    distances = [13, 15, 17, 19, 21, 23]
    for fac in ("cultivation", "15to1"):
        for d in distances:
            cfg_tag = f"LiH-{fac}-p0.001-d{d}"
            seed = zlib.crc32(cfg_tag.encode())
            print(f"Running QPE simulation: {cfg_tag}...")
            res = run_config(
                c,
                score,
                model,
                d=d,
                p=0.001,
                factory=fac,
                seed=seed,
                shots=500,
            )
            runs.append(
                {
                    "d": d,
                    "p": 0.001,
                    "factory": fac,
                    "scores": [float(s) for s in res["scores"]],
                    "score_sigmas": [float(s) for s in res["score_sigmas"]],
                    "expected_faults": float(res["expected_faults"]),
                    "rounds": int(res["rounds"]),
                    "seconds": float(res["seconds"]),
                    "physical_qubits": res["physical_qubits"],
                    "extrapolated": bool(res.get("extrapolated", d > 9)),
                }
            )

    result = {
        "molecule": "LiH",
        "r_equilibrium": 1.595,
        "ncas": 2,
        "qubits": prob.hamiltonian.n,
        "terms": len(prob.hamiltonian),
        "e_hf": prob.e_hf,
        "e_casci": prob.e_fci,
        "t_count": c.t_count(),
        "rotations_per_run": info.rotations_per_run,
        "noiseless": [float(s) for s in noiseless_score],
        "dissociation_curve": curve_data,
        "runs": runs,
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=2))
    print(f"Saved LiH results to {OUT}")


if __name__ == "__main__":
    main()
