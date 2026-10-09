"""Evaluates physical qubit footprint savings of the XZZX surface code under biased noise.

Outputs: data/results/biased_xzzx.json
"""

from __future__ import annotations

import json
import math
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "python"))

from ftalgo.biased import evaluate_bias_scaling
from ftalgo.calib.model import LogicalModel

OUT = ROOT / "data" / "results" / "biased_xzzx.json"


def main():
    print("Loading logical calibration model...")
    model = LogicalModel.load()

    etas = [1.0, 5.0, 10.0, 25.0, 50.0, 100.0, 250.0, 500.0, 1000.0, 2000.0]

    workloads = [
        {
            "id": "rsa2048_modern",
            "name": "RSA-2048 (Modern Arithmetic)",
            "qubits": 8199,
            "rounds": 3.5e10,
            "target_faults": 0.01,
        },
        {
            "id": "femoco",
            "name": "FeMoco Nitrogenase (Lee et al. 2021)",
            "qubits": 2142,
            "rounds": 2.12e8,
            "target_faults": 0.01,
        },
        {
            "id": "lih_qpe",
            "name": "LiH Active Space QPE (10-bit)",
            "qubits": 3,
            "rounds": 2.5e7,
            "target_faults": 0.01,
        },
        {
            "id": "h2_qpe",
            "name": "H2 Ground State QPE (10-bit)",
            "qubits": 3,
            "rounds": 1.5e7,
            "target_faults": 0.01,
        },
    ]

    results = {}
    for w in workloads:
        print(f"Evaluating biased noise scaling for {w['name']}...")
        rows = evaluate_bias_scaling(
            total_space_time_rounds=w["rounds"],
            num_logical_qubits=w["qubits"],
            target_faults=w["target_faults"],
            p_total=1e-3,
            etas=etas,
            model=model,
        )
        results[w["id"]] = {
            "name": w["name"],
            "logical_qubits": w["qubits"],
            "space_time_rounds": w["rounds"],
            "p_total": 1e-3,
            "points": [
                {
                    "eta": r.eta,
                    "dx": r.dx,
                    "dz": r.dz,
                    "tile_qubits": r.tile_qubits,
                    "symmetric_d": r.symmetric_d,
                    "symmetric_tile_qubits": r.symmetric_tile_qubits,
                    "qubit_reduction_ratio": round(r.qubit_reduction_ratio, 3),
                    "total_physical_qubits": r.total_physical_qubits,
                    "symmetric_physical_qubits": r.symmetric_physical_qubits,
                    "target_error_achieved": r.target_error_achieved,
                }
                for r in rows
            ],
        }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(results, fh, indent=2)
    print(f"Saved biased XZZX scaling results to {OUT}")


if __name__ == "__main__":
    main()
