"""Executes the sensitivity analysis campaign for the Gidney 2025 headline (d = 29).

Generates data/results/sensitivity.json.
"""

from __future__ import annotations

import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))

from ftalgo.calib.model import LogicalModel
from ftalgo.sensitivity import run_sensitivity_tornado

OUT = os.path.join(os.path.dirname(__file__), "..", "data", "results", "sensitivity.json")


def main():
    model = LogicalModel.load()
    print("Running sensitivity analysis campaign for RSA-2048 headline...")
    results = run_sensitivity_tornado(model)

    print("\n" + "=" * 80)
    print("SENSITIVITY ANALYSIS TORNADO TABLE (Gidney 2025 RSA-2048)")
    print("=" * 80)
    b = results["baseline"]
    print(f"BASELINE: d = {b['d']}, Physical Qubits = {b['physical_qubits']:.3g}, E[faults] = {b['expected_faults']:.3g}\n")

    print(f"{'Category':<16} {'Perturbation / Scenario':<42} {'d':<5} {'Delta d':<8} {'Physical Qubits':<16}")
    print("-" * 80)
    for exp in results["experiments"]:
        d_str = str(exp["d"]) if exp["d"] is not None else "None"
        delta_str = f"{exp['delta_d']:+d}" if exp["delta_d"] is not None else "N/A"
        pq_str = f"{exp['physical_qubits']:.3g}"
        print(f"{exp['category']:<16} {exp['name']:<42} {d_str:<5} {delta_str:<8} {pq_str:<16}")

    cov = results["fit_covariance"]
    print("\n" + "=" * 80)
    print("EXTRAPOLATION UNCERTAINTY BAND (idle_xy per-round logical error rate at p=0.001)")
    print("=" * 80)
    print(f"Slope: Lambda = {cov['lambda']:.2f} (1-sigma range: [{cov['lambda_lo']:.2f}, {cov['lambda_hi']:.2f}])")
    print(f"{'d':<5} {'Central rate':<15} {'1-sigma range':<32} {'2-sigma range':<32}")
    print("-" * 80)
    for d in [9, 11, 15, 21, 25, 27, 29, 31, 35]:
        entry = cov["band"][d]
        c_str = f"{entry['central']:.3e}"
        s1_str = f"[{entry['lo_1sigma']:.3e}, {entry['hi_1sigma']:.3e}]"
        s2_str = f"[{entry['lo_2sigma']:.3e}, {entry['hi_2sigma']:.3e}]"
        print(f"{d:<5} {c_str:<15} {s1_str:<32} {s2_str:<32}")

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(results, f, indent=1)
    print(f"\nSaved results to {OUT}")


if __name__ == "__main__":
    main()
