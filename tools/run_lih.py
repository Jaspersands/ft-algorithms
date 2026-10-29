"""Simulates Lithium Hydride (LiH) ab initio chemistry:
1. Potential energy dissociation curve (RHF, CASCI(2,2) on RHF orbitals, CASSCF(2,2), full FCI) across
   R in [1.0, 3.0] Angstrom.
2. 10-bit QPE ground-state energy estimation on the simulated surface-code machine across distances.

Output: data/results/lih.json
"""

from __future__ import annotations

import argparse
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
SHOTS = 400  # each LiH shot is ~0.5 s (4.2e7 T gates); stratified estimator, min 60 per stratum
MIN_SHOTS = 60


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default=None, help="a fitted model other than data/calibration/model.json")
    ap.add_argument("--out", default=str(OUT))
    args = ap.parse_args()
    out = pathlib.Path(args.out)
    print("Computing LiH potential energy dissociation curve...")
    r_sweep = np.linspace(1.0, 3.0, 15)
    curve_pts = lih_dissociation_curve(r_sweep)
    curve_data = [
        {
            "r_angstrom": round(p.r_angstrom, 4),
            "r_bohr": round(p.r_bohr, 4),
            "e_hf": p.e_hf,
            "e_casci": p.e_casci,
            "e_casscf": p.e_casscf,
            "e_fci": p.e_fci,
        }
        for p in curve_pts
    ]

    print("Setting up LiH QPE at equilibrium (R = 1.595 A)...")
    prob = lih_problem(1.595, ncas=2, orbitals="casscf")
    from pyscf import fci, gto, scf

    _mol = gto.M(atom="Li 0 0 0; H 0 0 1.595", basis="sto-3g", verbose=0)
    e_fci_full = float(fci.FCI(scf.RHF(_mol).run()).kernel()[0])
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
                np.abs(e - e_fci_full) < qpe.CHEMICAL_ACCURACY,
            ],
            axis=1,
        )

    # Noiseless sampling baseline
    p_noiseless = Program(c.to_text())
    noiseless_shots = p_noiseless.sample(300, seed=1, faults="none")
    noiseless_score = score(noiseless_shots).mean(axis=0)

    model = LogicalModel.load(args.model) if args.model else LogicalModel.load()
    runs = []
    if out.exists():
        try:
            prev = json.loads(out.read_text())
            # Only resume runs computed for the same Hamiltonian and shot budget.
            if prev.get("orbitals") == "casscf" and prev.get("shots") == SHOTS:
                runs = prev.get("runs", [])
            print(f"Loaded {len(runs)} existing runs from {out}", flush=True)
        except Exception:
            runs = []

    def save_current():
        result = {
            "molecule": "LiH",
            "r_equilibrium": 1.595,
            "ncas": 2,
            "qubits": prob.hamiltonian.n,
            "terms": len(prob.hamiltonian),
            "e_hf": prob.e_hf,
            "orbitals": "casscf",
            "e_active": prob.e_fci,
            "e_fci_full": e_fci_full,
            "scores_meaning": ["|E - E_active| < 1.6 mHa", "|E - E_active| < 6.4 mHa", "|E - E_FCI(full STO-3G)| < 1.6 mHa"],
            "shots": SHOTS,
            "min_shots": MIN_SHOTS,
            "t_count": c.t_count(),
            "rotations_per_run": info.rotations_per_run,
            "noiseless": [float(s) for s in noiseless_score],
            "dissociation_curve": curve_data,
            "runs": runs,
        }
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(result, indent=2))

    save_current()

    configs = []
    for d in [13, 15, 17, 19, 21, 23]:
        configs.append(("cultivation", d))
    for d in [15, 17, 19, 21]:
        configs.append(("15to1", d))

    done_keys = {(r["factory"], r["d"]) for r in runs}

    for fac, d in configs:
        if (fac, d) in done_keys:
            print(f"Skipping already computed LiH-{fac}-p0.001-d{d}", flush=True)
            continue
        cfg_tag = f"LiH-{fac}-p0.001-d{d}"
        seed = zlib.crc32(cfg_tag.encode())
        print(f"Running QPE simulation: {cfg_tag}...", flush=True)
        res = run_config(
            c,
            score,
            model,
            d=d,
            p=0.001,
            factory=fac,
            seed=seed,
            shots=SHOTS,
            min_shots=MIN_SHOTS,
        )
        run_record = {
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
            "wall_seconds": res.get("wall_seconds", 0),
        }
        runs.append(run_record)
        done_keys.add((fac, d))
        save_current()
        print(
            f"{cfg_tag}: chem.acc. {res['scores'][0]:.3f} ± {res['score_sigmas'][0]:.3f}  "
            f"E[faults] {res['expected_faults']:.3g}  ({res['wall_seconds']} s) [saved to {out}]",
            flush=True,
        )

    print(f"Completed all LiH runs! Total runs: {len(runs)}", flush=True)


if __name__ == "__main__":
    main()
