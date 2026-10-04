"""Phase estimation of molecular ground-state energies on the simulated machine: probability of
chemical accuracy (|E − E_FCI| < 1.6 mHa) vs d for each molecule and factory.

    python tools/run_qpe.py --mol H2 HeH+ --factory cultivation 15to1

One JSON per configuration in data/results/qpe/, plus a baseline file per molecule.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import zlib

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))

import numpy as np  # noqa: E402

from ftalgo import qpe  # noqa: E402
from ftalgo.calib.model import LogicalModel  # noqa: E402
from ftalgo.chem.qubit import molecule, qubit_problem  # noqa: E402
from ftalgo.engine import Program  # noqa: E402
from ftalgo.experiment import d_range, run_config  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), "..", "data", "results", "qpe")
PARAMS = {"m": 10, "tau": 2 * math.pi, "steps": 4, "order": 4}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mol", nargs="*", default=["H2", "HeH+"])
    ap.add_argument("--p", type=float, nargs="*", default=[0.001])
    ap.add_argument("--factory", nargs="*", default=["cultivation", "15to1"])
    ap.add_argument("--shots", type=int, default=1000)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    model = LogicalModel.load()
    for name in args.mol:
        prob = qubit_problem(molecule(name))
        circ, info = qpe.qpe_circuit(prob.hamiltonian, prob.hf_bits, e_ref=prob.e_hf, **PARAMS)

        def score(rows, info=info, prob=prob):
            e = info.energies(rows)
            return np.stack([np.abs(e - prob.e_fci) < qpe.CHEMICAL_ACCURACY, np.abs(e - prob.e_fci) < 4 * qpe.CHEMICAL_ACCURACY], axis=1)

        tag = name.replace("+", "p")
        base = os.path.join(OUT, f"{tag}-baselines.json")
        if not os.path.exists(base):
            ideal = qpe.ideal_distribution(prob.hamiltonian, prob.hf_bits, info)
            y = np.arange(1 << info.m)
            phi = np.where(y / (1 << info.m) >= 0.5, y / (1 << info.m) - 1, y / (1 << info.m))
            E = info.e_ref - 2 * math.pi * phi / info.tau
            noiseless = score(Program(circ.to_text()).sample(300, seed=1, faults="none")).mean(axis=0)
            json.dump({
                "molecule": name, "qubits": prob.hamiltonian.n, "terms": len(prob.hamiltonian),
                "e_hf": prob.e_hf, "e_fci": prob.e_fci, "params": PARAMS, "eps": info.eps,
                "rotations_per_run": info.rotations_per_run, "t_count": circ.t_count(), "counts": dict(circ.counts()),
                "resolution": info.resolution(),
                "ideal_chemical_accuracy": float(ideal[np.abs(E - prob.e_fci) < qpe.CHEMICAL_ACCURACY].sum()),
                "scores": ["chemical accuracy (1.6 mHa)", "within 6.4 mHa"],
                "noiseless": [float(x) for x in noiseless],
                "random": [float(np.mean(np.abs(E - prob.e_fci) < qpe.CHEMICAL_ACCURACY)), float(np.mean(np.abs(E - prob.e_fci) < 4 * qpe.CHEMICAL_ACCURACY))],
            }, open(base, "w"), indent=1)
        for p in args.p:
            for fac in args.factory:
                for d in d_range(circ, model, p, fac, lo_faults=12.0, hi_faults=0.01, dmax=35):
                    cfg = f"{tag}-{fac}-p{p}-d{d}"
                    path = os.path.join(OUT, cfg + ".json")
                    if os.path.exists(path):
                        continue
                    r = run_config(circ, score, model, d=d, p=p, factory=fac, seed=zlib.crc32(cfg.encode()), shots=args.shots)
                    r.update({"molecule": name, "algorithm": "qpe"})
                    with open(path + ".tmp", "w") as f:
                        json.dump(r, f, indent=1)
                    os.replace(path + ".tmp", path)
                    print(f"{cfg}: chem.acc. {r['scores'][0]:.3f} ± {r['score_sigmas'][0]:.3f}  E[faults] {r['expected_faults']:.3g}  "
                          f"{r['seconds']:.1f} s  {r['physical_qubits']['total']:.3g} qubits  ({r['wall_seconds']} s)", flush=True)


if __name__ == "__main__":
    main()
