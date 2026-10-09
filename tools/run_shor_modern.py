"""Runs the modern arithmetic Shor grid (Gidney measurement-based uncomputation + windowing).

Saves results to data/results/shor_modern/.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import zlib

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))

import numpy as np

from ftalgo import arith_modern as am, shor
from ftalgo.calib.model import LogicalModel
from ftalgo.engine import Program
from ftalgo.experiment import d_range, run_config

OUT = os.path.join(os.path.dirname(__file__), "..", "data", "results", "shor_modern")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, nargs="*", default=[15, 21, 35])
    ap.add_argument("--p", type=float, nargs="*", default=[0.001])
    ap.add_argument("--factory", nargs="*", default=["cultivation", "15to1"])
    ap.add_argument("--shots", type=int, default=2000)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    model = LogicalModel.load()

    for N in args.N:
        circ, info = am.order_finding_modern(N)
        table = shor.success_table(N, info.a, info.m)
        peaks = shor.peak_table(N, info.a, info.m)

        def score(rows, info=info, table=table, peaks=peaks):
            y = info.outcomes(rows)
            return np.stack([peaks[y], table[y]], axis=1)  # [peak, factors found]

        base_path = os.path.join(OUT, f"N{N}-baselines.json")
        if not os.path.exists(base_path):
            print(f"Computing noiseless baseline for N={N}...", flush=True)
            prog = Program(circ.to_text())
            noiseless = score(prog.sample(10000, seed=1, faults="none")).mean(axis=0)
            json.dump({
                "N": N, "a": info.a, "n": info.n, "m": info.m, "cutoff": info.cutoff, "eps": info.eps,
                "logical_qubits": circ.num_qubits, "counts": info.counts, "toffoli": circ.toffoli_count(),
                "t": circ.t_count(),
                "ideal_success_exact_qft": shor.ideal_success(N, info.a, info.m),
                "scores": ["peak", "factors"],
                "ideal_peak_exact_qft": shor.ideal_peak_probability(N, info.a, info.m),
                "noiseless": [float(x) for x in noiseless],
                "random": [float(peaks.mean()), float(table.mean())],
                "order": shor.multiplicative_order(info.a, N),
            }, open(base_path, "w"), indent=1)

        for p in args.p:
            for fac in args.factory:
                ds = d_range(circ, model, p, fac)
                for d in ds:
                    name = f"N{N}-{fac}-p{p}-d{d}"
                    path = os.path.join(OUT, name + ".json")
                    if os.path.exists(path):
                        continue
                    r = run_config(
                        circ, score, model, d=d, p=p, factory=fac,
                        seed=zlib.crc32(name.encode()), shots=args.shots
                    )
                    r.update({"N": N, "a": info.a, "algorithm": "shor_modern"})
                    with open(path + ".tmp", "w") as f:
                        json.dump(r, f, indent=1)
                    os.replace(path + ".tmp", path)
                    print(
                        f"Modern {name}: peak {r['scores'][0]:.3f} +/- {r['score_sigmas'][0]:.3f} "
                        f"factors {r['scores'][1]:.3f}  E[faults] {r['expected_faults']:.3g}  "
                        f"{r['seconds']*1e3:.0f} ms  {r['physical_qubits']['total']:.3g} qubits  ({r['wall_seconds']} s)",
                        flush=True,
                    )


if __name__ == "__main__":
    main()
