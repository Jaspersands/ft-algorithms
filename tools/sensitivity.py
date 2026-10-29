"""Checks the measured model's code distance for Gidney 2025 (RSA-2048) against his d = 25: the
per-round criterion he uses, one shot of his algorithm at his budget, and how that moves with the
fit, the extrapolation, the magic states and the budget (``ftalgo.sensitivity``).

Writes data/results/sensitivity.json.
"""

from __future__ import annotations

import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))

from ftalgo.calib.model import LogicalModel  # noqa: E402
from ftalgo.sensitivity import run  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), "..", "data", "results", "sensitivity.json")


def main():
    from ftalgo.calib.model import DATA

    mr_path = os.path.join(DATA, "model_mr.json")
    r = run(LogicalModel.load(), model_mr=LogicalModel.load(mr_path) if os.path.exists(mr_path) else None)
    pr = r["per_round"]
    print(f"per-round criterion ({pr['target']:.0e}): measured d = {pr['d_measured']} "
          f"(1σ {pr['d_band']['-1']}–{pr['d_band']['1']}), Gidney's projection d = {pr['d_gidney']}")
    for row in pr["rows"]:
        if row["d"] in (3, 7, 11, 21, 25, 29):
            print(f"  d = {row['d']:2d}: measured {row['measured']:.2e} [1σ {row['lo1']:.2e}, {row['hi1']:.2e}], "
                  f"projection {row['gidney']:.2e}, ratio {row['measured'] / row['gidney']:.1f}")
    print(f"one shot at E[faults] ≤ {r['target']:.3f}: measured d = {r['baseline']['d']}, "
          f"Gidney's projection d = {r['gidney_projection']['d']} (published {r['published_d']})")
    for s in r["scenarios"]:
        print(f"  {s['category']:<14} {s['name']:<48} d = {s['d']}")
    for k, m in (r["noise_models"] or {}).items():
        print(f"  memory model {k:<8} Λ = {m['lambda']['z']:.1f}, {m['lambda']['x']:.1f}  d (worse basis, either, one shot) = "
              f"{m['d_basis']}, {m['d_any']}, {m['d_shot']}")
    if r["mr"]:
        m = r["mr"]
        print(f"one measure-and-reset layer: per-round d = {m['d_round']} (1σ {m['d_band']['-1']}–{m['d_band']['1']}), "
              f"one shot d = {m['one_shot']['d']}, Λ = " + ", ".join(f"{f['lambda']:.1f}" for f in m["fits"].values()))
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(r, f, indent=1)


if __name__ == "__main__":
    main()
