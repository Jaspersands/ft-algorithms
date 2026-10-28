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
    r = run(LogicalModel.load())
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
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(r, f, indent=1)


if __name__ == "__main__":
    main()
