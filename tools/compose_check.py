"""Composition check over the grid: whole surgery experiments predicted (ftsim, logical level)
from calibrated parts vs measured (stabilizer-qec, circuit level). Writes
data/calibration/composition.json."""

import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))

from ftalgo.calib.compose import compare  # noqa: E402
from ftalgo.calib.model import load_raw  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), "..", "data", "calibration", "composition.json")

raw = load_raw()
res = []
for p in (0.001, 0.002, 0.003, 0.005):
    for d in (3, 5, 7):
        r = compare(raw, d, p, shots=1 << 21, seed=d * 1000 + int(p * 1e4))
        if r:
            res.append({"d": d, "p": p, "experiments": r})
            print(d, p, {k: round(v["any"]["predicted"] / v["any"]["measured"], 3) for k, v in r.items()}, flush=True)
with open(OUT, "w") as f:
    json.dump({"note": "predicted from calibrated parts by ftsim (logical level) vs measured by stabilizer-qec (circuit level); 2^21 predicted shots each", "results": res}, f, indent=1)
