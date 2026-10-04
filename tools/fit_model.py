"""Fits the logical model to the raw calibration data and writes data/calibration/model.json
(components, per-p and global fits) and a human-readable summary."""

import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))

from ftalgo.calib.model import DATA, LogicalModel  # noqa: E402

m = LogicalModel.from_raw()
with open(os.path.join(DATA, "model.json"), "w") as f:
    json.dump(m.to_json(), f, indent=1)
with open(os.path.join(DATA, "summary.json"), "w") as f:
    json.dump(m.summary(), f, indent=1)
for name, fs in sorted(m.per_p.items()):
    print(name, {p: round(2.718281828 ** (-f.b), 2) for p, f in fs.items()})
