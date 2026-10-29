"""Fits the logical model to the raw calibration data and writes data/calibration/model.json
(components, per-p and global fits) and a human-readable summary. With --mr, the calibration on
circuits with one measure-and-reset layer (raw_mr/) → model_mr.json and summary_mr.json."""

import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))

from ftalgo.calib.model import DATA, LogicalModel, load_raw  # noqa: E402

suffix = "_mr" if "--mr" in sys.argv else ""
m = LogicalModel.from_raw(load_raw(os.path.join(DATA, "raw" + suffix)))
with open(os.path.join(DATA, f"model{suffix}.json"), "w") as f:
    json.dump(m.to_json(), f, indent=1)
with open(os.path.join(DATA, f"summary{suffix}.json"), "w") as f:
    json.dump(m.summary(), f, indent=1)
for name, fs in sorted(m.per_p.items()):
    print(name, {p: round(2.718281828 ** (-f.b), 2) for p, f in fs.items()})
