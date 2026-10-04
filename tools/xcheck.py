"""Cross-check: a subset of calibration experiments rerun with Stim (sampling) and PyMatching
(correlated matching), compared with stabilizer-qec's recorded rates.
Writes data/calibration/xcheck.json."""

import json
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))

import stabilizer_qec as sq  # noqa: E402

from ftalgo.calib.model import load_raw  # noqa: E402
from ftalgo.calib.qec import Counts, run_stim  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), "..", "data", "calibration", "xcheck.json")
raw = load_raw()
cases = []
for p in (0.002, 0.003, 0.005):
    for d in (3, 5, 7):
        cases += [
            (f"mem-z-d{d}-p{p}-T{3 * d}", lambda d=d, p=p: sq.memory_circuit(distance=d, rounds=3 * d, p=p, basis="z")),
            (f"cnot-z-d{d}-p{p}", lambda d=d, p=p: sq.surgery.cnot(d, merged=d, p=p, inputs="z")),
            (f"cnot-x-d{d}-p{p}", lambda d=d, p=p: sq.surgery.cnot(d, merged=d, p=p, inputs="x")),
            (f"zz-z-d{d}-p{p}", lambda d=d, p=p: sq.surgery.zz_measurement(d, merged=d, p=p, basis="z")),
        ]
out = []
for name, build in cases:
    ours = Counts.from_json(raw[name]["counts"])
    rate = ours.failures / ours.shots
    shots = int(min(max(4000 / max(rate, 1e-9), 1 << 16), 1 << 23))
    theirs = run_stim(str(build()), shots=shots, seed=12345)
    r2 = theirs.failures / theirs.shots
    sig = math.sqrt(rate * (1 - rate) / ours.shots + r2 * (1 - r2) / theirs.shots)
    z = (r2 - rate) / sig
    out.append({"name": name, "stabilizer_qec": {"failures": ours.failures, "shots": ours.shots}, "stim_pymatching": {"failures": theirs.failures, "shots": theirs.shots}, "z": z})
    print(f"{name}: ours {rate:.4e}  stim+pymatching {r2:.4e}  z={z:+.2f}", flush=True)
with open(OUT, "w") as f:
    json.dump({"note": "any-observable failure rate; correlated matching on both sides", "cases": out}, f, indent=1)
