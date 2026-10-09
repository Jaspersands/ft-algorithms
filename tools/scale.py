"""Scaling estimates: the textbook Shor circuit from the simulated sizes to RSA-2048, and the
same error model applied to Gidney 2025 (RSA-2048) and Lee et al. 2021 (FeMoco).
Writes data/results/scaling.json."""

import dataclasses
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))

from ftalgo import scale  # noqa: E402
from ftalgo.calib.model import LogicalModel  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), "..", "data", "results", "scaling.json")

model = LogicalModel.load()
st = scale.fit_structure(model)
st_mod = scale.fit_modern_structure(model)
print("textbook structure", st)
print("modern structure", st_mod)
rows = []
rows_mod = []
for factory in ("cultivation", "15to1"):
    for n in (4, 5, 6, 8, 10, 16, 32, 64, 128, 256, 512, 1024, 2048):
        e = scale.estimate(scale.textbook_shor(n, st), model, factory=factory)
        rows.append({"kind": "textbook_shor", "n": n, "factory": factory, **dataclasses.asdict(e)})
        print(f"textbook {factory} n={n}: d={e.d}, {e.physical_qubits:.3g} qubits, {e.seconds:.3g} s")

        e_mod = scale.estimate(scale.modern_shor(n, st_mod), model, factory=factory)
        rows_mod.append({"kind": "modern_shor", "n": n, "factory": factory, **dataclasses.asdict(e_mod)})
        print(f"modern   {factory} n={n}: d={e_mod.d}, {e_mod.physical_qubits:.3g} qubits, {e_mod.seconds:.3g} s")

published = []
for w in (scale.GIDNEY_2025, scale.FEMOCO_THC):
    e = scale.estimate(w, model)
    published.append({"source": w.source, **dataclasses.asdict(w), "estimate": dataclasses.asdict(e)})
    print(w.name, e.d, f"{e.physical_qubits:.3g} qubits", f"{e.seconds / 86400:.2f} days")
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w") as f:
    json.dump({
        "structure": dataclasses.asdict(st),
        "structure_modern": dataclasses.asdict(st_mod),
        "target_expected_faults": 0.1,
        "p": 1e-3,
        "textbook": rows,
        "modern": rows_mod,
        "published": published
    }, f, indent=1)
