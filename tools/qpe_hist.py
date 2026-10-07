"""Energy estimates of individual phase-estimation runs (H2, cultivated states, p = 0.1%) at a
few distances: the distribution of the answer itself. Writes data/results/qpe_hist.json."""

import json
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))

from ftalgo import qpe  # noqa: E402
from ftalgo.arch import FACTORIES, Architecture  # noqa: E402
from ftalgo.calib.model import LogicalModel  # noqa: E402
from ftalgo.chem.qubit import molecule, qubit_problem  # noqa: E402
from ftalgo.engine import Program  # noqa: E402
from ftalgo.schedule import compile_noisy  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), "..", "data", "results", "qpe_hist.json")
model = LogicalModel.load()
prob = qubit_problem(molecule("H2"))
circ, info = qpe.qpe_circuit(prob.hamiltonian, prob.hf_bits, m=10, tau=2 * math.pi, steps=4, e_ref=prob.e_hf, order=4)
out = {"molecule": "H2", "e_fci": prob.e_fci, "e_hf": prob.e_hf, "resolution": info.resolution(), "runs": []}
for d in (13, 15, 17, 21):
    nz = compile_noisy(circ, Architecture(model, d, 0.001, FACTORIES["cultivation"]))
    rows = Program(nz.text).sample(400, seed=100 + d)
    e = info.energies(rows)
    out["runs"].append({"d": d, "energies": [float(x) for x in e]})
    print(d, f"{(abs(e - prob.e_fci) < qpe.CHEMICAL_ACCURACY).mean():.3f}", flush=True)
noiseless = Program(circ.to_text()).sample(400, seed=7, faults="none")
out["noiseless"] = [float(x) for x in info.energies(noiseless)]
with open(OUT, "w") as f:
    json.dump(out, f)
