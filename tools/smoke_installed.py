"""Smoke test of an installed ftalgo (a wheel, not the repository): the bundled model loads, Shor
for N = 15 compiles onto the machine and runs, and the estimator answers. Run it from a directory
without the repository, so the import cannot fall back to python/ftalgo."""

import ftalgo
from ftalgo import scale, shor
from ftalgo.arch import FACTORIES, Architecture
from ftalgo.calib.model import LogicalModel
from ftalgo.engine import Program
from ftalgo.schedule import compile_noisy

model = LogicalModel.load()
c, info = shor.order_finding(15)
noisy = compile_noisy(c, Architecture(model, 11, 1e-3, FACTORIES["cultivation"]))
outcomes = info.outcomes(Program(noisy.text).sample(200, seed=1))
peak = shor.peak_table(15, info.a, info.m)
share = sum(bool(peak[o]) for o in outcomes) / len(outcomes)
d = scale.estimate(scale.GIDNEY_2025, model, target=scale.GIDNEY_2025_TARGET).d
print(f"ftalgo {ftalgo.engine_version} at {ftalgo.__file__}: N = 15, d = 11 peak probability {share:.3f}; Gidney 2025 d = {d}")
assert share > 0.8 and d is not None
