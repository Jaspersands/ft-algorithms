"""Gathers the committed data into what the site draws: site/data/site.json, and the browser
demo's precompiled noisy programs (Shor, N = 15) in site/data/demo/."""

from __future__ import annotations

import glob
import gzip
import json
import math
import os
import sys

ROOT = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, os.path.join(ROOT, "python"))

from ftalgo import shor  # noqa: E402
from ftalgo.arch import FACTORIES, Architecture  # noqa: E402
from ftalgo.calib.model import LogicalModel  # noqa: E402
from ftalgo.engine import Program  # noqa: E402
from ftalgo.schedule import compile_noisy  # noqa: E402

SITE = os.path.join(ROOT, "site", "data")
DEMO_D = [7, 9, 11, 13, 15]
DEMO_P = [0.001, 0.002]
DEMO_F = ["cultivation", "15to1", "injected"]


def clean(o):
    """JSON has no Infinity or NaN: non-finite numbers become null."""
    if isinstance(o, float) and not math.isfinite(o):
        return None
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    return o


def load_dir(path):
    out = []
    for f in sorted(glob.glob(os.path.join(path, "*.json"))):
        with open(f) as fh:
            out.append(json.load(fh))
    return out


def calibration(model: LogicalModel) -> dict:
    comps = {}
    for name, pts in model.points.items():
        comps[name] = {
            "points": [{"d": q.d, "p": q.p, "v": q.value, "s": q.sigma, "n": q.events} for q in sorted(pts, key=lambda q: (q.p, q.d))],
            "fits": {repr(p): {"a": f.a, "b": f.b, "dmax": f.dmax} for p, f in model.per_p.get(name, {}).items()},
        }
    return comps


def main():
    os.makedirs(os.path.join(SITE, "demo"), exist_ok=True)
    model = LogicalModel.load()
    out = {"calibration": calibration(model)}
    for name in ("composition", "xcheck"):
        path = os.path.join(ROOT, "data", "calibration", f"{name}.json")
        if os.path.exists(path):
            out[name] = json.load(open(path))
    raw = sorted(glob.glob(os.path.join(ROOT, "data", "calibration", "raw", "*.json")))
    out["calibration_meta"] = {
        "experiments": len(raw),
        "shots": sum(json.load(open(f))["counts"]["shots"] for f in raw),
        "failures": sum(sum(json.load(open(f))["counts"]["patterns"].values()) for f in raw),
    }
    shor_runs = [r for r in load_dir(os.path.join(ROOT, "data", "results", "shor")) if "d" in r]
    for r in shor_runs:
        r.pop("estimate", None)
    out["shor"] = {"runs": shor_runs, "baselines": [r for r in load_dir(os.path.join(ROOT, "data", "results", "shor")) if "d" not in r]}
    qpe_runs = [r for r in load_dir(os.path.join(ROOT, "data", "results", "qpe")) if "d" in r]
    for r in qpe_runs:
        r.pop("estimate", None)
    out["qpe"] = {"runs": qpe_runs, "baselines": [r for r in load_dir(os.path.join(ROOT, "data", "results", "qpe")) if "d" not in r]}
    sc = os.path.join(ROOT, "data", "results", "scaling.json")
    if os.path.exists(sc):
        out["scaling"] = json.load(open(sc))
    out["factories"] = {k: {"label": f.label, "source": f.source, "eps_t_1e-3": f.eps_t(1e-3), "eps_ccz_1e-3": f.eps_ccz(1e-3)} for k, f in FACTORIES.items()}

    # The browser demo: Shor, N = 15, every (d, p, factory) of the grid, gzipped.
    c, info = shor.order_finding(15)
    demo = {"N": 15, "a": info.a, "m": info.m, "records": info.records, "peaks": shor.peak_table(15, info.a, info.m).astype(int).tolist(),
            "factors": shor.success_table(15, info.a, info.m).astype(int).tolist(), "configs": [],
            "ideal": shor.ideal_distribution(15, info.a, info.m).tolist(), "logical_qubits": c.num_qubits,
            "toffoli": c.toffoli_count(), "t": c.t_count()}
    for p in DEMO_P:
        for fac in DEMO_F:
            for d in DEMO_D:
                nz = compile_noisy(c, Architecture(model, d, p, FACTORIES[fac]))
                text = nz.text
                fname = f"N15-{fac}-p{p}-d{d}.txt.gz"
                with open(os.path.join(SITE, "demo", fname), "wb") as f:
                    f.write(gzip.compress(text.encode(), 9, mtime=0))
                prog = Program(text)
                demo["configs"].append({
                    "file": fname, "d": d, "p": p, "factory": fac, "expected_faults": prog.expected_faults,
                    "p0": float(prog.fault_count_distribution(0)[0]), "rounds": nz.rounds, "qubits": nz.physical_qubits(),
                    "budget": dict(nz.budget), "extrapolated": bool(nz.arch.extrapolated),
                })
    out["demo"] = demo
    with open(os.path.join(SITE, "site.json"), "w") as f:
        json.dump(clean(out), f, separators=(",", ":"), allow_nan=False, default=lambda o: float(o) if isinstance(o, (int, float)) else str(o))
    print("site.json", os.path.getsize(os.path.join(SITE, "site.json")) // 1024, "KiB;", len(demo["configs"]), "demo programs")


if __name__ == "__main__":
    main()
