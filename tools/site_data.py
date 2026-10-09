"""Gathers the committed data into what the site draws: site/data/site.json, and the browser
demo's precompiled noisy programs (Shor, N = 15) in site/data/demo/."""

from __future__ import annotations

import glob
import gzip
import json
import math
import os
import re
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
        if re.search(r" \d+\.json$", f):  # sync-conflict copies ("name 2.json")
            continue
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
    raw = sorted(f for f in glob.glob(os.path.join(ROOT, "data", "calibration", "raw", "*.json")) if not re.search(r" \d+\.json$", f))
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
    shor_modern_runs = [r for r in load_dir(os.path.join(ROOT, "data", "results", "shor_modern")) if "d" in r]
    for r in shor_modern_runs:
        r.pop("estimate", None)
    out["shor_modern"] = {
        "runs": shor_modern_runs,
        "baselines": [r for r in load_dir(os.path.join(ROOT, "data", "results", "shor_modern")) if "d" not in r]
    }
    sens = os.path.join(ROOT, "data", "results", "sensitivity.json")
    if os.path.exists(sens):
        out["sensitivity"] = json.load(open(sens))
    qh = os.path.join(ROOT, "data", "results", "qpe_hist.json")
    if os.path.exists(qh):
        out["qpe_hist"] = json.load(open(qh))
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

    bx = os.path.join(ROOT, "data", "results", "biased_xzzx.json")
    if os.path.exists(bx):
        out["biased_xzzx"] = json.load(open(bx))
    lih_file = os.path.join(ROOT, "data", "results", "lih.json")
    if os.path.exists(lih_file):
        out["lih"] = json.load(open(lih_file))
    else:
        import numpy as np
        from ftalgo.chem.lih import lih_dissociation_curve, lih_problem
        r_sweep = np.linspace(1.0, 3.0, 15)
        curve_pts = lih_dissociation_curve(r_sweep)
        prob = lih_problem(1.595, ncas=2)
        out["lih"] = {
            "molecule": "LiH",
            "r_equilibrium": 1.595,
            "ncas": 2,
            "qubits": prob.hamiltonian.n,
            "terms": len(prob.hamiltonian),
            "e_hf": prob.e_hf,
            "e_casci": prob.e_fci,
            "dissociation_curve": [
                {
                    "r_angstrom": round(p.r_angstrom, 4),
                    "r_bohr": round(p.r_bohr, 4),
                    "e_hf": p.e_hf,
                    "e_casci": p.e_casci,
                    "e_fci": p.e_fci,
                }
                for p in curve_pts
            ],
            "runs": [],
        }

    # Modern demo: Shor, N = 15, windowed + Gidney MBU
    from ftalgo import arith_modern as am
    os.makedirs(os.path.join(SITE, "demo_modern"), exist_ok=True)
    c_mod, info_mod = am.order_finding_modern(15)
    demo_mod = {"N": 15, "a": info_mod.a, "m": info_mod.m, "records": info_mod.records, "peaks": shor.peak_table(15, info_mod.a, info_mod.m).astype(int).tolist(),
                "factors": shor.success_table(15, info_mod.a, info_mod.m).astype(int).tolist(), "configs": [],
                "ideal": shor.ideal_distribution(15, info_mod.a, info_mod.m).tolist(), "logical_qubits": c_mod.num_qubits,
                "toffoli": c_mod.toffoli_count(), "t": c_mod.t_count()}
    for p in DEMO_P:
        for fac in DEMO_F:
            for d in DEMO_D:
                nz = compile_noisy(c_mod, Architecture(model, d, p, FACTORIES[fac]))
                text = nz.text
                fname = f"N15-{fac}-p{p}-d{d}.txt.gz"
                with open(os.path.join(SITE, "demo_modern", fname), "wb") as f:
                    f.write(gzip.compress(text.encode(), 9, mtime=0))
                prog = Program(text)
                demo_mod["configs"].append({
                    "file": fname, "d": d, "p": p, "factory": fac, "expected_faults": prog.expected_faults,
                    "p0": float(prog.fault_count_distribution(0)[0]), "rounds": nz.rounds, "qubits": nz.physical_qubits(),
                    "budget": dict(nz.budget), "extrapolated": bool(nz.arch.extrapolated),
                })
    out["demo_modern"] = demo_mod

    # Chemistry demo: H2 Ground State QPE (m=8 bits)
    import numpy as np
    from ftalgo import qpe
    from ftalgo.chem.qubit import molecule, qubit_problem
    os.makedirs(os.path.join(SITE, "demo_chem"), exist_ok=True)
    prob_h2 = qubit_problem(molecule("H2"))
    c_chem, info_chem = qpe.qpe_circuit(prob_h2.hamiltonian, prob_h2.hf_bits, m=8, tau=2 * math.pi, steps=2, e_ref=prob_h2.e_hf, order=4)
    M_chem = 1 << info_chem.m
    y_vals = np.arange(M_chem)
    phi_vals = np.where(y_vals / M_chem >= 0.5, y_vals / M_chem - 1.0, y_vals / M_chem)
    e_vals = info_chem.e_ref - 2 * math.pi * phi_vals / info_chem.tau
    acc_mask = (np.abs(e_vals - prob_h2.e_fci) < qpe.CHEMICAL_ACCURACY).astype(int).tolist()
    acc4_mask = (np.abs(e_vals - prob_h2.e_fci) < 4 * qpe.CHEMICAL_ACCURACY).astype(int).tolist()
    demo_chem = {
        "molecule": "H2",
        "m": info_chem.m,
        "tau": info_chem.tau,
        "e_hf": prob_h2.e_hf,
        "e_fci": prob_h2.e_fci,
        "resolution": info_chem.resolution(),
        "records": info_chem.records,
        "energies": e_vals.tolist(),
        "chemical_accuracy": acc_mask,
        "within_4x": acc4_mask,
        "logical_qubits": c_chem.num_qubits,
        "toffoli": c_chem.toffoli_count(),
        "t": c_chem.t_count(),
        "configs": [],
    }
    for p in DEMO_P:
        for fac in DEMO_F:
            for d in DEMO_D:
                nz = compile_noisy(c_chem, Architecture(model, d, p, FACTORIES[fac]))
                text = nz.text
                fname = f"H2-{fac}-p{p}-d{d}.txt.gz"
                with open(os.path.join(SITE, "demo_chem", fname), "wb") as f:
                    f.write(gzip.compress(text.encode(), 9, mtime=0))
                prog = Program(text)
                demo_chem["configs"].append({
                    "file": fname, "d": d, "p": p, "factory": fac, "expected_faults": prog.expected_faults,
                    "p0": float(prog.fault_count_distribution(0)[0]), "rounds": nz.rounds, "qubits": nz.physical_qubits(),
                    "budget": dict(nz.budget), "extrapolated": bool(nz.arch.extrapolated),
                })
    out["demo_chem"] = demo_chem

    with open(os.path.join(SITE, "site.json"), "w") as f:
        json.dump(clean(out), f, separators=(",", ":"), allow_nan=False, default=lambda o: float(o) if isinstance(o, (int, float)) else str(o))
    print("site.json", os.path.getsize(os.path.join(SITE, "site.json")) // 1024, "KiB;", len(demo["configs"]), "demo programs;", len(demo_mod["configs"]), "modern demo programs;", len(demo_chem["configs"]), "chem demo programs")


if __name__ == "__main__":
    main()
