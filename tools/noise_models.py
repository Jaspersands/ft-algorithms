"""Surface-code memories at p = 0.1% under two more circuit noise models, to set beside this
project's SD6 calibration and locate the gap to Gidney 2025's projection of surface-code memory:

- si1000: the superconducting-inspired model of the yoked-surface-code simulations that projection
  descends from, on stabilizer-qec's memory rewritten over CZ;
- uniform: Stim's generated memory with every noise parameter at p ("uniform depolarizing", the
  label of Gidney 2025's Figure 6), which unlike SD6 leaves qubits idle during gate layers noiseless.

Each is held for T = d and 3d rounds in both bases and decoded by correlated matching, exactly as the
calibration's memories are (``ftalgo.calib.circuit_noise``); one JSON file per experiment in
data/calibration/noise_models/, resumable. ``--summary`` fits all three models and writes
data/calibration/noise_models.json.

    python tools/noise_models.py               # run everything missing
    python tools/noise_models.py --summary
"""

from __future__ import annotations

import argparse
import fnmatch
import glob
import json
import math
import os
import sys
import time
import zlib

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))

from ftalgo.calib.qec import Counts, run  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
RAW = os.path.join(ROOT, "data", "calibration", "noise_models")
CAL = os.path.join(ROOT, "data", "calibration", "raw")
OUT = os.path.join(ROOT, "data", "calibration", "noise_models.json")
P = 0.001
DS = (3, 5, 7, 9, 11)
MODELS = ("si1000", "uniform")
TARGET = 300
CAP_SECONDS = 1800


def experiments():
    from ftalgo.calib.circuit_noise import memory

    exps = []
    for d in DS:
        for model in MODELS:
            for b in ("z", "x"):
                for T in (d, 3 * d):
                    exps.append((f"{model}-mem-{b}-d{d}-p{P}-T{T}", dict(kind="memory", model=model, basis=b, d=d, p=P, rounds=T),
                                 lambda d=d, T=T, b=b, model=model: memory(d, T, b, model, P)))
    return exps


def per_round(rec_d: dict, rec_3d: dict) -> tuple[float, float]:
    """(per-round flip probability, its 1σ) from memories of d and 3d rounds: the ratio of their
    fidelities removes preparation and readout, as in the calibration."""
    out = []
    for r in (rec_d, rec_3d):
        c = Counts.from_json(r["counts"])
        q = c.failures / c.shots
        out.append((1 - 2 * q, math.sqrt(q * (1 - q) / c.shots)))
    (f1, s1), (f3, s3) = out
    T = rec_3d["rounds"] - rec_d["rounds"]
    fr = (f3 / f1) ** (1 / T)
    # σ of log(f3/f1)/T, from the binomial σ of each flip probability.
    sl = math.sqrt((2 * s1 / f1) ** 2 + (2 * s3 / f3) ** 2) / T
    return (1 - fr) / 2, fr * sl / 2


def load(path_glob: str) -> dict[str, dict]:
    out = {}
    for f in sorted(glob.glob(path_glob)):
        with open(f) as fh:
            r = json.load(fh)
        if os.path.basename(f) == r["name"] + ".json":
            out[r["name"]] = r
    return out


def summary() -> dict:
    new = load(os.path.join(RAW, "*.json"))
    sd = load(os.path.join(CAL, "mem-*.json"))
    models = {}
    for model, recs, prefix in (("sd6", sd, "mem"), ("si1000", new, "si1000-mem"), ("uniform", new, "uniform-mem")):
        rows = []
        for d in DS:
            row = {"d": d}
            for b in ("z", "x"):
                a, c = recs.get(f"{prefix}-{b}-d{d}-p{P}-T{d}"), recs.get(f"{prefix}-{b}-d{d}-p{P}-T{3 * d}")
                if a and c:
                    q, s = per_round(a, c)
                    row[b] = {"q": q, "sigma": s, "failures": [Counts.from_json(r["counts"]).failures for r in (a, c)],
                              "shots": [Counts.from_json(r["counts"]).shots for r in (a, c)]}
            if "z" in row and "x" in row:
                rows.append(row)
        fits = {}
        for b in ("z", "x"):
            pts = [(r["d"], r[b]["q"], min(r[b]["failures"])) for r in rows if min(r[b]["failures"]) >= 10]
            if len(pts) >= 2:
                xs = [(d + 1) / 2 for d, _, _ in pts]
                ys = [math.log(q) for _, q, _ in pts]
                ws = [w for _, _, w in pts]
                sw = sum(ws)
                mx, my = sum(w * x for w, x in zip(ws, xs)) / sw, sum(w * y for w, y in zip(ws, ys)) / sw
                slope = sum(w * (x - mx) * (y - my) for w, x, y in zip(ws, xs, ys)) / sum(w * (x - mx) ** 2 for w, x in zip(ws, xs))
                fits[b] = {"a": my - slope * mx, "b": slope, "lambda": math.exp(-slope), "dmax": max(d for d, _, _ in pts)}
        models[model] = {"rows": rows, "fits": fits}
    return {"p": P, "decoder": "correlated matching (stabilizer_qec.Matching)",
            "basis_note": "z: Z-basis memory, flips of the Z logical (X or Y errors); x: the X logical (Z or Y)",
            "models": models}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--summary", action="store_true")
    ap.add_argument("--only", default="*")
    ap.add_argument("--target", type=int, default=TARGET)
    ap.add_argument("--cap", type=float, default=CAP_SECONDS)
    args = ap.parse_args()
    if args.summary:
        with open(OUT, "w") as f:
            json.dump(summary(), f, indent=1)
        return
    import stabilizer_qec as sq
    import subprocess

    os.makedirs(RAW, exist_ok=True)
    sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True, cwd=ROOT).stdout.strip()
    vers = {"stabilizer_qec": sq.__version__, "ftalgo_git": sha, "decoder": "stabilizer_qec.Matching(enable_correlations=True)"}
    for name, meta, build in experiments():
        if not fnmatch.fnmatch(name, args.only):
            continue
        path = os.path.join(RAW, name + ".json")
        if args.list or os.path.exists(path):
            print(("done " if os.path.exists(path) else "todo ") + name)
            continue
        seed = zlib.crc32(name.encode())
        circuit = str(build())
        t0 = time.time()
        total, batch = None, 1 << 14
        while True:
            res = run(circuit, shots=batch, seed=(seed + (total.shots if total else 0)) & 0xFFFFFFFF, batch=batch)
            total = res if total is None else total.merge(res)
            if total.failures >= args.target or time.time() - t0 > args.cap:
                break
            batch = min(batch * 2, 1 << 22)
        rec = {"name": name, **meta, "counts": total.to_json(), "seed": seed, "seed_note": "batches seeded seed + shots done so far",
               "versions": vers, "command": "python tools/noise_models.py", "wall_seconds": round(time.time() - t0, 1)}
        with open(path + ".tmp", "w") as f:
            json.dump(rec, f, indent=1)
        os.replace(path + ".tmp", path)
        print(f"{name}: {total.failures} failures / {total.shots} shots ({time.time() - t0:.0f} s)", flush=True)


if __name__ == "__main__":
    main()
