"""The algorithm grid (textbook and windowed Shor, phase estimation, LiH) under a chosen calibration,
split into jobs that run in parallel (.github/workflows/grid.yml runs it on GitHub Actions).

    python tools/grid.py plan --model data/calibration/model_mr.json --out data/results_mr --jobs 60 --parts parts > plan.json
    python tools/grid.py run --plan plan.json --job 3 --parts parts     # one job's chunks, one JSON each
    python tools/grid.py merge --plan plan.json --parts parts           # chunks -> one result per configuration

The grid has the (N or molecule, p, factory) groups of data/results/, with each group's distances
chosen afresh under the model (experiment.d_range, as the run_* tools do). A configuration whose
run is expected to be long is sampled in chunks (seeds differ, shots split evenly) that merge into
one result (experiment.merge_runs). The expected cost of a run comes from the data/results/ run of
the same group with the nearest expected number of faults. Chunks already in --parts are skipped,
so a second plan after a partial run plans only what is missing.
"""

from __future__ import annotations

import argparse
import glob
import json
import math
import os
import shutil
import subprocess
import sys
import zlib

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "python"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from ftalgo.arch import FACTORIES, Architecture  # noqa: E402
from ftalgo.calib.model import LogicalModel  # noqa: E402
from ftalgo.engine import Program  # noqa: E402
from ftalgo.experiment import d_range, merge_runs, run_config  # noqa: E402
from ftalgo.schedule import compile_noisy  # noqa: E402

KINDS = {"shor": ("run_shor", 4000, 200), "shor_modern": ("run_shor_modern", 2000, 200), "qpe": ("run_qpe", 1000, 200)}
CHUNK_SECONDS = 1800.0  # expected cost of one chunk, in the seconds of the data/results/ runs
LIH_SECONDS = 3000.0


def builder(kind):
    return __import__(KINDS[kind][0]).build


def old_runs(kind):
    out = []
    for f in glob.glob(os.path.join(ROOT, "data", "results", kind, "*-d*.json")):
        r = json.load(open(f))
        out.append(r)
    return out


def groups(kind):
    return sorted({(r.get("N") or r.get("molecule"), r["p"], r["factory"]) for r in old_runs(kind)}, key=str)


def cost(old, key, ef):
    same = [r for r in old if (r.get("N") or r.get("molecule"), r["p"], r["factory"]) == key]
    if not same:
        return 600.0
    near = min(same, key=lambda r: abs(math.log(max(r["expected_faults"], 1e-9)) - math.log(max(ef, 1e-9))))
    return float(near["wall_seconds"])


def plan(args):
    model = LogicalModel.load(args.model)
    items = []
    for kind in KINDS:
        old = old_runs(kind)
        build = builder(kind)
        built = {}
        for key, p, fac in groups(kind):
            if key not in built:
                built[key] = build(key)
            circ = built[key][0]
            kw = __import__(KINDS[kind][0]).D_RANGE if kind == "qpe" else {}
            for d in d_range(circ, model, p, fac, **kw):
                ef = Program(compile_noisy(circ, Architecture(model, d, p, FACTORIES[fac])).text).expected_faults
                c = cost(old, (key, p, fac), ef)
                tag = str(key).replace("+", "p")
                name = f"N{key}-{fac}-p{p}-d{d}" if kind != "qpe" else f"{tag}-{fac}-p{p}-d{d}"
                items.append({"kind": kind, "key": key, "p": p, "factory": fac, "d": d, "name": name,
                              "expected_faults": ef, "cost": c, "chunks": max(1, math.ceil(c / CHUNK_SECONDS))})
                print(f"{kind} {name}: E[faults] {ef:.3g}, cost {c:.0f} s, {items[-1]['chunks']} chunk(s)", file=sys.stderr, flush=True)
    items.append({"kind": "lih", "name": "lih", "cost": LIH_SECONDS, "chunks": 1})
    units = [(i, c, it["cost"] / it["chunks"]) for i, it in enumerate(items) for c in range(it["chunks"])
             if not os.path.exists(part_path(args.parts, it, c))]
    # Longest first, each to the least loaded job.
    jobs = [[] for _ in range(max(1, min(args.jobs, len(units))))]
    load = [0.0] * len(jobs)
    for i, c, w in sorted(units, key=lambda u: -u[2]):
        j = load.index(min(load))
        jobs[j].append([i, c])
        load[j] += w
    print(f"{len(items)} configurations, {len(units)} chunks to run in {len(jobs)} jobs; "
          f"heaviest job {max(load) / 3600:.1f} h, total {sum(load) / 3600:.1f} h (data/results/ seconds)", file=sys.stderr)
    json.dump({"model": args.model, "out": args.out, "items": items, "jobs": jobs}, sys.stdout, indent=1)


def part_path(parts, item, c):
    return os.path.join(parts, f"{item['kind']}-{item['name']}-c{c}.json")


def run(args):
    pl = json.load(open(args.plan))
    model = LogicalModel.load(pl["model"])
    os.makedirs(args.parts, exist_ok=True)
    built = {}
    for i, c in pl["jobs"][args.job]:
        it = pl["items"][i]
        path = part_path(args.parts, it, c)
        if os.path.exists(path):
            continue
        if it["kind"] == "lih":
            # run_lih.py saves after each distance; only a finished run counts as done.
            subprocess.run([sys.executable, os.path.join(ROOT, "tools", "run_lih.py"), "--model", pl["model"], "--out", path + ".partial"], check=True)
            os.replace(path + ".partial", path)
            continue
        k = (it["kind"], it["key"])
        if k not in built:
            built[k] = builder(it["kind"])(it["key"])
        circ, score = built[k][0], built[k][2]
        _, shots, min_shots = KINDS[it["kind"]]
        n = it["chunks"]
        r = run_config(circ, score, model, d=it["d"], p=it["p"], factory=it["factory"], seed=zlib.crc32(f"{it['name']}#c{c}".encode()),
                       shots=math.ceil(shots / n), min_shots=math.ceil(min_shots / n))
        r.update({"N": it["key"]} if it["kind"] != "qpe" else {"molecule": it["key"]})
        if it["kind"] != "qpe":
            r["a"] = built[k][1].a
        r["algorithm"] = it["kind"]
        with open(path + ".tmp", "w") as f:
            json.dump(r, f, indent=1)
        os.replace(path + ".tmp", path)
        print(f"{it['kind']} {it['name']} chunk {c + 1}/{n}: {r['scores'][0]:.3f} ± {r['score_sigmas'][0]:.3f} ({r['wall_seconds']} s)", flush=True)


def merge(args):
    pl = json.load(open(args.plan))
    out = args.out or pl["out"]
    missing = 0
    for kind in KINDS:
        os.makedirs(os.path.join(out, kind), exist_ok=True)
        # Baselines (noiseless and random scores) do not depend on the model.
        for f in glob.glob(os.path.join(ROOT, "data", "results", kind, "*-baselines.json")):
            shutil.copy(f, os.path.join(out, kind))
    for it in pl["items"]:
        paths = [part_path(args.parts, it, c) for c in range(it["chunks"])]
        have = [p for p in paths if os.path.exists(p)]
        if len(have) < len(paths):
            missing += 1
            continue
        if it["kind"] == "lih":
            shutil.copy(have[0], os.path.join(out, "lih.json"))
            continue
        r = merge_runs([json.load(open(p)) for p in have])
        with open(os.path.join(out, it["kind"], it["name"] + ".json"), "w") as f:
            json.dump(r, f, indent=1)
    print(f"merged {len(pl['items']) - missing} of {len(pl['items'])} configurations into {out}; {missing} incomplete")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("plan")
    a.add_argument("--model", default=os.path.join("data", "calibration", "model_mr.json"))
    a.add_argument("--out", default=os.path.join("data", "results_mr"))
    a.add_argument("--jobs", type=int, default=60)
    a.add_argument("--parts", default="parts")
    a = sub.add_parser("run")
    a.add_argument("--plan", required=True)
    a.add_argument("--job", type=int, required=True)
    a.add_argument("--parts", default="parts")
    a = sub.add_parser("merge")
    a.add_argument("--plan", required=True)
    a.add_argument("--parts", default="parts")
    a.add_argument("--out", default=None, help="instead of the plan's output directory")
    args = ap.parse_args()
    {"plan": plan, "run": run, "merge": merge}[args.cmd](args)


if __name__ == "__main__":
    main()
