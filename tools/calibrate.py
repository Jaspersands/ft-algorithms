"""Runs the calibration grid on stabilizer-qec and records raw flip-pattern counts.

Each experiment is one JSON file in data/calibration/raw/, written when it finishes, so the run
can be interrupted and resumed (finished experiments are skipped). Cheap experiments run first.

    python tools/calibrate.py               # run everything missing
    python tools/calibrate.py --list        # show the plan
    python tools/calibrate.py --only 'mem-*-d3-*'
"""

from __future__ import annotations

import argparse
import fnmatch
import json
import os
import subprocess
import sys
import time
import zlib

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))

from ftalgo.calib.choi import choi_memory  # noqa: E402
from ftalgo.calib.qec import run  # noqa: E402

RAW = os.path.join(os.path.dirname(__file__), "..", "data", "calibration", "raw")
PS = [0.001, 0.002, 0.003, 0.005]
TARGET = 500        # failures (shots with any flipped observable)
CAP_SECONDS = 300   # per experiment
MAX_SHOTS = 1 << 31


def experiments():
    import stabilizer_qec as sq

    exps = []
    for p in PS:
        for d in (3, 5, 7, 9, 11):
            for b in ("z", "x"):
                for T in (d, 2 * d, 3 * d, 4 * d):
                    exps.append((f"mem-{b}-d{d}-p{p}-T{T}", dict(kind="memory", basis=b, d=d, p=p, rounds=T),
                                 lambda d=d, T=T, p=p, b=b: sq.memory_circuit(distance=d, rounds=T, p=p, basis=b)))
            if d <= 7:
                for T in (d, 3 * d):
                    exps.append((f"choi-d{d}-p{p}-T{T}", dict(kind="choi", d=d, p=p, rounds=T),
                                 lambda d=d, T=T, p=p: choi_memory(d, T, p)))
            if d <= 9 or p >= 0.003:
                for inp in ("z", "x"):
                    exps.append((f"cnot-{inp}-d{d}-p{p}", dict(kind="cnot", inputs=inp, d=d, p=p, merged=d, rounds=4 * d),
                                 lambda d=d, p=p, inp=inp: sq.surgery.cnot(d, merged=d, p=p, inputs=inp)))
                for b in ("z", "x"):
                    exps.append((f"zz-{b}-d{d}-p{p}", dict(kind="zz", basis=b, d=d, p=p, merged=d, rounds=3 * d),
                                 lambda d=d, p=p, b=b: sq.surgery.zz_measurement(d, merged=d, p=p, basis=b)))
                    exps.append((f"xx-{b}-d{d}-p{p}", dict(kind="xx", basis=b, d=d, p=p, merged=d, rounds=3 * d),
                                 lambda d=d, p=p, b=b: sq.surgery.xx_measurement(d, merged=d, p=p, basis=b)))
            if d <= 7:
                exps.append((f"rep-k3-d{d}-p{p}", dict(kind="repeated_zz", k=3, d=d, p=p, merged=d, rounds=d + 3 * (d + 1) + d - 1),
                             lambda d=d, p=p: sq.surgery.repeated_zz(d, k=3, merged=d, p=p)))
                exps.append((f"line-n3-d{d}-p{p}", dict(kind="line", n=3, d=d, p=p, merged=d, rounds=3 * d),
                             lambda d=d, p=p: sq.surgery.line(d, n=3, merged=d, p=p)))
    # Cheap first: small d, high p.
    exps.sort(key=lambda e: (e[1]["d"], -e[1]["p"], e[1].get("rounds", 0)))
    return exps


def versions() -> dict:
    import stabilizer_qec as sq

    try:
        sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True, cwd=os.path.dirname(__file__)).stdout.strip()
    except OSError:
        sha = ""
    return {"stabilizer_qec": sq.__version__, "ftalgo_git": sha, "decoder": "stabilizer_qec.Matching(enable_correlations=True)"}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--only", default="*")
    ap.add_argument("--target", type=int, default=TARGET)
    ap.add_argument("--cap", type=float, default=CAP_SECONDS)
    args = ap.parse_args()
    os.makedirs(RAW, exist_ok=True)
    vers = versions()
    todo = [e for e in experiments() if fnmatch.fnmatch(e[0], args.only)]
    for name, meta, build in todo:
        path = os.path.join(RAW, name + ".json")
        done = os.path.exists(path)
        if args.list:
            print(("done " if done else "todo ") + name)
            continue
        if done:
            continue
        seed = zlib.crc32(name.encode())
        circuit = build()
        t0 = time.time()
        # Grow the batch so a low-rate experiment amortizes overhead; stop at the time cap.
        total = None
        batch = 1 << 14
        while True:
            res = run(circuit, shots=batch, seed=(seed + (total.shots if total else 0)) & 0xFFFFFFFF, batch=batch)
            total = res if total is None else total.merge(res)
            if total.failures >= args.target or time.time() - t0 > args.cap or total.shots >= MAX_SHOTS:
                break
            batch = min(batch * 2, 1 << 22)
        rec = {
            "name": name,
            **meta,
            "counts": total.to_json(),
            "seed": seed,
            "seed_note": "batches seeded seed + shots done so far",
            "versions": vers,
            "command": "python tools/calibrate.py",
            "wall_seconds": round(time.time() - t0, 1),
        }
        with open(path + ".tmp", "w") as f:
            json.dump(rec, f, indent=1)
        os.replace(path + ".tmp", path)
        print(f"{name}: {total.failures} failures / {total.shots} shots ({time.time() - t0:.0f} s)", flush=True)


if __name__ == "__main__":
    main()
