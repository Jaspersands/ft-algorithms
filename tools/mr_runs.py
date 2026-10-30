"""The algorithm runs under both calibrations, side by side: for each workload (cultivated magic
states), the smallest simulated distance whose success reaches 90% of the noiseless value, from
data/results/ (separate measure and reset layers) and data/results_mr/ (one measure-and-reset
layer; tools/grid.py on GitHub Actions). Used by readme_tables.py, site_data.py and report.py."""

from __future__ import annotations

import glob
import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent
CAL = {"separate": ROOT / "data" / "results", "mr": ROOT / "data" / "results_mr"}


def _runs(root, kind):
    out, base = {}, {}
    for f in glob.glob(str(root / kind / "*.json")):
        r = json.load(open(f))
        key = r.get("N") or r.get("molecule")
        if "d" not in r:
            base[key] = r
        else:
            out.setdefault((key, r["p"], r["factory"]), []).append(r)
    return out, base


def _d90(runs, noiseless):
    hit = next((r for r in sorted(runs, key=lambda r: r["d"]) if r["scores"][0] >= 0.9 * noiseless), None)
    return hit["d"] if hit else None


def rows() -> list[dict]:
    """One row per workload and p: noiseless success and the d reaching 90% of it under each
    calibration (None where no simulated distance does). Empty if data/results_mr/ is absent."""
    if not CAL["mr"].is_dir():
        return []
    out = []
    for kind, label in (("shor", "Shor, textbook, N = {}"), ("shor_modern", "Shor, windowed + MBU, N = {}"), ("qpe", "{}")):
        runs = {c: _runs(root, kind) for c, root in CAL.items()}
        base = runs["separate"][1]
        for key, p, fac in sorted(runs["separate"][0], key=lambda k: (str(type(k[0])), k[0], k[1])):
            if fac != "cultivation":
                continue
            nl = base[key]["noiseless"][0]
            d = {c: _d90(runs[c][0].get((key, p, fac), []), nl) for c in CAL}
            out.append({"workload": label.format(key), "p": p, "noiseless": nl, "d": d["separate"], "d_mr": d["mr"]})
    lih = {c: json.load(open(root / "lih.json")) for c, root in CAL.items()}
    nl = lih["separate"]["noiseless"][0]
    d = {c: _d90([r for r in lih[c]["runs"] if r["factory"] == "cultivation" and r["p"] == 0.001], nl) for c in CAL}
    out.append({"workload": "LiH, CAS(2, 2)", "p": 0.001, "noiseless": nl, "d": d["separate"], "d_mr": d["mr"]})
    return out


def summary(rs: list[dict]) -> dict:
    both = [r for r in rs if r["d"] is not None and r["d_mr"] is not None]
    return {"rows": len(rs), "compared": len(both), "lower": sum(r["d_mr"] < r["d"] for r in both),
            "same": sum(r["d_mr"] == r["d"] for r in both), "higher": sum(r["d_mr"] > r["d"] for r in both),
            "step": sorted({r["d"] - r["d_mr"] for r in both if r["d_mr"] < r["d"]})}


if __name__ == "__main__":
    rs = rows()
    for r in rs:
        print(r)
    print(summary(rs))
