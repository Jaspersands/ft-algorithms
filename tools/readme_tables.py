"""Fills the README's results section from the committed data (between the results markers).
With --check, fails if the README is out of date instead of writing it."""

import glob
import json
import math
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent


def real(paths):
    """Drops sync-conflict copies ("name 2.json")."""
    return [f for f in paths if not re.search(r" \d+\.json$", f)]


def jl(p):
    return json.loads(pathlib.Path(p).read_text())


def fmt_time(s):
    for lim, div, unit in ((1, 1e-3, "ms"), (120, 1, "s"), (7200, 60, "min"), (2 * 86400, 3600, "h"), (400 * 86400, 86400, "days")):
        if s < lim:
            return f"{s / div:.3g} {unit}"
    return f"{s / 86400 / 365.25:.3g} years"


def fmt_q(x):
    return f"{x / 1e6:.3g} M" if x >= 1e6 else f"{x / 1e3:.3g} k"


def sci(x):
    e = math.floor(math.log10(x))
    return f"{x / 10**e:.1f}e{e}"


def build() -> str:
    out = []
    raw = real(glob.glob(str(ROOT / "data/calibration/raw/*.json")))
    comp = jl(ROOT / "data/calibration/composition.json")
    ratios = sorted(e["any"]["predicted"] / e["any"]["measured"] for r in comp["results"] for e in r["experiments"].values())
    xc = jl(ROOT / "data/calibration/xcheck.json")
    shots = sum(jl(f)["counts"]["shots"] for f in raw)
    out.append(f"**Calibration**: {len(raw)} circuit-level experiments, {shots / 1e9:.2f} billion shots. "
               f"Whole surgery experiments predicted from calibrated parts: {ratios[0]:.2f}–{ratios[-1]:.2f}× the measured "
               f"failure rate (median {ratios[len(ratios) // 2]:.2f}×). Stim + PyMatching cross-check: {len(xc['cases'])} of "
               f"{len(xc['cases'])} within {max(abs(c['z']) for c in xc['cases']):.1f}σ.\n")
    rs = [jl(f) for f in real(glob.glob(str(ROOT / "data/results/shor/*.json")))]
    base = {r["N"]: r for r in rs if "d" not in r}
    runs = [r for r in rs if "d" in r]
    out.append("**Shor** (p = 0.1%, cultivated magic states; peak probability = share of runs on one of the r ideal peaks):\n")
    out.append("| N | logical qubits | Toffolis | noiseless peak | d for 90% of it | physical qubits | run time |\n|---:|---:|---:|---:|---:|---:|---:|")
    for N in sorted(base):
        b = base[N]
        rr = sorted((r for r in runs if r["N"] == N and r["factory"] == "cultivation" and r["p"] == 0.001), key=lambda r: r["d"])
        hit = next((r for r in rr if r["scores"][0] >= 0.9 * b["noiseless"][0]), None)
        out.append(f"| {N} | {b['logical_qubits']} | {b['toffoli']:,} | {b['noiseless'][0]:.3f} | {hit['d'] if hit else '—'} | "
                   f"{fmt_q(hit['physical_qubits']['total']) if hit else '—'} | {fmt_time(hit['seconds']) if hit else '—'} |")
    q = [jl(f) for f in real(glob.glob(str(ROOT / "data/results/qpe/*.json")))]
    qb = {r["molecule"]: r for r in q if "d" not in r}
    if qb:
        out.append("\n**Phase estimation** (chemical accuracy, 1.6 mHa from FCI; p = 0.1%, cultivated):\n")
        out.append("| molecule | T gates per run | noiseless | d for 90% of it | run time |\n|---|---:|---:|---:|---:|")
        for m, b in qb.items():
            rr = sorted((r for r in q if "d" in r and r["molecule"] == m and r["factory"] == "cultivation"), key=lambda r: r["d"])
            hit = next((r for r in rr if r["scores"][0] >= 0.9 * b["noiseless"][0]), None)
            out.append(f"| {m} | {sci(b['t_count'])} | {b['noiseless'][0]:.2f} | {hit['d'] if hit else '—'} | {fmt_time(hit['seconds']) if hit else '—'} |")
    sc = jl(ROOT / "data/results/scaling.json")
    tb = next(r for r in sc["textbook"] if r["n"] == 2048 and r["factory"] == "cultivation")
    g = next(p for p in sc["published"] if p["name"].startswith("RSA"))
    out.append(f"\n**To scale** (≤ 0.1 expected faults, p = 0.1%): textbook-arithmetic RSA-2048 needs d = {tb['d']}, "
               f"{fmt_q(tb['physical_qubits'])} physical qubits and {fmt_time(tb['seconds'])}; Gidney's 2025 counts need "
               f"d = {g['estimate']['d']} under our measured error model (he assumed d = 25).")
    return "\n".join(out)


def main():
    path = ROOT / "README.md"
    text = path.read_text()
    a, b = text.index("<!-- results:start -->"), text.index("<!-- results:end -->")
    new = text[: a + len("<!-- results:start -->")] + "\n" + build() + "\n" + text[b:]
    if "--check" in sys.argv:
        sys.exit(0 if new == text else 1)
    path.write_text(new)


if __name__ == "__main__":
    main()
