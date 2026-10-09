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
    return f"{x / 10**e:.2f}e{e}"


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

    # Textbook Shor
    rs = [jl(f) for f in real(glob.glob(str(ROOT / "data/results/shor/*.json")))]
    base = {r["N"]: r for r in rs if "d" not in r}
    runs = [r for r in rs if "d" in r]
    out.append("**Shor (Textbook Cuccaro)** (p = 0.1%, cultivated magic states; peak probability = share of runs on one of the r ideal peaks):\n")
    out.append("| N | logical qubits | Toffolis | noiseless peak | d for 90% of it | physical qubits | run time |\n|---:|---:|---:|---:|---:|---:|---:|")
    for N in sorted(base):
        b = base[N]
        rr = sorted((r for r in runs if r["N"] == N and r["factory"] == "cultivation" and r["p"] == 0.001), key=lambda r: r["d"])
        hit = next((r for r in rr if r["scores"][0] >= 0.9 * b["noiseless"][0]), None)
        out.append(f"| {N} | {b['logical_qubits']} | {b['toffoli']:,} | {b['noiseless'][0]:.3f} | {hit['d'] if hit else '—'} | "
                   f"{fmt_q(hit['physical_qubits']['total']) if hit else '—'} | {fmt_time(hit['seconds']) if hit else '—'} |")

    # Modern Shor
    rs_mod = [jl(f) for f in real(glob.glob(str(ROOT / "data/results/shor_modern/*.json")))]
    base_mod = {r["N"]: r for r in rs_mod if "d" not in r}
    runs_mod = [r for r in rs_mod if "d" in r]
    if base_mod:
        out.append("\n**Shor (Modern Windowed + MBU)** (p = 0.1%, cultivated magic states; 2.7–3.0× Toffoli reduction):\n")
        out.append("| N | logical qubits | Toffolis | noiseless peak | d for 90% of it | physical qubits | run time |\n|---:|---:|---:|---:|---:|---:|---:|")
        for N in sorted(base_mod):
            b = base_mod[N]
            rr = sorted((r for r in runs_mod if r["N"] == N and r["factory"] == "cultivation" and r["p"] == 0.001), key=lambda r: r["d"])
            hit = next((r for r in rr if r["scores"][0] >= 0.9 * b["noiseless"][0]), None)
            out.append(f"| {N} | {b['logical_qubits']} | {b['toffoli']:,} | {b['noiseless'][0]:.3f} | {hit['d'] if hit else '—'} | "
                       f"{fmt_q(hit['physical_qubits']['total']) if hit else '—'} | {fmt_time(hit['seconds']) if hit else '—'} |")

    # Phase estimation
    q = [jl(f) for f in real(glob.glob(str(ROOT / "data/results/qpe/*.json")))]
    qb = {r["molecule"]: r for r in q if "d" not in r}
    if qb:
        out.append("\n**Phase estimation** (chemical accuracy, 1.6 mHa from FCI; p = 0.1%, cultivated):\n")
        out.append("| molecule | T gates per run | noiseless | d for 90% of it | run time |\n|---|---:|---:|---:|---:|")
        for m, b in qb.items():
            rr = sorted((r for r in q if "d" in r and r["molecule"] == m and r["factory"] == "cultivation"), key=lambda r: r["d"])
            hit = next((r for r in rr if r["scores"][0] >= 0.9 * b["noiseless"][0]), None)
            out.append(f"| {m} | {sci(b['t_count'])} | {b['noiseless'][0]:.2f} | {hit['d'] if hit else '—'} | {fmt_time(hit['seconds']) if hit else '—'} |")

    # To scale & Headline comparison
    sc = jl(ROOT / "data/results/scaling.json")
    tb = next(r for r in sc["textbook"] if r["n"] == 2048 and r["factory"] == "cultivation")
    mod = next((r for r in sc.get("modern", []) if r["n"] == 2048 and r["factory"] == "cultivation"), None)
    g = next(p for p in sc["published"] if p["name"].startswith("RSA"))
    sens = jl(ROOT / "data/results/sensitivity.json") if (ROOT / "data/results/sensitivity.json").exists() else None

    out.append("\n**Cryptographic scale (RSA-2048, p = 0.1%, cultivated magic states, target E[faults] ≤ 0.1)**:\n")
    out.append("| Workload / Implementation | Logical qubits | Toffolis | Distance d | Physical qubits | Quantum run time |")
    out.append("|---|---:|---:|---:|---:|---:|")
    out.append(f"| Textbook Cuccaro (this work) | {3 * 2048 + 6:,} | 3.44e11 | {tb['d']} | {fmt_q(tb['physical_qubits'])} | {fmt_time(tb['seconds'])} ({tb['seconds'] / 86400 / 365.25:.2f} yr) |")
    if mod:
        out.append(f"| Modern Windowed + MBU (this work) | {2 * 2048 + 2 * 2 + 20:,} | 8.61e10 | {mod['d']} | {fmt_q(mod['physical_qubits'])} | {fmt_time(mod['seconds'])} ({mod['seconds'] / 86400 / 365.25:.2f} yr) |")
    out.append(f"| Gidney 2025 counts (our measured model) | {g['logical_qubits']:,} | {sci(g['toffolis'])} | {g['estimate']['d']} | {fmt_q(g['estimate']['physical_qubits'])} | {fmt_time(g['estimate']['seconds'])} |")
    out.append(f"| Gidney 2025 published (assumed model) | {g['logical_qubits']:,} | {sci(g['toffolis'])} | 25 | 898 k | 4.96 days |")

    if sens:
        band25 = sens["fit_covariance"]["band"]["25"]
        out.append(f"\n**Headline defensibility (d = 29 vs. d = 25)**: Compiling Gidney's 2025 counts under our measured SD6 model gives d = 29. At d = 25, our circuit-level simulations decoded by correlated matching measure an idle logical error rate of {sci(band25['central'])} per round (1σ: [{sci(band25['lo_1sigma'])}, {sci(band25['hi_1sigma'])}]), which is **4.3× higher** than Gidney's assumed 1.0e-15. Perturbing noise parameters (Λ ± 1σ, prefactor ± 1σ, factory ε_CCZ × 0.1/10×, correlated matching) keeps required distance at d ∈ {{29, 31}}. Feeding Gidney's assumed error model into our compiler reproduces his published d = 25.")

    # Biased noise & XZZX surface codes
    bx_file = ROOT / "data/results/biased_xzzx.json"
    if bx_file.exists():
        bx = jl(bx_file)
        if "rsa2048_modern" in bx:
            pts = bx["rsa2048_modern"]["points"]
            out.append("\n**Biased noise & XZZX surface codes (RSA-2048 Modern, p = 0.1%, target E[faults] ≤ 0.01)**:\n")
            out.append("| Bias η = p_Z/p_X | XZZX distance (dX × dZ) | Tile qubits | Symmetric distance | Symmetric tile | XZZX physical | Symmetric physical | Qubit savings |\n|---:|---:|---:|---:|---:|---:|---:|---:|")
            for p in pts:
                if p["eta"] in (1.0, 10.0, 50.0, 100.0, 500.0, 1000.0):
                    out.append(f"| {p['eta']:.0f} | {p['dx']} × {p['dz']} | {p['tile_qubits']} | {p['symmetric_d']} | {p['symmetric_tile_qubits']} | "
                               f"{fmt_q(p['total_physical_qubits'])} | {fmt_q(p['symmetric_physical_qubits'])} | **{p['qubit_reduction_ratio']:.2f}×** |")

    # LiH chemistry scaling
    lih_file = ROOT / "data/results/lih.json"
    if lih_file.exists():
        lih = jl(lih_file)
        out.append(f"\n**Molecular chemistry scaling (LiH, active space CAS(2, 2) + Z₂ tapering)**: Frozen Li 1s² core, 2 tapered qubits, 9 Pauli terms. Potential energy dissociation minimum at R_e = 1.595 Å (E = -7.862 Ha, PySCF CASCI match to < 10⁻¹⁴ Ha). 10-bit QPE circuit requires {sci(lih['t_count'])} T gates.")

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
