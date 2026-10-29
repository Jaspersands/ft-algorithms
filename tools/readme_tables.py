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
    """Drops sync-conflict copies ("name 2.json"), in a fixed order: glob's order is the
    filesystem's, and the tables must not depend on it."""
    return sorted(f for f in paths if not re.search(r" \d+\.json$", f))


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
        cut = [base[N]["toffoli"] / base_mod[N]["toffoli"] for N in base_mod if N in base]
        out.append(f"\n**Shor (Modern Windowed + MBU)** (p = 0.1%, cultivated magic states; {min(cut):.1f}–{max(cut):.1f}× fewer Toffolis):\n")
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
    lih = jl(ROOT / "data/results/lih.json")
    if qb:
        out.append("\n**Phase estimation** (p = 0.1%, cultivated magic states; success = energy within chemical accuracy, "
                   "1.6 mHa, of the exact ground state of the Hamiltonian simulated):\n")
        out.append("| molecule | system qubits | T gates per run | noiseless | d for 90% of it | run time |\n|---|---:|---:|---:|---:|---:|")
        for m, b in qb.items():
            rr = sorted((r for r in q if "d" in r and r["molecule"] == m and r["factory"] == "cultivation"), key=lambda r: r["d"])
            hit = next((r for r in rr if r["scores"][0] >= 0.9 * b["noiseless"][0]), None)
            out.append(f"| {m} | {b['qubits']} | {sci(b['t_count'])} | {b['noiseless'][0]:.2f} | {hit['d'] if hit else '—'} | "
                       f"{fmt_time(hit['seconds']) if hit else '—'} |")
        rr = sorted((r for r in lih["runs"] if r["factory"] == "cultivation" and r["p"] == 0.001), key=lambda r: r["d"])
        hit = next((r for r in rr if r["scores"][0] >= 0.9 * lih["noiseless"][0]), None)
        out.append(f"| LiH, CAS(2, 2) | {lih['qubits']} | {sci(lih['t_count'])} | {lih['noiseless'][0]:.2f} | {hit['d'] if hit else '—'} | "
                   f"{fmt_time(hit['seconds']) if hit else '—'} |")
        out.append(f"\nH₂ and HeH⁺ in STO-3G are simulated whole (tapered to {qb['H2']['qubits']} qubits). LiH is a deliberately small "
                   f"active space: Li 1s² frozen, two electrons in two CASSCF orbitals, {lih['terms']} Pauli terms on {lih['qubits']} "
                   f"qubits; its exact energy is {abs(lih['e_active'] - lih['e_fci_full']) * 1e3:.1f} mHa above full STO-3G FCI at "
                   f"R = {lih['r_equilibrium']} Å.")

    # To scale
    sc = jl(ROOT / "data/results/scaling.json")
    out.append(f"\n**RSA-2048, this work's circuits** (p = 0.1%, cultivated magic states, E[faults] ≤ {sc['target_expected_faults']:g}; "
               "operation counts exact, schedule fitted on compiled instances, error rates extrapolated):\n")
    cal = [("separate measure and reset layers", sc)] + ([("one measure-and-reset layer", sc["mr"])] if sc.get("mr") else [])
    head = "| arithmetic | calibration | logical qubits | Toffolis | d | physical qubits | run time |\n|---|---|---:|---:|---:|---:|---:|"
    out.append(head if len(cal) > 1 else head.replace(" calibration |", "").replace("|---|---|", "|---|"))
    for lab, k in (("Textbook Cuccaro", "textbook"), ("Modern windowed + MBU", "modern")):
        for clab, src in cal:
            r = next(x for x in src[k] if x["n"] == 2048 and x["factory"] == "cultivation")
            c = f" {clab} |" if len(cal) > 1 else ""
            out.append(f"| {lab} |{c} {r['logical_qubits']:,} | {sci(r['toffolis'])} | {r['d']} | {fmt_q(r['physical_qubits'])} | {fmt_time(r['seconds'])} |")

    # Published estimates, each at its authors' budget
    sens = jl(ROOT / "data/results/sensitivity.json")
    out.append("\n**Published estimates under the measured model** (p = 0.1%; each at its authors' failure budget; distances "
               "only, since their qubit totals use layouts this model does not include):\n")
    mrp = {g["name"]: g for g in (sc.get("mr") or {}).get("published", [])}
    if mrp:
        out.append("| workload | failure budget | published d | d, measured model | d, recalibrated with one measure-and-reset layer |\n|---|---:|---:|---:|---:|")
    else:
        out.append("| workload | failure budget | published d | d, measured model |\n|---|---:|---:|---:|")
    for g in sc["published"]:
        extra = f" {mrp[g['name']]['estimate']['d']} |" if mrp else ""
        out.append(f"| {g['name']} | {(1 - math.exp(-g['target'])) * 100:.1f}% | {g['published_d']} | {g['estimate']['d']} |{extra}")
    pr = sens["per_round"]
    rows25 = {r["d"]: r for r in pr["rows"]}
    meas = [r["measured"] / r["gidney"] for r in pr["rows"] if not r["extrapolated"]]
    lam = sorted(f["lambda"] for f in sens["fits"].values())
    r25 = rows25[sens["published_d"]]
    budget90 = next(s for s in sens["scenarios"] if s["category"] == "Budget" and "90%" in s["name"])
    out.append(
        f"\nGidney 2025 picks d = {sens['published_d']} because his projection of surface-code memory reaches "
        f"{pr['target']:.0e} errors per patch per round there. The measured model (SD6 noise, correlated matching) "
        f"reaches it at d = {pr['d_measured']} (d = {pr['d_band']['-1']}–{pr['d_band']['1']} within the fit's 1σ band). "
        f"At the measured distances its rates are already {min(meas):.1f}–{max(meas):.1f}× his projection, and the gap "
        f"widens with d because the fitted suppression per two units of distance (Λ = {lam[0]:.1f} and {lam[-1]:.1f} for "
        f"the two bases) is below the projection's {sens['gidney_lambda']:.2f}; at d = {sens['published_d']} the fit extrapolates to "
        f"{sci(r25['measured'])} per patch-round (1σ: {sci(r25['lo1'])}–{sci(r25['hi1'])}). One shot of his algorithm "
        f"through the estimator needs d = {sens['baseline']['d']} at his budget (d = {budget90['d']} at 90%); the same "
        f"estimator with his projection in place of the measured model returns d = {sens['gidney_projection']['d']}. "
        f"The table below shows how much of that gap is the circuit and noise model rather than his design.")

    nm = sens.get("noise_models")
    if nm:
        lab = {"sd6": "SD6, as calibrated: separate measure and reset layers", "sd6_mr": "SD6, one measure-and-reset layer",
               "si1000": "SI1000 on CZ, one measure-and-reset layer", "uniform": "uniform depolarizing, Stim's generated circuit"}
        out.append("\n**Which memory model?** Memories at p = 0.1% under other circuits and noise models, each decoded by "
                   "correlated matching and fitted the same way (per-round rates from d and 3d rounds), against the projection:\n")
        out.append("| memory model | Λ (Z, X) | measured to d | rate ÷ projection there | d, worse basis ≤ 1e-15 | d, either basis ≤ 1e-15 | d, one shot at his budget |")
        out.append("|---|---:|---:|---:|---:|---:|---:|")
        for k in ("sd6", "sd6_mr", "si1000", "uniform"):
            if k not in nm:
                continue
            m = nm[k]
            top = m["measured"][-1]
            out.append(f"| {lab[k]} | {m['lambda']['z']:.1f}, {m['lambda']['x']:.1f} | {top['d']} | {(top['z'] + top['x']) / top['gidney']:.1f}× | "
                       f"{m['d_basis'] or '> 61'} | {m['d_any'] or '> 61'} | {m['d_shot'] or '> 61'} |")
        out.append(f"| Gidney 2025 projection, 3.5^−d/40 | {sens['gidney_lambda']:.2f} | — | 1× | {pr['d_gidney']} | {pr['d_gidney']} | "
                   f"{sens['gidney_projection']['d']} |")
        ds = [nm[k]["d_any"] for k in nm if nm[k]["d_any"]]
        top = {k: nm[k]["measured"][-1] for k in nm}
        mr = {b: top["sd6"][b] / top["sd6_mr"][b] for b in ("z", "x")}
        si = nm["si1000"]
        out.append(
            f"\nThe distance that reaches 1e-15 per patch-round at p = 0.1% runs from d = {min(ds)} to d = {max(ds)} across these "
            f"circuit-level models, and Gidney's d = {sens['published_d']} lies inside that range. SI1000, the noise model of the "
            f"yoked-surface-code paper's simulations (arXiv:2312.04522), has the projection's slope here (Λ = "
            f"{min(si['lambda'].values()):.1f}–{max(si['lambda'].values()):.1f}) but a "
            f"{(top['si1000']['z'] + top['si1000']['x']) / top['si1000']['gidney']:.1f}× higher rate at d = {top['si1000']['d']}, and "
            f"needs d = {si['d_any']}. Gidney 2025 labels the projection uniform depolarizing; Stim's generated circuit with every "
            f"parameter at p, one reading of that which leaves qubits idle during gate layers noiseless, is below it. The "
            f"calibration's circuit is the most pessimistic of the four: it "
            f"depolarizes idle qubits in every layer and measures and resets the ancillas in separate layers, so data qubits idle "
            f"through two of them each round. Merging them alone lowers the d = {top['sd6']['d']} rates "
            f"{min(mr.values()):.1f}–{max(mr.values()):.1f}× and raises Λ to ≈ {min(nm['sd6_mr']['lambda'].values()):.0f}, so every "
            f"result here leans pessimistic by about that. (The SD6 row is the calibration's own memories fitted like the others, "
            f"over d ≤ {top['sd6']['d']}; with all its memory lengths and d = 11 the calibrated model needs d = {pr['d_measured']}.)")
    mr = sens.get("mr")
    if mr:
        lam_mr = sorted(f["lambda"] for f in mr["fits"].values())
        out.append(
            f"\n**Recalibrated with one measure-and-reset layer.** Every calibration experiment was rerun on the same "
            f"circuits with each reset layer folded into the measurement layer before it (`tools/calibrate.py --mr`; "
            f"data/calibration/raw_mr/, model_mr.json). Idling then has Λ = {lam_mr[0]:.1f}–{lam_mr[-1]:.1f}, reaches 1e-15 per "
            f"patch-round at d = {mr['d_round']} (d = {mr['d_band']['-1']}–{mr['d_band']['1']} within 1σ), and one shot of "
            f"Gidney's algorithm needs d = {mr['one_shot']['d']} at his budget. The tables above show both calibrations; the "
            f"simulated algorithm runs (Shor, phase estimation) use the original one.")

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
