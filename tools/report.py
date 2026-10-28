"""Builds the technical report from the committed data: report/report.md, report/report.html,
report/report.pdf and the figures in report/figures/. Every number is computed here from
data/; none is typed by hand.

    .venv/bin/python tools/report.py            # everything
    .venv/bin/python tools/report.py --no-pdf
"""

from __future__ import annotations

import argparse
import datetime
import glob
import json
import math
import os
import pathlib
import re
import subprocess
import sys
import tempfile
import shutil

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "python"))

import matplotlib  # noqa: E402

matplotlib.use("svg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from ftalgo.arch import FACTORIES  # noqa: E402
from ftalgo.calib.model import LogicalModel  # noqa: E402

REPORT = ROOT / "report"
FIGS = REPORT / "figures"
CHROME = os.environ.get("CHROME", "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome")
C = dict(s1="#2f5fd0", s2="#c0392b", s3="#2c7a4b", s4="#7b4bb5", s5="#a8620b", s6="#1b8a99", s7="#4a4c50",
         ink="#1c1d1f", ink2="#4a4c50", ink3="#6b6d71", rule="#ded9cd", paper="#fffefb")
NCOL = [C["s1"], C["s2"], C["s3"], C["s4"], C["s5"], C["s6"], C["s7"]]
FCOL = {"cultivation": C["s3"], "15to1": C["s1"], "injected": C["s2"]}
FLAB = {"cultivation": "cultivated", "15to1": "15-to-1 distilled", "injected": "injected"}

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 8.5, "axes.edgecolor": C["ink"], "axes.linewidth": 0.8,
    "axes.labelcolor": C["ink2"], "xtick.color": C["ink3"], "ytick.color": C["ink3"], "axes.grid": True,
    "grid.color": C["rule"], "grid.linewidth": 0.6, "legend.frameon": False, "legend.fontsize": 7.5,
    "svg.hashsalt": "ftalgo", "figure.dpi": 100, "axes.spines.top": False, "axes.spines.right": False,
})


def jload(p):
    return json.loads(pathlib.Path(p).read_text())


def runs_of(kind):
    rs = [jload(f) for f in sorted(glob.glob(str(ROOT / "data" / "results" / kind / "*.json"))) if not re.search(r" \d+\.json$", f)]
    return [r for r in rs if "d" in r], [r for r in rs if "d" not in r]


def sci(x, digits=2):
    if x == 0:
        return "0"
    e = math.floor(math.log10(abs(x)))
    m = x / 10**e
    if abs(round(m, digits - 1)) >= 10:  # 9.95 rounds to 10.0: carry into the exponent
        e, m = e + 1, m / 10
    if -2 <= e <= 3:
        return f"{x:.{max(digits - 1 - e, 0)}f}" if e < 0 else f"{x:,.{max(digits - 1 - e, 0)}f}"
    return f"{m:.{digits - 1}f}×10^{e}^"


def fmt_time(s):
    if s < 1:
        return f"{s * 1e3:.3g} ms"
    if s < 120:
        return f"{s:.3g} s"
    if s < 7200:
        return f"{s / 60:.3g} min"
    if s < 2 * 86400:
        return f"{s / 3600:.3g} h"
    if s < 400 * 86400:
        return f"{s / 86400:.3g} days"
    return f"{s / 86400 / 365.25:.3g} years"


def fmt_q(x):
    return f"{x / 1e6:.3g} M" if x >= 1e6 else f"{x / 1e3:.3g} k" if x >= 1e3 else f"{x:.0f}"


def table(head, rows, align=None):
    align = align or ["l"] + ["r"] * (len(head) - 1)
    sep = ["---:" if a == "r" else ":---" for a in align]
    out = ["| " + " | ".join(head) + " |", "| " + " | ".join(sep) + " |"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def savefig(fig, name):
    FIGS.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGS / f"{name}.svg", bbox_inches="tight", metadata={"Date": None})
    plt.close(fig)
    return f"figures/{name}.svg"


def fitv(f, d):
    return math.exp(f.a + f.b * (d + 1) / 2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-pdf", action="store_true")
    args = ap.parse_args()
    model = LogicalModel.load()
    raw = [jload(f) for f in sorted(glob.glob(str(ROOT / "data" / "calibration" / "raw" / "*.json"))) if not re.search(r" \d+\.json$", f)]
    comp = jload(ROOT / "data" / "calibration" / "composition.json")
    xc = jload(ROOT / "data" / "calibration" / "xcheck.json")
    shor_runs, shor_base = runs_of("shor")
    shor_mod_runs, shor_mod_base = runs_of("shor_modern")
    qpe_runs, qpe_base = runs_of("qpe")
    sc = jload(ROOT / "data" / "results" / "scaling.json")
    sens = jload(ROOT / "data" / "results" / "sensitivity.json")
    lih_path = ROOT / "data" / "results" / "lih.json"
    lih = jload(lih_path) if lih_path.exists() else None
    V = {}

    # -- calibration ---------------------------------------------------------------------------
    V["n_exp"] = len(raw)
    V["shots"] = sum(r["counts"]["shots"] for r in raw)
    V["fails"] = sum(sum(r["counts"]["patterns"].values()) for r in raw)
    lam_rows = []
    for p in (0.001, 0.002, 0.003, 0.005):
        def lam(name):
            f = model.per_p.get(name, {}).get(p)
            return f"{math.exp(-f.b):.2f}" if f else "—"
        lam_rows.append([f"{p * 100:.1f}%", lam("idle_xy"), lam("idle_zy"), lam("cnot_x2"), lam("cnot_z1"), lam("zz_o"), lam("zz_j1")])
    V["lam_table"] = table(["p", "idle (X∪Y)", "idle (Z∪Y)", "CNOT X on target", "CNOT Z on control", "Z⊗Z patch Z", "Z⊗Z wrong outcome"], lam_rows)
    f_idle = model.per_p["idle_xy"][0.001]
    V["lam_idle"] = math.exp(-f_idle.b)
    V["idle_d9"] = fitv(f_idle, 9)
    V["idle_d25"] = fitv(f_idle, 25)
    yp = {q.d: q.value for q in model.points["idle_y"] if q.p == 0.001}
    xp = {q.d: q.value for q in model.points["idle_xy"] if q.p == 0.001}
    zp = {q.d: q.value for q in model.points["idle_zy"] if q.p == 0.001}
    V["y_frac_d3"] = yp[3] / xp[3]
    V["y_vs_indep_d3"] = yp[3] / (xp[3] * zp[3])
    ratios = sorted(e["any"]["predicted"] / e["any"]["measured"] for r in comp["results"] for e in r["experiments"].values())
    V["comp_lo"], V["comp_hi"], V["comp_med"] = ratios[0], ratios[-1], ratios[len(ratios) // 2]
    V["comp_n"] = len(ratios)
    zs = [abs(c["z"]) for c in xc["cases"]]
    V["xc_n"], V["xc_maxz"] = len(zs), max(zs)

    # Figure: idle per round.
    fig, ax = plt.subplots(figsize=(5.6, 3.2))
    for p, col in zip((0.001, 0.002, 0.003, 0.005), NCOL):
        pts = [q for q in model.points["idle_xy"] if q.p == p and q.events >= 10 and q.value > 0]
        ax.errorbar([q.d for q in pts], [q.value for q in pts], yerr=[q.sigma for q in pts], fmt="o", ms=3.5, color=col, label=f"p = {p * 100:.1f}%", capsize=0)
        f = model.per_p["idle_xy"].get(p)
        if f:
            dd = np.arange(3, f.dmax + 0.01, 0.25)
            de = np.arange(f.dmax, 15.01, 0.25)
            ax.plot(dd, [fitv(f, d) for d in dd], color=col, lw=1.1)
            ax.plot(de, [fitv(f, d) for d in de], color=col, lw=1.1, ls="--")
    ax.set_yscale("log")
    ax.set_xlabel("code distance d")
    ax.set_ylabel("logical error per round")
    ax.set_xticks(range(3, 16, 2))
    ax.legend(ncol=4, loc="lower left")
    V["fig_idle"] = savefig(fig, "idle")

    # Figure: composition.
    fig, ax = plt.subplots(figsize=(4.2, 3.8))
    kinds = {"cnot-z": ("CNOT, Z inputs", C["s1"]), "cnot-x": ("CNOT, X inputs", C["s6"]), "repeated_zz_k3": ("three Z⊗Z", C["s2"]), "line_n3": ("three-patch merge", C["s3"])}
    for k, (lab, col) in kinds.items():
        xs, ys = [], []
        for r in comp["results"]:
            e = r["experiments"].get(k, {})
            for name, v in e.items():
                if name != "any" and v["measured"] > 0:
                    xs.append(v["measured"])
                    ys.append(v["predicted"])
        ax.scatter(xs, ys, s=9, color=col, label=lab, zorder=3)
    g = np.array([1e-4, 1])
    ax.plot(g, g, color=C["ink"], lw=0.8)
    ax.plot(g, g * 1.25, color=C["ink3"], lw=0.7, ls="--")
    ax.plot(g, g / 1.25, color=C["ink3"], lw=0.7, ls="--")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1e-4, 1)
    ax.set_ylim(1e-4, 1)
    ax.set_xlabel("measured (circuit level)")
    ax.set_ylabel("predicted from parts (logical level)")
    ax.legend(loc="upper left")
    V["fig_comp"] = savefig(fig, "composition")

    # -- Shor Textbook -------------------------------------------------------------------------
    base = {b["N"]: b for b in shor_base}
    Ns = sorted(base)
    sel = lambda N, f, p: sorted((r for r in shor_runs if r["N"] == N and r["factory"] == f and r["p"] == p), key=lambda r: r["d"])  # noqa: E731
    rows = []
    V["shor"] = {}
    for N in Ns:
        b = base[N]
        rs = sel(N, "cultivation", 0.001)
        hit = next((r for r in rs if r["scores"][0] >= 0.9 * b["noiseless"][0]), None)
        V["shor"][N] = hit
        rows.append([N, b["n"], b["logical_qubits"], f"{b['toffoli']:,}", f"{b['t']:,}", f"{b['noiseless'][0]:.3f}", f"{b['random'][0]:.4f}", f"{b['random'][1]:.3f}",
                     hit["d"] if hit else "—", fmt_q(hit["physical_qubits"]["total"]) if hit else "—", fmt_time(hit["seconds"]) if hit else "—"])
    V["shor_table"] = table(["N", "n", "logical qubits", "Toffolis", "T", "noiseless peak", "random peak", "random factors", "d for 90%", "physical qubits", "run time"], rows)
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.9))
    for i, N in enumerate(Ns):
        rs = sel(N, "cultivation", 0.001)
        if not rs:
            continue
        col = NCOL[i % len(NCOL)]
        for k, ax in enumerate(axs):
            d = [r["d"] for r in rs]
            y = [r["scores"][k] for r in rs]
            e = [2 * r["score_sigmas"][k] for r in rs]
            ax.errorbar(d, y, yerr=e, color=col, ms=3.5, fmt="-o", lw=1.1, label=f"N = {N}")
            ax.axhline(base[N]["noiseless"][0] if k == 0 else base[N]["random"][1], color=col, lw=0.7, ls=":")
    axs[0].set_ylabel("peak probability")
    axs[1].set_ylabel("factors found")
    for ax in axs:
        ax.set_xlabel("code distance d")
        ax.set_ylim(-0.02, 1.03)
    axs[0].legend(loc="lower right", ncol=2)
    V["fig_shor"] = savefig(fig, "shor")

    fN = next((N for N in (35, 21, 15) if all(sel(N, f, 0.001) for f in FACTORIES)), Ns[0])
    V["fN"] = fN
    fig, ax = plt.subplots(figsize=(5.2, 3.0))
    frows = []
    for f in ("cultivation", "15to1", "injected"):
        for p in (0.001, 0.002):
            rs = sel(fN, f, p)
            if not rs:
                continue
            ax.errorbar([r["d"] for r in rs], [r["scores"][0] for r in rs], yerr=[2 * r["score_sigmas"][0] for r in rs], color=FCOL[f], fmt="o", ms=3.2, lw=1.1, ls="--" if p == 0.002 else "-", label=f"{FLAB[f]}, p = {p * 100:.1f}%")
            best = max(rs, key=lambda r: r["scores"][0])
            hit = next((r for r in rs if r["scores"][0] >= 0.9 * base[fN]["noiseless"][0]), None)
            frows.append([FLAB[f], f"{p * 100:.1f}%", sci(FACTORIES[f].eps_t(p)), sci(FACTORIES[f].eps_ccz(p)), f"{best['scores'][0]:.3f} (d = {best['d']})", hit["d"] if hit else "never"])
    ax.set_xlabel("code distance d")
    ax.set_ylabel(f"peak probability, N = {fN}")
    ax.set_ylim(-0.02, 1.03)
    ax.legend(fontsize=6.8, loc="center right")
    V["fig_fact"] = savefig(fig, "factories")
    V["fact_table"] = table(["magic states", "p", "ε_T †", "ε_CCZ †", "best peak", "d for 90%"], frows)

    hb = V["shor"].get(fN)
    if hb:
        tot = sum(hb["budget"].values())
        V["budget_rows"] = table(["class", "expected faults", "share", "score of one-fault runs"],
                                 [[k, sci(v), f"{100 * v / tot:.1f}%", f"{hb['error_budget'][k]['single_fault_score']:.3f}" if hb["error_budget"].get(k, {}).get("single_fault_score") is not None else "—"]
                                  for k, v in sorted(hb["budget"].items(), key=lambda kv: -kv[1])])
        V["budget_d"] = hb["d"]
        V["idle_share"] = hb["budget"].get("idle", 0) / tot

    # -- Shor Modern ---------------------------------------------------------------------------
    base_mod = {b["N"]: b for b in shor_mod_base}
    sel_mod = lambda N, f, p: sorted((r for r in shor_mod_runs if r["N"] == N and r["factory"] == f and r["p"] == p), key=lambda r: r["d"])  # noqa: E731
    mod_rows = []
    for N in (15, 21, 35):
        if N not in base or N not in base_mod:
            continue
        tb_b = base[N]
        mod_b = base_mod[N]
        tb_rs = sel(N, "cultivation", 0.001)
        mod_rs = sel_mod(N, "cultivation", 0.001)
        tb_hit = next((r for r in tb_rs if r["scores"][0] >= 0.9 * tb_b["noiseless"][0]), None)
        mod_hit = next((r for r in mod_rs if r["scores"][0] >= 0.9 * mod_b["noiseless"][0]), None)
        tb_d11 = next((r for r in tb_rs if r["d"] == 11), None)
        mod_d11 = next((r for r in mod_rs if r["d"] == 11), None)
        mod_rows.append([N, "textbook (Cuccaro)", tb_b["logical_qubits"], f"{tb_b['toffoli']:,}", f"{tb_d11['rounds']:,}" if tb_d11 else "—",
                         f"{tb_d11['scores'][0]:.3f}" if tb_d11 else "—", tb_hit["d"] if tb_hit else "—", fmt_q(tb_hit["physical_qubits"]["total"]) if tb_hit else "—", fmt_time(tb_hit["seconds"]) if tb_hit else "—"])
        mod_rows.append([N, "windowed (k = 2, MBU)", mod_b["logical_qubits"], f"{mod_b['toffoli']:,}", f"{mod_d11['rounds']:,}" if mod_d11 else "—",
                         f"{mod_d11['scores'][0]:.3f}" if mod_d11 else "—", mod_hit["d"] if mod_hit else "—", fmt_q(mod_hit["physical_qubits"]["total"]) if mod_hit else "—", fmt_time(mod_hit["seconds"]) if mod_hit else "—"])
    V["shor_mod_table"] = table(["N", "arithmetic", "logical qubits", "Toffolis", "rounds (d=11)", "peak (d=11)", "d for 90%", "physical qubits", "run time"], mod_rows)
    mod_Ns_ok = [N for N in (15, 21, 35) if N in base and N in base_mod]
    cut = [base[N]["toffoli"] / base_mod[N]["toffoli"] for N in mod_Ns_ok]
    d11 = {N: (next(r for r in sel(N, "cultivation", 0.001) if r["d"] == 11), next(r for r in sel_mod(N, "cultivation", 0.001) if r["d"] == 11)) for N in mod_Ns_ok}
    rcut = [t["rounds"] / m["rounds"] for t, m in d11.values()]
    hit90 = {}
    for N in mod_Ns_ok:
        hit90[N] = tuple(next((r["d"] for r in rs if r["scores"][0] >= 0.9 * b["noiseless"][0]), None)
                         for rs, b in ((sel(N, "cultivation", 0.001), base[N]), (sel_mod(N, "cultivation", 0.001), base_mod[N])))
    V["mod"] = {"cut_lo": min(cut), "cut_hi": max(cut), "rcut_lo": min(rcut), "rcut_hi": max(rcut), "d11": d11, "hit90": hit90,
                "extra_qubits": {N: base_mod[N]["logical_qubits"] - base[N]["logical_qubits"] for N in mod_Ns_ok}}

    # Figure: Shor modern vs textbook
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.9))
    mod_Ns = [15, 21, 35]
    for i, N in enumerate(mod_Ns):
        tb_rs = sel(N, "cultivation", 0.001)
        mod_rs = sel_mod(N, "cultivation", 0.001)
        col = NCOL[i % len(NCOL)]
        if tb_rs:
            axs[0].plot([r["d"] for r in tb_rs], [r["scores"][0] for r in tb_rs], "--o", ms=3.2, color=col, alpha=0.55, label=f"N = {N}, textbook")
            axs[1].plot([r["d"] for r in tb_rs], [r["rounds"] for r in tb_rs], "--o", ms=3.2, color=col, alpha=0.55)
        if mod_rs:
            axs[0].plot([r["d"] for r in mod_rs], [r["scores"][0] for r in mod_rs], "-s", ms=3.5, color=col, lw=1.2, label=f"N = {N}, windowed")
            axs[1].plot([r["d"] for r in mod_rs], [r["rounds"] for r in mod_rs], "-s", ms=3.5, color=col, lw=1.2)
        if N in base_mod:
            axs[0].axhline(base_mod[N]["noiseless"][0], color=col, lw=0.6, ls=":")
    axs[0].set_xlabel("code distance d")
    axs[0].set_ylabel("peak probability")
    axs[0].set_ylim(-0.02, 1.03)
    axs[0].legend(loc="lower right", fontsize=6.5, ncol=2)
    axs[1].set_xlabel("code distance d")
    axs[1].set_ylabel("circuit rounds")
    axs[1].set_yscale("log")
    V["fig_shor_modern"] = savefig(fig, "shor_modern")

    # -- QPE ----------------------------------------------------------------------------------
    qb = {b["molecule"]: b for b in qpe_base}
    qrows = []
    fig, ax = plt.subplots(figsize=(5.2, 3.0))
    mcol = {"H2": C["s1"], "HeH+": C["s2"]}
    for m, b in qb.items():
        for f in ("cultivation", "15to1"):
            rs = sorted((r for r in qpe_runs if r["molecule"] == m and r["factory"] == f and r["p"] == 0.001), key=lambda r: r["d"])
            if not rs:
                continue
            ax.errorbar([r["d"] for r in rs], [r["scores"][0] for r in rs], yerr=[2 * r["score_sigmas"][0] for r in rs], color=mcol.get(m, C["s3"]), fmt="o", ms=3.2, lw=1.1, ls="-" if f == "cultivation" else "--", label=f"{m}, {FLAB[f]}")
        ax.axhline(b["noiseless"][0], color=mcol.get(m, C["s3"]), lw=0.7, ls=":")
        rs = sorted((r for r in qpe_runs if r["molecule"] == m and r["factory"] == "cultivation" and r["p"] == 0.001), key=lambda r: r["d"])
        hit = next((r for r in rs if r["scores"][0] >= 0.9 * b["noiseless"][0]), None)
        best = max(rs, key=lambda r: r["scores"][0]) if rs else None
        qrows.append([m, b["qubits"] + 1, b["terms"], sci(b["t_count"]), f"{b['e_hf']:.6f}", f"{b['e_fci']:.6f}", f"{b['noiseless'][0]:.2f}",
                      f"{best['scores'][0]:.2f} (d = {best['d']})" if best else "—", hit["d"] if hit else "—", fmt_time(hit["seconds"]) if hit else "—", fmt_q(hit["physical_qubits"]["total"]) if hit else "—"])
    ax.set_xlabel("code distance d")
    ax.set_ylabel("P(chemical accuracy)")
    ax.set_ylim(-0.02, 1.03)
    if ax.get_legend_handles_labels()[0]:
        ax.legend(fontsize=7)
    V["fig_qpe"] = savefig(fig, "qpe")
    V["qpe_table"] = table(["molecule", "logical qubits", "terms", "T gates", "E_HF (Ha)", "E_FCI (Ha)", "noiseless", "best", "d for 90%", "run time", "physical qubits"], qrows)
    V["h2"] = qb.get("H2")

    # Figure: the energies individual runs return (H2, cultivated, p = 0.1%).
    qh_path = ROOT / "data" / "results" / "qpe_hist.json"
    V["fig_qpe_hist"] = None
    if qh_path.exists():
        qh = jload(qh_path)
        hrows = [(f"d = {r['d']}", r["energies"]) for r in qh["runs"]] + [("noiseless", qh["noiseless"])]
        lim, ca, bw = 12.0, 1.6, qh["resolution"] * 1000
        edges = np.arange(-lim, lim + bw / 2, bw)
        fig, axs = plt.subplots(len(hrows), 1, figsize=(5.6, 0.62 * len(hrows) + 0.5), sharex=True)
        for ax, (lab, es) in zip(axs, hrows):
            err = (np.asarray(es) - qh["e_fci"]) * 1000
            ok = float(np.mean(np.abs(err) < ca))
            off = int(np.sum(np.abs(err) >= lim))
            ax.axvspan(-ca, ca, color="#dff0e5", lw=0, zorder=0)
            ax.hist(err[np.abs(err) < lim], bins=edges, color=C["ink3"] if lab == "noiseless" else C["s1"], zorder=2)
            ax.set_yticks([])
            ax.grid(False)
            ax.spines["left"].set_visible(False)
            ax.set_ylabel(lab, rotation=0, ha="right", va="center", fontsize=7.5)
            ax.text(1.01, 0.5, f"{100 * ok:.0f}% ok" + (f"\n{off} off-scale" if off else ""), transform=ax.transAxes, fontsize=6.8, color=C["ink3"], va="center")
        axs[-1].set_xlim(-lim, lim)
        axs[-1].set_xlabel("E − E_FCI (mHa)")
        V["fig_qpe_hist"] = savefig(fig, "qpe_hist")
        V["qpe_hist_n"] = len(qh["noiseless"])

    # -- scaling -------------------------------------------------------------------------------
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.9))
    for f in ("cultivation", "15to1"):
        rr = [r for r in sc["textbook"] if r["factory"] == f and r["d"]]
        rr_mod = [r for r in sc.get("modern", []) if r["factory"] == f and r["d"]]
        axs[0].plot([r["n"] for r in rr], [r["physical_qubits"] for r in rr], "-o", ms=2.6, color=FCOL[f], alpha=0.5, label=f"textbook, {FLAB[f]}")
        axs[0].plot([r["n"] for r in rr_mod], [r["physical_qubits"] for r in rr_mod], "--s", ms=2.6, color=FCOL[f], label=f"windowed, {FLAB[f]}")
        axs[1].plot([r["n"] for r in rr], [r["seconds"] for r in rr], "-o", ms=2.6, color=FCOL[f], alpha=0.5)
        axs[1].plot([r["n"] for r in rr_mod], [r["seconds"] for r in rr_mod], "--s", ms=2.6, color=FCOL[f])
    axs[0].plot([2048], [897864], "D", color=C["s4"], ms=5, label="Gidney 2025 (published, yoked storage)")
    axs[1].plot([2048], [4.96 * 86400], "D", color=C["s4"], ms=5)
    for ax in axs:
        ax.set_xscale("log", base=2)
        ax.set_yscale("log")
        ax.set_xlabel("modulus size n (bits)")
    axs[0].set_ylabel("physical qubits")
    axs[1].set_ylabel("run time (s)")
    axs[0].legend(fontsize=6.5)
    V["fig_scale"] = savefig(fig, "scaling")

    tb = next(r for r in sc["textbook"] if r["n"] == 2048 and r["factory"] == "cultivation")
    mod2048 = next(r for r in sc["modern"] if r["n"] == 2048 and r["factory"] == "cultivation")
    V["scale_table"] = table(["arithmetic", "logical qubits", "Toffolis", "d", "physical qubits", "run time"],
                             [[lab, f"{r['logical_qubits']:,}", sci(r["toffolis"]), r["d"], fmt_q(r["physical_qubits"]), fmt_time(r["seconds"])]
                              for lab, r in (("textbook", tb), ("windowed + MBU", mod2048))])
    V["pub_table"] = table(["workload", "logical qubits", "Toffolis", "run time", "failure budget", "published d", "d, measured model"],
                           [[p["name"], f"{p['logical_qubits']:,}", sci(p["toffolis"]), fmt_time(p["estimate"]["seconds"]),
                             f"{100 * (1 - math.exp(-p['target'])):.1f}%", p["published_d"], p["estimate"]["d"]] for p in sc["published"]])
    V["tb2048"] = tb
    V["mod2048"] = mod2048
    V["femoco"] = next(p for p in sc["published"] if p["name"].startswith("FeMoco"))
    st = sc["structure"]
    st_mod = sc.get("structure_modern", st)

    # -- Gidney 2025's distance ----------------------------------------------------------------
    pr = sens["per_round"]
    meas = [r for r in pr["rows"] if not r["extrapolated"]]
    V["g"] = {"ratio_lo": min(r["measured"] / r["gidney"] for r in meas), "ratio_hi": max(r["measured"] / r["gidney"] for r in meas),
              "dmax": max(r["d"] for r in meas), "lams": sorted(f["lambda"] for f in sens["fits"].values()),
              "r25": next(r for r in pr["rows"] if r["d"] == sens["published_d"]),
              "b90": next(x for x in sens["scenarios"] if x["category"] == "Budget" and "90%" in x["name"])}
    V["per_round_table"] = table(["d", "measured model", "1σ band", "2σ band", "Gidney's projection", "ratio"],
                                 [[f"{r['d']}{'' if not r['extrapolated'] else ' *'}", sci(r["measured"]), f"{sci(r['lo1'])} – {sci(r['hi1'])}",
                                   f"{sci(r['lo2'])} – {sci(r['hi2'])}", sci(r["gidney"]), f"{r['measured'] / r['gidney']:.1f}×"]
                                  for r in pr["rows"] if r["d"] in (3, 7, 11, 21, 23, 25, 27, 29, 31)])
    srows = [["baseline", f"one shot at his budget (E[faults] ≤ {sens['target']:.3f})", sens["baseline"]["d"], "0"]]
    srows += [[x["category"], x["name"], x["d"] if x["d"] is not None else "> 61", f"{x['delta_d']:+d}" if x["delta_d"] is not None else "—"] for x in sens["scenarios"]]
    srows.append(["his model", "Gidney's projection in place of the measured model", sens["gidney_projection"]["d"], f"{sens['gidney_projection']['d'] - sens['baseline']['d']:+d}"])
    V["scen_table"] = table(["varied", "scenario", "d", "Δd"], srows, ["l", "l", "r", "r"])

    fig, ax = plt.subplots(figsize=(5.6, 3.3))
    rows_ = [r for r in pr["rows"] if r["d"] <= 31]
    dd = [r["d"] for r in rows_]
    dm = V["g"]["dmax"]
    ax.fill_between([r["d"] for r in rows_ if r["d"] >= dm], [r["lo1"] for r in rows_ if r["d"] >= dm], [r["hi1"] for r in rows_ if r["d"] >= dm], color=C["s1"], alpha=0.15, lw=0)
    ax.plot([d for d in dd if d <= dm], [r["measured"] for r in rows_ if r["d"] <= dm], "-", color=C["s1"], lw=1.2)
    ax.plot([d for d in dd if d >= dm], [r["measured"] for r in rows_ if r["d"] >= dm], "--", color=C["s1"], lw=1.2, label="measured model (SD6, correlated matching)")
    ax.plot([r["d"] for r in meas], [r["measured"] for r in meas], "o", color=C["s1"], ms=3.5)
    ax.plot(dd, [r["gidney"] for r in rows_], "-", color=C["s4"], lw=1.2, label="Gidney 2025 projection, $3.5^{-d}/40$")
    ax.axhline(pr["target"], color=C["ink3"], lw=0.8, ls=":")
    ax.text(3.3, pr["target"] * 2, "10⁻¹⁵", fontsize=7, color=C["ink3"])
    ax.set_yscale("log")
    ax.set_xlabel("code distance d")
    ax.set_ylabel("logical error per patch per round")
    ax.set_xticks(range(3, 32, 4))
    ax.legend(fontsize=7, loc="upper right")
    V["fig_gidney"] = savefig(fig, "gidney")

    # Memory models: SD6 (as calibrated), Stim's uniform depolarizing, SI1000, against the projection.
    nm = sens.get("noise_models") or {}
    NM_LAB = {"sd6": "SD6 (this calibration)", "uniform": "uniform depolarizing (Stim)", "si1000": "SI1000 (CZ circuits)"}
    NM_COL = {"sd6": C["s1"], "uniform": C["s3"], "si1000": C["s2"]}
    V["nm"] = nm
    if nm:
        fig, ax = plt.subplots(figsize=(5.6, 3.3))
        for k in ("sd6", "uniform", "si1000"):
            if k not in nm:
                continue
            m = nm[k]
            ax.plot([r["d"] for r in m["measured"]], [r["z"] + r["x"] for r in m["measured"]], "o", ms=3.5, color=NM_COL[k])
            dd = [d for d in range(3, 32, 2)]
            ys = [m["fit"][str(d)]["z"] + m["fit"][str(d)]["x"] for d in dd]
            ax.plot([d for d in dd if d <= m["dmax"]], [y for d, y in zip(dd, ys) if d <= m["dmax"]], "-", color=NM_COL[k], lw=1.2)
            ax.plot([d for d in dd if d >= m["dmax"]], [y for d, y in zip(dd, ys) if d >= m["dmax"]], "--", color=NM_COL[k], lw=1.2, label=NM_LAB[k])
        dd = list(range(3, 32, 2))
        ax.plot(dd, [3.5 ** -d / 40 for d in dd], "-", color=C["s4"], lw=1.2, label="Gidney 2025 projection")
        ax.axhline(pr["target"], color=C["ink3"], lw=0.8, ls=":")
        ax.set_yscale("log")
        ax.set_ylim(1e-20, 1e-1)
        ax.set_xlabel("code distance d")
        ax.set_ylabel("logical error per patch per round")
        ax.set_xticks(range(3, 32, 4))
        ax.legend(fontsize=7, loc="upper right")
        V["fig_nm"] = savefig(fig, "noise_models")
        nrows = []
        for k in ("sd6", "uniform", "si1000"):
            if k not in nm:
                continue
            m = nm[k]
            top = m["measured"][-1]
            nrows.append([NM_LAB[k], f"{m['lambda']['z']:.1f}, {m['lambda']['x']:.1f}", top["d"], f"{(top['z'] + top['x']) / top['gidney']:.1f}×",
                          m["d_basis"] or "> 61", m["d_any"] or "> 61", m["d_shot"] or "> 61"])
        nrows.append(["Gidney 2025 projection", f"{sens['gidney_lambda']:.2f}", "—", "1×", sens["per_round"]["d_gidney"], sens["per_round"]["d_gidney"],
                      sens["gidney_projection"]["d"]])
        V["nm_table"] = table(["memory model (p = 0.1%)", "Λ (Z, X)", "largest d", "rate / projection there", "d: worse basis ≤ 10⁻¹⁵",
                               "d: either basis ≤ 10⁻¹⁵", "d: one shot at his budget"], nrows, ["l", "r", "r", "r", "r", "r", "r"])

    # -- LiH Chemistry -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(5.6, 3.2))
    if lih is not None:
        curve = lih["dissociation_curve"]
        r_pts = [p["r_angstrom"] for p in curve]
        e_hf = [p["e_hf"] for p in curve]
        e_casci = [p["e_casci"] for p in curve]
        e_casscf = [p.get("e_casscf", p["e_casci"]) for p in curve]
        e_fci = [p["e_fci"] for p in curve]
    else:
        from ftalgo.chem.lih import lih_dissociation_curve
        r_sweep = np.linspace(1.0, 3.0, 15)
        c_pts = lih_dissociation_curve(r_sweep)
        r_pts = [round(p.r_angstrom, 4) for p in c_pts]
        e_hf = [p.e_hf for p in c_pts]
        e_casci = [p.e_casci for p in c_pts]
        e_casscf = [p.e_casscf for p in c_pts]
        e_fci = [p.e_fci for p in c_pts]

    ax.plot(r_pts, e_hf, "--", color=C["s2"], lw=1.2, label="RHF (restricted HF)")
    ax.plot(r_pts, e_casci, ":", color=C["ink3"], lw=1.0, label="CASCI (canonical RHF orbitals)")
    ax.plot(r_pts, e_casscf, "-o", color=C["s1"], ms=3.5, lw=1.2, label="CASSCF (orbital-optimised active space)")
    ax.plot(r_pts, e_fci, "-.", color=C["s3"], lw=1.2, label="Full FCI (4e in 6 orb, STO-3G)")
    ax.axvline(1.595, color=C["ink3"], ls=":", lw=0.9, label="R_e = 1.595 Å")
    ax.set_xlabel("internuclear separation R (Å)")
    ax.set_ylabel("ground state energy (Ha)")
    ax.legend(fontsize=7.2, loc="upper right")
    V["fig_lih_curve"] = savefig(fig, "lih_curve")

    lrows = []
    for r, ehf, ecas, ecasscf, efci in zip(r_pts, e_hf, e_casci, e_casscf, e_fci):
        lrows.append([f"{r:.3f}", f"{ehf:.6f}", f"{ecas:.6f}", f"{ecasscf:.6f}", f"{efci:.6f}", f"{(ecasscf - efci) * 1000:+.3f}"])
    V["lih_table"] = table(["R (Å)", "E_RHF (Ha)", "E_CASCI (Ha)", "E_CASSCF (Ha)", "E_FCI (Ha)", "Δ(CASSCF − FCI) (mHa)"], lrows)
    gaps = [(a - b) * 1000 for a, b in zip(e_casscf, e_fci)]
    gaps_ci = [(a - b) * 1000 for a, b in zip(e_casci, e_fci)]
    V["lih"] = {"gap_lo": min(gaps), "gap_hi": max(gaps), "ci_lo": min(gaps_ci), "ci_hi": max(gaps_ci),
                "rhf_end": (e_hf[-1] - e_fci[-1]) * 1000, "r_end": r_pts[-1]}
    if lih is not None:
        V["lih"].update({"corr": (lih["e_hf"] - lih["e_active"]) * 1000, "corr_full": (lih["e_hf"] - lih["e_fci_full"]) * 1000,
                         "gap_re": (lih["e_active"] - lih["e_fci_full"]) * 1000})
        lr = sorted((r for r in lih["runs"] if r["factory"] == "cultivation"), key=lambda r: r["d"])
        ld = sorted((r for r in lih["runs"] if r["factory"] == "15to1"), key=lambda r: r["d"])
        V["lih"].update({"best": lr[-1], "dist_best": max(ld, key=lambda r: r["scores"][0]) if ld else None})

    # LiH QPE simulation table
    qpe_lrows = []
    if lih and "runs" in lih:
        for r in lih["runs"]:
            fac_name = "cultivated" if r["factory"] == "cultivation" else "15-to-1 distilled"
            acc_active = f"{100 * r['scores'][0]:.1f}% ± {100 * r['score_sigmas'][0]:.1f}%"
            acc_fci = f"{100 * r['scores'][2]:.1f}%" if len(r["scores"]) >= 3 else "—"
            ef = r["expected_faults"]
            ef_str = f"{ef:.2e}" if ef < 0.1 else f"{ef:.2f}"
            qpe_lrows.append([
                f"d = {r['d']}",
                fac_name,
                f"{r['p'] * 100:.1f}%",
                acc_active,
                acc_fci,
                ef_str,
                fmt_q(r["physical_qubits"]["total"]),
                fmt_time(r["seconds"]),
            ])
    V["lih_qpe_table"] = table(
        ["distance", "factory", "noise p", "chem. acc. (CASSCF)", "chem. acc. (full FCI)", "E[faults]", "physical qubits", "quantum time"],
        qpe_lrows
    )

    # -- the text ------------------------------------------------------------------------------
    sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True, cwd=ROOT).stdout.strip()
    s15 = V["shor"].get(15)
    md = f"""---
title: Factoring and chemistry on a simulated fault-tolerant quantum computer
subtitle: Algorithms run end to end on a surface-code machine whose logical error rates were measured by circuit-level simulation
author: Jasper Sands
date: {datetime.date.today().isoformat()}
version: 0.3
commit: {sha}
description: Technical report of the ft-algorithms project.
toc: true
abstract: |
  Shor's algorithm (full modular exponentiation, with textbook ripple-carry arithmetic and with windowed
  arithmetic using measurement-based uncomputation; semiclassical approximate QFT) and iterative phase
  estimation of molecular ground-state energies (H₂, HeH⁺ and a small LiH active space) are run end to end
  on a simulated rotated-surface-code machine. Every logical error channel, idling, lattice-surgery merges
  and the surgery CNOT, was measured in {V['n_exp']} circuit-level experiments ({V['shots'] / 1e9:.2f}×10⁹ shots)
  with stabilizer-qec under SD6 noise and extracted with baselines propagated through each operation; magic
  states come from published factory models. Whole surgery experiments predicted from their parts agree with
  circuit-level measurements to {V['comp_lo']:.2f}–{V['comp_hi']:.2f}× (median {V['comp_med']:.2f}×). Windowed
  arithmetic cuts Toffolis by {V['mod']['cut_lo']:.1f}–{V['mod']['cut_hi']:.1f}× and rounds by
  {V['mod']['rcut_lo']:.1f}–{V['mod']['rcut_hi']:.1f}× on the simulated moduli. Extrapolated with the same error
  model to RSA-2048, textbook arithmetic needs d = {V['tb2048']['d']}, {fmt_q(V['tb2048']['physical_qubits'])}
  physical qubits and {fmt_time(V['tb2048']['seconds'])}; windowed arithmetic {fmt_time(V['mod2048']['seconds'])}
  with {fmt_q(V['mod2048']['physical_qubits'])}. At their authors' own failure budgets, the model gives Lee
  et al.'s FeMoco estimate its published d = {V['femoco']['estimate']['d']}, and Gidney's 2025 RSA-2048
  estimate d = {sens['baseline']['d']} rather than {sens['published_d']}: within the measured distances SD6
  memories already fail {V['g']['ratio_lo']:.1f}–{V['g']['ratio_hi']:.1f}× more often than the projection his
  distance rests on, and their fitted suppression is weaker.
---

# Introduction

Resource estimates for fault-tolerant algorithms usually multiply a gate count by an assumed
logical error rate. This work runs the algorithms instead. A logical-level simulator (ftsim,
written for this project) executes the actual circuits, with every operation followed by its
measured logical channel and every waiting patch accumulating measured idle error. The results
are success probabilities as functions of code distance, physical error rate and magic-state
source, with error budgets that say which operations cause the failures.

The approach rests on one fact and one approximation. The fact: under Pauli circuit noise, for a
Clifford circuit decoded into a Pauli-frame correction, the logical action of one surface-code
operation is exactly a Pauli channel, since every syndrome history ends in some logical Pauli. The
approximation: consecutive operations' channels are composed as independent. Section 3.2
measures how good that is.

Contributions:

- A calibrated logical model of a rotated surface code under SD6 noise: idle (including logical Y),
  Z⊗Z and X⊗X merges and the surgery CNOT, from {V['n_exp']} experiments, with fits in d and p.
- A composition test: whole lattice-surgery experiments predicted from parts at the logical level
  and compared with circuit-level measurement.
- Shor's algorithm with honest arithmetic (no use of the period anywhere in the circuit) run on
  that model for N up to {max(Ns)}, scored by peak probability because, for small N, uniformly random
  outcomes already "factor" with probability {base[15]['random'][1]:.2f} (N = 15); and the same with
  windowed arithmetic and measurement-based uncomputation.
- Phase estimation of molecular ground-state energies to chemical accuracy, from Hamiltonians computed
  from scratch (H₂, HeH⁺) and from a frozen-core CASSCF active space (LiH).
- Extrapolation to RSA-2048 and FeMoco, and a check of two published estimates' code distances against
  the measured model at their own failure budgets.

# Methods

## Circuit-level calibration

All physical-level simulation uses stabilizer-qec 1.2 (the author's engine, which matches Stim and
PyMatching) with SD6 noise at p ∈ {{0.1, 0.2, 0.3, 0.5}}% and correlated matching. Experiments: memories
of d, 2d, 3d and 4d rounds in both bases (d = 3 … 11), reference-qubit memories (d ≤ 7), Z⊗Z and X⊗X
merges of d rounds with d rounds apart before and after, the surgery CNOT (two merges and an ancilla
patch, 4d rounds in all), and, for the composition test, three Z⊗Z measurements in a row and a
three-patch merge. Each ran to 500 failures or a 300-second cap; every record keeps its seed.

**Channels through fidelities.** A distribution over flip patterns has Pauli fidelities
f(s) = Σₓ P(x)(−1)^(s·x), which multiply under composition. The per-round idle channel is the slope
of log f over the memory lengths (preparation and readout cancel). An operation's own channel is the
measured one divided by its baseline, where the baseline (preparation plus d rounds before, d rounds
plus readout after, each the square root of a 2d-round memory's fidelity) is first pushed through the
ideal operation: a CNOT copies an X error on the control onto the target and a Z on the target onto the
control, and an error before a merge changes its outcome. Dividing naively gives negative
probabilities for "X on the control only". Uncertainties are bootstrap (300 resamples).

**Logical Y.** X- and Z-basis memories see X∪Y and Z∪Y. A reference-qubit experiment gives the full
channel: the patch is projected noiselessly (every Z check by a noiseless Pauli-product measurement,
which round 1 is compared with), entangled with a perfect qubit by a noiseless Z_R Z_L measurement,
run, and ended by a perfect round and noiseless X_R X_L and Z_R Z_L measurements. Without the initial
projection an X error before the first Z-check readout flips Z_R Z_L unseen.

**Fits.** Each component ε(d, p) is fitted per p as log ε = a + b(d + 1)/2 by weighted least squares
over points with at least ten events, and globally as A(p/p*)^((d+1)/2). Values beyond the largest
measured distance are extrapolations and are marked as such in every figure.

## The machine

{table(["operation", "rounds", "logical channel"], [
        ["X, Y, Z", "0", "none (Pauli frame)"],
        ["H", "0", "none: transversal, orientation tracked in the fast block †"],
        ["S, S†", "d", "Z⊗Z merge channel (|Y⟩ resource); wrong outcome → Z"],
        ["T, T†", "d", "merge channel + Z w.p. ε_T † + wrong outcome → Z w.p. ½ (twirled S)"],
        ["CX, CZ", "2d", "measured surgery-CNOT channel, idling included"],
        ["CCX, CCZ", "2d", "3 merges + Z-type w.p. ε_CCZ † + wrong outcomes → twirled CZ; d rounds idle"],
        ["M, R", "1", "one round of idling"],
        ["idle", "per round", "measured (p_X, p_Y, p_Z)"],
        ["feed-forward", "10", "decoder reaction time before a TABLE †"]], ["l", "r", "l"])}

† marks cited inputs. Magic states: injection (ε_T = p, conservative against Li 2015's ≈ 0.4p), 15-to-1
distillation (Litinski 2019, Table 1: 4.5×10⁻⁸ per T and 5.2×10⁻¹¹ per CCZ at p = 0.1%) and cultivation
(Gidney, Shutty and Jones 2024: 2×10⁻⁹ at p = 0.1%, 4×10⁻¹¹ at 0.05%; |CCZ⟩ by 8T-to-CCZ, 28 ε_T²). Physical
qubits: Litinski's fast block, 2n + ⌈√(8n)⌉ + 1 tiles of 2(d + 1)² qubits, plus factories sized by their
published qubit·rounds per state. One round is 1 µs.

The noise compiler gives every qubit a clock. An operation starts when its qubits are free (a TABLE
also waits for its records plus the reaction time); each waiting qubit receives the composed idle
channel for its gap; the operation's channel follows it with a class label for attribution. REPEAT
blocks are barriers, scheduled once; TABLE cases are padded to the longest, so the schedule does not
depend on outcomes. Channel composition is written with log1p/expm1: per-round rates reach 10⁻²⁰ at
large d and a naive 1 − (1 − 2q)ⁿ cancels to zero.

## The simulator

ftsim (Rust) runs programs in a Stim-like format: Clifford + T + Toffoli gates, Pauli channels on
one to three qubits, measurements, REPEAT, and TABLE (classical feed-forward choosing a block by
earlier records). Backends: a sparse state (hash map of basis states; runs of permutation and phase
gates fused per basis state), which handles 36-qubit Shor circuits in milliseconds per shot, and a
dense vector. Each shot's RNG derives from (seed, shot), so results do not depend on threads.

**Stratified sampling.** With K the number of static faults, P(success) = Σₖ P(K = k)·Sₖ. P(K = k) is
computed exactly (a truncated generating polynomial over the site tree, REPEAT by exponentiation), and Sₖ
is estimated from shots with exactly k faults, drawn by sampling k sites with weights q/(1 − q) and
rejecting duplicates, which is exactly the conditional law. Strata beyond the tail cutoff carry their
mass as an interval. The one-fault stratum, with faults logged, gives each operation class's harm.

## Algorithms

**Textbook Shor.** Cuccaro ripple-carry adders; modular addition of a classical constant with a sign flag
(Vedral–Barenco–Ekert, Beauregard); controlled modular multiplication with uncomputation by the inverse
and a controlled swap; 3n + 6 qubits and 2n(2n(10n + 12) + n) Toffolis. One control qubit is recycled
through 2n rounds of a semiclassical approximate QFT whose correction rotations, chosen by TABLE from
the previous ⌈log₂ 2n⌉ + 2 records, are synthesized by gridsynth to ε = 10⁻³/2n and verified as matrices.
The base is a = 2 for every N. Every multiplier is the generic circuit, even when its constant is 1.

**Windowed Shor with measurement-based uncomputation.** The adder computes each carry with a Toffoli
into a fresh qubit and, instead of uncomputing it with a second Toffoli, measures it in the X basis
and applies CZ to the carry's two inputs when the outcome is 1 (Gidney 2018), a TABLE in the program;
subtraction is addition under bitwise complement. The exponent is read two bits at a time: for each
window the multiple a^(2^j·m) mod N (m ∈ {{1, 2, 3}}) is looked up into an (n + 1)-qubit register and
added under modular reduction, halving the number of controlled multiplications. The circuit has
4n + 7 qubits; its Toffoli count is a closed form equal to the compiled circuits' for every simulated N.

**Chemistry.** STO-3G integrals of s Gaussians in closed form (Boys F₀), restricted Hartree–Fock with
DIIS, Jordan–Wigner with interleaved spins, and Z₂ tapering of both spin parities (Bravyi et al.).
Iterative phase estimation reads 10 bits of U = exp(−i(H − E_HF)τ) with τ = 2π, approximated by four
fourth-order Suzuki steps; controlled Pauli rotations use only uncontrolled synthesized Rz's, so their
global phases stay global.

**LiH active space.** LiH in STO-3G has six spatial orbitals. Integrals come from PySCF; the Li 1s² core
is frozen into an effective one-body term, h_pq + Σ_c [2(pq|cc) − (pc|cq)], and the two valence electrons
are placed in two orbitals optimised by CASSCF(2, 2), which after tapering leaves a 2-qubit, 9-term
Hamiltonian. Its ground energy equals PySCF's CASSCF energy to 10⁻⁹ Ha. It is a deliberately small
problem: its value here is a third molecule through the same pipeline, not chemical scale.

# Results

## The calibrated machine

{V['n_exp']} experiments, {V['shots'] / 1e9:.2f}×10⁹ shots and {V['fails']:,} logical failures. Suppression
factors Λ (per two units of distance):

{V['lam_table']}

![Per-round logical error of an idle patch (X∪Y). Points measured (1σ), lines per-p fits, dashed beyond the largest measured distance.]({V['fig_idle']})

At p = 0.1% an idle patch errs {sci(V['idle_d9'])} per round at d = 9 (X∪Y), and the fit gives {sci(V['idle_d25'])}
at d = 25 for X∪Y alone ({sci(V['g']['r25']['measured'])} for any Pauli; Section 3.8). Logical Y is {100 * V['y_frac_d3']:.0f}% of the X∪Y rate at d = 3, {V['y_vs_indep_d3']:.0f}× what
independent X and Z would give: Y errors are included as measured.

**Cross-check.** {V['xc_n']} experiments rerun with Stim and PyMatching (correlated) agree with
stabilizer-qec within statistics (largest |z| = {V['xc_maxz']:.1f}).

## Composition

![Observable failure rates of whole surgery experiments: predicted from calibrated parts by ftsim vs measured at circuit level; dashed ±25%.]({V['fig_comp']})

Across {V['comp_n']} experiment–noise–distance combinations the predicted failure rate is
{V['comp_lo']:.2f}–{V['comp_hi']:.2f}× the measured one, median {V['comp_med']:.2f}×. The independent-composition
model is accurate to about 10% and leans pessimistic.

## Shor's algorithm

At p = 0.1% with cultivated magic states:

{V['shor_table']}

![Left: peak probability vs d (cultivated states, p = 0.1%); dotted: noiseless. Right: the same runs scored by factors found; dotted: random-outcome baseline.]({V['fig_shor']})

Factoring success is a poor score for small N. Classical post-processing that tries the
convergents' small multiples finds the factors from most outcomes, and the noisiest runs look best.
The peak probability falls to its random baseline (r/2^m) instead.

![N = {fN}: peak probability vs d for each magic-state source; solid p = 0.1%, dashed p = 0.2%.]({V['fig_fact']})

{V['fact_table']}

**Error budget** (N = {fN}, d = {V.get('budget_d', '—')}, cultivated, p = 0.1%):

{V.get('budget_rows', '')}

Idling accounts for {100 * V.get('idle_share', 0):.0f}% of the expected faults. The textbook circuit runs one Toffoli at a
time (the schedule averages {st['rounds_per_toffoli']:.2f}·d rounds per Toffoli with {st['busy']:.1f} patches busy)
while every other patch waits. The last column is the peak probability of runs with exactly one fault
of that class: a single logical fault is far from always fatal.

## Windowed arithmetic

{V['shor_mod_table']}

![Windowed (solid) vs textbook (dashed) Shor on the simulated machine (cultivated magic states, p = 0.1%). Left: peak probability vs d. Right: rounds.]({V['fig_shor_modern']})

On N = {', '.join(str(N) for N in V['mod']['hit90'])} the windowed circuit uses {V['mod']['cut_lo']:.1f}–{V['mod']['cut_hi']:.1f}× fewer
Toffolis and {V['mod']['rcut_lo']:.1f}–{V['mod']['rcut_hi']:.1f}× fewer rounds at d = 11, for
{min(V['mod']['extra_qubits'].values())}–{max(V['mod']['extra_qubits'].values())} more logical qubits. Rounds fall less than
Toffolis because each measurement-based uncomputation waits for the decoder before its CZ correction.
Shorter runs idle less: the distance for 90% of the noiseless peak is
{'; '.join(f"N = {N}: {t} → {m}" for N, (t, m) in V['mod']['hit90'].items())}.

## Phase estimation

{V['qpe_table']}

![Probability of chemical accuracy (1.6 mHa from FCI) vs d at p = 0.1%; dotted: noiseless.]({V['fig_qpe']})

Each run uses ≈ {sci(V['h2']['t_count']) if V['h2'] else '—'} T gates for H₂, so the magic-state error sets a floor
that no distance removes. With cultivated states (2×10⁻⁹) that floor is a few percent; with 15-to-1
distilled states (4.5×10⁻⁸) most runs fail.
{f'''
![The answer itself: energies returned by {V.get('qpe_hist_n')} individual H₂ runs at several distances (cultivated states, p = 0.1%), as the error from the exact energy; shaded: chemical accuracy; grey: the same circuit without noise.]({V['fig_qpe_hist']})
''' if V['fig_qpe_hist'] else ''}

## LiH

{V['lih_table']}

![LiH in STO-3G: restricted Hartree–Fock, CASCI(2, 2) on canonical RHF orbitals, CASSCF(2, 2) with the 1s² core frozen, and full FCI. Dotted: R = 1.595 Å, where phase estimation is run.]({V['fig_lih_curve']})

With the core frozen, CASSCF(2, 2) stays within {V['lih']['gap_lo']:.1f}–{V['lih']['gap_hi']:.1f} mHa of full STO-3G FCI along the
curve; CASCI on canonical RHF orbitals is {V['lih']['ci_lo']:.0f}–{V['lih']['ci_hi']:.0f} mHa above it, and RHF, which cannot
dissociate the bond correctly, is {V['lih']['rhf_end']:.0f} mHa above it at R = {V['lih']['r_end']:.1f} Å. At R = 1.595 Å the active space
recovers {V['lih'].get('corr', float('nan')):.1f} of the {V['lih'].get('corr_full', float('nan')):.1f} mHa correlation energy and lies
{V['lih'].get('gap_re', float('nan')):.1f} mHa above FCI, inside chemical accuracy but using most of it.

Phase estimation of the active-space Hamiltonian ({sci(lih['t_count']) if lih else '—'} T gates per run; 400 shots per
configuration, at least 60 per fault stratum):

{V['lih_qpe_table']}

The noiseless circuit is within chemical accuracy of the active-space energy {100 * lih['noiseless'][0] if lih else float('nan'):.0f}% of the
time. With cultivated states the best distance reaches {100 * V['lih']['best']['scores'][0]:.0f}% (d = {V['lih']['best']['d']});
with 15-to-1 distilled states the success saturates near {100 * V['lih']['dist_best']['scores'][0]:.0f}%, the magic-state
floor of ≈ {V['lih']['dist_best']['expected_faults']:.1f} expected faults per run.

## To scale

Operation counts of both RSA-2048 circuits are exact; the schedule structure is fitted on compiled
instances (rounds per Toffoli: {st['rounds_per_toffoli']:.2f} textbook vs. {st_mod['rounds_per_toffoli']:.2f} windowed; patches busy:
{st['busy']:.2f} vs. {st_mod['busy']:.2f}). The smallest odd d with at most {sc['target_expected_faults']:g} expected faults
(p = 0.1%, cultivated magic states):

{V['scale_table']}

![Shor scaling from the simulated sizes to RSA-2048 (p = 0.1%). Solid: textbook; dashed: windowed. The diamond is Gidney 2025's published estimate, which stores idle qubits in yoked surface codes.]({V['fig_scale']})

Windowed arithmetic shortens the RSA-2048 run from {fmt_time(V['tb2048']['seconds'])} to {fmt_time(V['mod2048']['seconds'])} but needs
more logical qubits, so it is not smaller. Both are far from Gidney's 2025 design, whose arithmetic,
layout and storage are different in kind; this section measures the error model, not that design.

**Published estimates.** For published algorithms the logical counts and run time come from the paper and
only the error model is ours. Each is compared at its authors' own failure budget and by distance only:
their physical-qubit totals depend on layouts (yoked storage, hand-made floor plans) not modelled here.
Gidney 2025 is taken one shot at a time (12.07 h, at most 1,409 logical qubits): his budget is per shot,
and a failed shot is repeated.

{V['pub_table']}

## Gidney 2025's distance under the measured model

Gidney chooses d = {sens['published_d']} because his projection of surface-code memory (the legend of his Figure 6:
p~L~ ≈ 3.5^−d^/40 per patch per round, from the yoked-surface-code simulations of Gidney, Newman, Brooks
and Jones) first reaches 10⁻¹⁵ per patch per round there; over a 12-hour shot of up to 1,600
logical qubits that leaves 93.3% of shots free of logical errors. Asked the same question, the measured
model answers d = {pr['d_measured']} (d = {pr['d_band']['-1']}–{pr['d_band']['1']} within the fit's 1σ band):

{V['per_round_table']}

\\* extrapolated beyond the measured distances (d ≤ {V['g']['dmax']}). Bands are the weighted least-squares fits'
parameter covariance, both bases moved together.

![Logical error per patch per round of an idle patch at p = 0.1%: the measured model (points measured, line fitted, dashed beyond the largest measured distance, shaded 1σ) and Gidney's projection. Dotted: the 10⁻¹⁵ he designs to.]({V['fig_gidney']})

The gap is not only extrapolation. At the measured distances SD6 memories decoded by correlated matching
already fail {V['g']['ratio_lo']:.1f}–{V['g']['ratio_hi']:.1f}× more often than the projection, and the fitted suppression per two
units of distance (Λ = {V['g']['lams'][0]:.1f} and {V['g']['lams'][-1]:.1f} for the two bases) is below its
{sens['gidney_lambda']:.2f}, so the ratio grows to {V['g']['r25']['measured'] / V['g']['r25']['gidney']:.0f}× at d = {sens['published_d']}. The two
are fits to different simulations (circuits, noise-model details and decoders all differ); this work does
not isolate which difference matters.

Running one shot of his algorithm through the estimator of this report:

{V['scen_table']}

At his budget the measured model needs d = {sens['baseline']['d']} (d = {V['g']['b90']['d']} if 90% of shots may
succeed). The fit's uncertainty and a 15% steeper extrapolation move it by at most two; the magic states
do not move it at all, since idling is {100 * sens['baseline']['breakdown']['idle'] / sum(sens['baseline']['breakdown'].values()):.0f}% of the expected faults. Given his projection
instead, the same estimator returns his d = {sens['gidney_projection']['d']}. Lee et al.'s FeMoco estimate assumes a suppression
of 10 per two units of distance, close to the measured {V['g']['lams'][0]:.1f}–{V['g']['lams'][-1]:.1f}, and under the measured model it
comes out at their d = {V['femoco']['estimate']['d']}. What the comparison
shows is how much a distance chosen by extrapolation depends on the extrapolation; it is not a measurement of
Gidney's design, which would need his layout, its yoked storage and his decoder simulated directly.

# Validation

- ftsim: dense = sparse on random circuits; 100 random circuits equal Qiskit's statevector to 10⁻¹⁰;
  noisy distributions equal Aer's exact density matrices (χ², plain and stratified); WebAssembly = native.
- Arithmetic, textbook and windowed, verified on every input by an independent bit-level simulator,
  including the measurement-based uncomputation and its feed-forward corrections; Shor's noiseless
  outcomes match the exact distribution; synthesized rotations verified as matrices; QPE matches its
  exact Trotterized distribution.
- Scaling: the Toffoli and qubit counts used at RSA-2048 equal the compiled circuits' for every simulated N;
  the estimator, given Gidney 2025's projection and budget, returns his d = 25.
- Chemistry: integrals, HF and FCI equal PySCF to 10⁻⁸ Ha; Szabo and Ostlund's H₂ reproduced; tapering
  preserves the spectrum; the LiH active-space Hamiltonian's ground energy equals PySCF's CASSCF(2, 2) to
  10⁻⁹ Ha (and CASCI on RHF orbitals to 10⁻¹¹ Ha).
- Calibration: Stim + PyMatching cross-check; reference-qubit experiment against plain memories;
  composition test.

# Limitations

- Composition is accurate to about 10% (Section 3.2), measured only for d ≤ 7.
- Logical channels are measured to d = 11 (idling) and d = 9 (surgery) at p = 0.1%. Most algorithm
  results, and every RSA-2048 and FeMoco number, lie beyond that, along exponential fits; they are marked.
- Magic-state error rates and footprints, the H and S operation model and the timing constants are
  cited inputs, not measured here.
- Published algorithms enter only through their gate counts and run times; their layouts (Gidney 2025's
  yoked storage and parallel lookups, Lee et al.'s floor plan) are not modelled, so only distances are compared.
- The LiH problem is a 2-qubit active space; chemistry at scale is represented only by FeMoco's published counts.
- Noise is SD6 throughout; biased noise and the codes tailored to it (XZZX) are not modelled.

# Reproduction

```
tools/build.sh                     # the engine (Rust → Python module)
python tools/calibrate.py          # circuit-level experiments (resumable)
python tools/fit_model.py          # data/calibration/model.json
python tools/compose_check.py; python tools/xcheck.py
python tools/noise_models.py; python tools/noise_models.py --summary
python tools/run_shor.py; python tools/run_shor_modern.py
python tools/run_qpe.py; python tools/run_lih.py
python tools/scale.py; python tools/sensitivity.py
python tools/site_data.py; python tools/report.py
```
"""
    (REPORT / "report.md").write_text(md)
    subprocess.run(["pandoc", str(REPORT / "report.md"), "--from", "markdown+pipe_tables+superscript", "--to", "html5", "--standalone",
                    "--template", str(REPORT / "template.html"), "--toc", "--toc-depth=2", "--number-sections", "-o", str(REPORT / "report.html")],
                   check=True, cwd=REPORT)
    if not args.no_pdf:
        import tempfile
        tmp_profile = tempfile.mkdtemp()
        try:
            subprocess.run([
                CHROME, "--headless=new", "--disable-gpu", "--no-pdf-header-footer",
                f"--user-data-dir={tmp_profile}", "--no-first-run", "--no-default-browser-check",
                f"--print-to-pdf={REPORT / 'report.pdf'}", (REPORT / "report.html").as_uri()
            ], check=True, capture_output=True, timeout=30)
        except Exception as e:
            print(f"Warning: PDF generation failed or timed out: {e}")
        finally:
            shutil.rmtree(tmp_profile, ignore_errors=True)
    # The site serves the report from site/report/.
    dst = ROOT / "site" / "report"
    if dst.exists():
        shutil.rmtree(dst)
    dst.mkdir(parents=True)
    for name in ("report.html", "report.css", "report.pdf"):
        if (REPORT / name).exists():
            shutil.copy(REPORT / name, dst / name)
    shutil.copytree(FIGS, dst / "figures")
    print("report.md, report.html" + ("" if args.no_pdf else ", report.pdf") + " (copied to site/report/)")


if __name__ == "__main__":
    main()
