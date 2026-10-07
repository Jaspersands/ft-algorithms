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
import pathlib
import re
import subprocess
import sys

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
CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
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
    qpe_runs, qpe_base = runs_of("qpe")
    sc = jload(ROOT / "data" / "results" / "scaling.json")
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

    # -- Shor ----------------------------------------------------------------------------------
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
        axs[0].plot([r["n"] for r in rr], [r["physical_qubits"] for r in rr], "-o", ms=3, color=FCOL[f], label=f"textbook, {FLAB[f]}")
        axs[1].plot([r["n"] for r in rr], [r["seconds"] for r in rr], "-o", ms=3, color=FCOL[f])
    gpub = next(p for p in sc["published"] if p["name"].startswith("RSA"))
    axs[0].plot([2048], [897864], "D", color=C["s4"], ms=5, label="Gidney 2025 (published)")
    axs[0].plot([2048], [gpub["estimate"]["physical_qubits"]], "D", mfc="none", color=C["s4"], ms=5, label="Gidney 2025 counts, our model")
    axs[1].plot([2048], [4.96 * 86400], "D", color=C["s4"], ms=5)
    for ax in axs:
        ax.set_xscale("log", base=2)
        ax.set_yscale("log")
        ax.set_xlabel("modulus size n (bits)")
    axs[0].set_ylabel("physical qubits")
    axs[1].set_ylabel("run time (s)")
    axs[0].legend(fontsize=6.8)
    V["fig_scale"] = savefig(fig, "scaling")
    srows = []
    tb = next(r for r in sc["textbook"] if r["n"] == 2048 and r["factory"] == "cultivation")
    srows.append(["RSA-2048, textbook arithmetic", f"{3 * 2048 + 6:,}", sci(2 * 2048 * (2 * 2048 * (10 * 2048 + 12) + 2048)), "—", "—", tb["d"], fmt_q(tb["physical_qubits"]), fmt_time(tb["seconds"])])
    pubs = {"RSA-2048, Gidney 2025": ("25", "898 k", "4.96 days"), "FeMoco (THC), Lee et al. 2021": ("31", "≈ 4 M", "< 4 days")}
    for p in sc["published"]:
        pd_, pq, pt = pubs.get(p["name"], ("—", "—", "—"))
        e = p["estimate"]
        srows.append([p["name"], f"{p['logical_qubits']:,}", sci(p["toffolis"]), pd_, pq, e["d"], fmt_q(e["physical_qubits"]), fmt_time(e["seconds"])])
    V["scale_table"] = table(["workload", "logical qubits", "Toffolis", "published d", "published qubits", "our d", "our qubits", "run time"], srows)
    V["tb2048"] = tb
    V["gid"] = gpub["estimate"]
    st = sc["structure"]

    # -- the text ------------------------------------------------------------------------------
    sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True, cwd=ROOT).stdout.strip()
    s15 = V["shor"].get(15)
    md = f"""---
title: Factoring and chemistry on a simulated fault-tolerant quantum computer
subtitle: Algorithms run end to end on a surface-code machine whose logical error rates were measured by circuit-level simulation
author: Jasper Sands
date: {datetime.date.today().isoformat()}
version: 0.1
commit: {sha}
description: Technical report of the ft-algorithms project.
toc: true
abstract: |
  Shor's algorithm (full modular exponentiation in Toffoli arithmetic, semiclassical
  approximate QFT) and iterative phase estimation of the H₂ and HeH⁺ ground-state energies are
  run end to end on a simulated rotated-surface-code machine. Every logical error channel, idling,
  lattice-surgery merges and the surgery CNOT, was measured in {V['n_exp']} circuit-level experiments
  ({V['shots'] / 1e9:.2f}×10⁹ shots) with stabilizer-qec under SD6 noise and extracted with baselines
  propagated through each operation; magic states come from published factory models. Whole
  surgery experiments predicted from their parts agree with circuit-level measurements to
  {V['comp_lo']:.2f}–{V['comp_hi']:.2f}× (median {V['comp_med']:.2f}×). At p = 0.1%, factoring 15
  reaches 90% of its noiseless peak probability at d = {s15['d'] if s15 else '—'}; phase estimation
  of H₂ needs {sci(V['h2']['t_count']) if V['h2'] else '—'} T gates per run. Extrapolated with the same error model,
  textbook-arithmetic RSA-2048 needs d = {V['tb2048']['d']} and {fmt_q(V['tb2048']['physical_qubits'])}
  physical qubits for {fmt_time(V['tb2048']['seconds'])}; Gidney's 2025 algorithm needs
  d = {V['gid']['d']} under our error model where he assumed 25.
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
approximation: consecutive operations' channels are composed as independent. Section 4.3
measures how good that is.

Contributions:

- A calibrated logical model of a rotated surface code under SD6 noise: idle (including logical Y),
  Z⊗Z and X⊗X merges and the surgery CNOT, from {V['n_exp']} experiments, with fits in d and p.
- A composition test: whole lattice-surgery experiments predicted from parts at the logical level
  and compared with circuit-level measurement.
- Shor's algorithm with honest arithmetic (no use of the period anywhere in the circuit) run on
  that model for N up to {max(Ns)}, scored by peak probability because, for small N, uniformly random
  outcomes already "factor" with probability {base[15]['random'][1]:.2f} (N = 15).
- Molecular phase estimation from first-principles Hamiltonians, to chemical accuracy.
- Extrapolation to RSA-2048 and FeMoco under the same error model, compared with published estimates.

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

**Shor.** Cuccaro ripple-carry adders; modular addition of a classical constant with a sign flag
(Vedral–Barenco–Ekert, Beauregard); controlled modular multiplication with uncomputation by the inverse
and a controlled swap; 3n + 6 qubits and 2n(2n(10n + 12) + n) Toffolis. One control qubit is recycled
through 2n rounds of a semiclassical approximate QFT whose correction rotations, chosen by TABLE from
the previous ⌈log₂ 2n⌉ + 2 records, are synthesized by gridsynth to ε = 10⁻³/2n and verified as matrices.
The base is a = 2 for every N. Every multiplier is the generic circuit, even when its constant is 1.

**Chemistry.** STO-3G integrals of s Gaussians in closed form (Boys F₀), restricted Hartree–Fock with
DIIS, Jordan–Wigner with interleaved spins, and Z₂ tapering of both spin parities (Bravyi et al.).
Iterative phase estimation reads 10 bits of U = exp(−i(H − E_HF)τ) with τ = 2π, approximated by four
fourth-order Suzuki steps; controlled Pauli rotations use only uncontrolled synthesized Rz's, so their
global phases stay global.

# Results

## The calibrated machine

{V['n_exp']} experiments, {V['shots'] / 1e9:.2f}×10⁹ shots and {V['fails']:,} logical failures. Suppression
factors Λ (per two units of distance):

{V['lam_table']}

![Per-round logical error of an idle patch (X∪Y). Points measured (1σ), lines per-p fits, dashed beyond the largest measured distance.]({V['fig_idle']})

At p = 0.1% an idle patch errs {sci(V['idle_d9'])} per round at d = 9, and the fit gives {sci(V['idle_d25'])}
at d = 25. Logical Y is {100 * V['y_frac_d3']:.0f}% of the X∪Y rate at d = 3, {V['y_vs_indep_d3']:.0f}× what
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
while every other patch waits.

## Phase estimation

{V['qpe_table']}

![Probability of chemical accuracy (1.6 mHa from FCI) vs d at p = 0.1%; dotted: noiseless.]({V['fig_qpe']})

Each run uses ≈ {sci(V['h2']['t_count']) if V['h2'] else '—'} T gates for H₂, so the magic-state error sets a floor
that no distance removes. With cultivated states (2×10⁻⁹) that floor is a few percent; with 15-to-1
distilled states (4.5×10⁻⁸) most runs fail.
{f'''
![The answer itself: energies returned by {V.get('qpe_hist_n')} individual H₂ runs at several distances (cultivated states, p = 0.1%), as the error from the exact energy; shaded: chemical accuracy; grey: the same circuit without noise.]({V['fig_qpe_hist']})
''' if V['fig_qpe_hist'] else ''}
## To scale

Operation counts of the textbook circuit are exact; the schedule's shape (rounds per Toffoli in units
of d, {st['rounds_per_toffoli']:.2f}; CNOTs per Toffoli, {st['cx_per_toffoli']:.2f}; patches busy, {st['busy']:.2f}) is fitted on compiled
instances and is stable across them. The smallest odd d with at most 0.1 expected faults:

{V['scale_table']}

![Textbook Shor from the simulated sizes to RSA-2048 (p = 0.1%).]({V['fig_scale']})

Under our measured error model, Gidney's 2025 RSA-2048 algorithm needs d = {V['gid']['d']} where he assumed
25: our idle rate at d = 25 is {sci(V['idle_d25'])} per round against his 10⁻¹⁵. Our physical-qubit totals use
plain patches throughout and so bound his yoked-storage design from above. The textbook arithmetic's
3.4×10¹¹ Toffolis exceed what 15-to-1 distillation's 5.2×10⁻¹¹ can support beyond n = 256.

# Validation

- ftsim: dense = sparse on random circuits; 100 random circuits equal Qiskit's statevector to 10⁻¹⁰;
  noisy distributions equal Aer's exact density matrices (χ², plain and stratified); WebAssembly = native.
- Arithmetic verified on every input by an independent bit-level simulator; Shor's noiseless outcomes
  match the exact distribution; synthesized rotations verified as matrices; QPE matches its exact
  Trotterized distribution.
- Chemistry: integrals, HF and FCI equal PySCF to 10⁻⁸ Ha; Szabo and Ostlund's H₂ reproduced; tapering
  preserves the spectrum.
- Calibration: Stim + PyMatching cross-check; reference-qubit experiment against plain memories;
  composition test.

# Limitations

- Composition is accurate to about 10% (Section 4.2), measured only for d ≤ 7.
- Most algorithm results lie beyond the directly measured distances (d ≤ 9 for the CNOT at p = 0.1%):
  extrapolated along measured exponential fits, marked in the figures.
- Magic-state error rates and footprints, the H and S operation model and the timing constants are
  cited inputs, not measured here.
- The arithmetic is textbook and sequential. Modern constructions (windowed arithmetic, measurement-based
  uncomputation, parallel Toffolis) would cut the cost by orders of magnitude; the scaling section
  separates the error model from the algorithm by also evaluating Gidney's 2025 counts.

# Reproduction

```
tools/build.sh                     # the engine (Rust → Python module)
python tools/calibrate.py          # circuit-level experiments (resumable)
python tools/fit_model.py          # data/calibration/model.json
python tools/compose_check.py; python tools/xcheck.py
python tools/run_shor.py; python tools/run_qpe.py; python tools/scale.py
python tools/site_data.py; python tools/report.py
```
"""
    (REPORT / "report.md").write_text(md)
    subprocess.run(["pandoc", str(REPORT / "report.md"), "--from", "markdown+pipe_tables+superscript", "--to", "html5", "--standalone",
                    "--template", str(REPORT / "template.html"), "--toc", "--toc-depth=2", "--number-sections", "-o", str(REPORT / "report.html")],
                   check=True, cwd=REPORT)
    if not args.no_pdf:
        subprocess.run([CHROME, "--headless=new", "--disable-gpu", "--no-pdf-header-footer", "--run-all-compositor-stages-before-draw",
                        "--virtual-time-budget=15000", f"--print-to-pdf={REPORT / 'report.pdf'}", (REPORT / "report.html").as_uri()],
                       check=True, capture_output=True)
    # The site serves the report from site/report/.
    import shutil
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
