import { plot, stacked, fmtSci, cssVar } from "./plot.js";
import { initDemo } from "./demo.js";

const P_COLORS = { 0.001: "var(--s1)", 0.002: "var(--s2)", 0.003: "var(--s3)", 0.005: "var(--s4)" };
const N_COLORS = ["var(--s1)", "var(--s2)", "var(--s3)", "var(--s4)", "var(--s5)", "var(--s6)", "var(--s7)"];
const F_COLORS = { cultivation: "var(--s3)", "15to1": "var(--s1)", injected: "var(--s2)" };
const F_LABEL = { cultivation: "cultivated", "15to1": "15-to-1 distilled", injected: "injected" };
const BUDGET_COLORS = { idle: "var(--s7)", CNOT: "var(--s1)", CCZ: "var(--s2)", T: "var(--s3)", S: "var(--s5)", measure: "var(--s4)" };

const $ = (id) => document.getElementById(id);
const pct = (x) => `${(100 * x).toFixed(x >= 0.995 || x < 0.1 ? 1 : 0)}%`;
const fmtInt = (x) => Math.round(x).toLocaleString("en-US");
const fmtQ = (x) => (x >= 1e6 ? `${(x / 1e6).toPrecision(3)} M` : x >= 1e3 ? `${(x / 1e3).toPrecision(3)} k` : fmtInt(x));
function fmtTime(s) {
  if (s < 1) return `${(s * 1e3).toPrecision(3)} ms`;
  if (s < 120) return `${s.toPrecision(3)} s`;
  if (s < 7200) return `${(s / 60).toPrecision(3)} min`;
  if (s < 2 * 86400) return `${(s / 3600).toPrecision(3)} h`;
  if (s < 400 * 86400) return `${(s / 86400).toPrecision(3)} days`;
  return `${(s / 86400 / 365.25).toPrecision(3)} years`;
}

function setK(key, text) {
  document.querySelectorAll(`[data-k="${key}"]`).forEach((e) => (e.textContent = text));
}

function fitValue(fit, d) {
  return Math.exp(fit.a + fit.b * (d + 1) / 2);
}

function table(container, head, rows, highlight = () => false) {
  if (!container) return;
  const t = document.createElement("table");
  t.innerHTML = `<thead><tr>${head.map((h) => `<th>${h}</th>`).join("")}</tr></thead>`;
  const tb = document.createElement("tbody");
  rows.forEach((r, i) => {
    const tr = document.createElement("tr");
    if (highlight(r, i)) tr.className = "hl";
    tr.innerHTML = r.map((c, j) => `<td${j === 0 && typeof c === "string" && c.length > 18 ? ' class="t"' : ""}>${c}</td>`).join("");
    tb.appendChild(tr);
  });
  t.appendChild(tb);
  container.innerHTML = "";
  container.appendChild(t);
}

// -- calibration ---------------------------------------------------------------------------
function figIdle(cal) {
  const comp = cal.idle_xy;
  const series = [];
  for (const p of [0.001, 0.002, 0.003, 0.005]) {
    const pts = comp.points.filter((q) => q.p === p && q.n >= 10 && q.v > 0);
    series.push({ label: `p = ${(p * 100).toFixed(1)}%`, color: P_COLORS[p], line: false,
      points: pts.map((q) => ({ x: q.d, y: q.v, lo: q.v - q.s, hi: q.v + q.s, title: `d = ${q.d}: ${fmtSci(q.v)} ± ${fmtSci(q.s)} (${Math.round(q.n)} failures)` })) });
    const fit = comp.fits[String(p)] || comp.fits[p.toString()];
    if (fit) {
      const line = [];
      for (let d = 3; d <= 15; d += 1) line.push({ x: d, y: fitValue(fit, d), dashed: d > fit.dmax });
      series.push({ color: P_COLORS[p], marker: false, points: line, width: 1.4 });
    }
  }
  plot($("fig-idle"), { x: { label: "code distance d", min: 3, max: 15, ticks: [3, 5, 7, 9, 11, 13, 15] }, y: { label: "logical error per round", log: true, min: 1e-10, max: 1e-1 }, series });
}

function sumPoints(cal, names, p) {
  const by = {};
  for (const n of names) for (const q of cal[n].points.filter((q) => q.p === p)) {
    by[q.d] = by[q.d] || { v: 0, s2: 0, n: 0 };
    by[q.d].v += q.v;
    by[q.d].s2 += q.s * q.s;
    by[q.d].n += q.n;
  }
  return Object.entries(by).map(([d, o]) => ({ x: +d, y: o.v, lo: o.v - Math.sqrt(o.s2), hi: o.v + Math.sqrt(o.s2) })).filter((q) => q.y > 0);
}

function sumFit(cal, names, p, d) {
  let s = 0;
  for (const n of names) {
    const f = cal[n].fits[String(p)];
    if (f) s += fitValue(f, d);
  }
  return s;
}

function figOps(cal) {
  const p = 0.001;
  const groups = [
    ["CNOT, X part", ["cnot_x1", "cnot_x2", "cnot_x3"], "var(--s1)"],
    ["CNOT, Z part", ["cnot_z1", "cnot_z2", "cnot_z3"], "var(--s2)"],
    ["Z⊗Z merge, wrong outcome", ["zz_j1", "zz_j3", "zz_j5", "zz_j7"], "var(--s5)"],
  ];
  const series = [];
  for (const [label, names, color] of groups) {
    series.push({ label, color, line: false, points: sumPoints(cal, names, p) });
    const dmax = Math.max(...names.map((n) => cal[n].fits[String(p)]?.dmax || 0));
    series.push({ color, marker: false, width: 1.3, points: [3, 5, 7, 9, 11, 13, 15].map((d) => ({ x: d, y: sumFit(cal, names, p, d), dashed: d > dmax })) });
  }
  series.push({ label: "idle patch, 2d rounds", color: "var(--s7)", marker: false, dash: true,
    points: [3, 5, 7, 9, 11, 13, 15].map((d) => ({ x: d, y: 2 * d * (fitValue(cal.idle_xy.fits["0.001"], d) + fitValue(cal.idle_zy.fits["0.001"], d)) })) });
  plot($("fig-ops"), { x: { label: "code distance d", min: 3, max: 15, ticks: [3, 5, 7, 9, 11, 13, 15] }, y: { label: "logical error per operation", log: true, min: 1e-9, max: 1e-1 }, series });
}

function xcheckTable(xc) {
  if (!xc) return;
  const rows = xc.cases.map((c) => {
    const a = c.stabilizer_qec, b = c.stim_pymatching;
    return [c.name, fmtSci(a.failures / a.shots), fmtSci(b.failures / b.shots), `${c.z >= 0 ? "+" : "−"}${Math.abs(c.z).toFixed(1)}`];
  });
  table($("tab-xcheck"), ["experiment", "stabilizer-qec", "Stim + PyMatching", "z"], rows);
  const maxz = Math.max(...xc.cases.map((c) => Math.abs(c.z)));
  setK("xc.summary", `${xc.cases.length} of ${xc.cases.length} within ${maxz < 3 ? "3" : maxz.toFixed(1)}σ (largest |z| = ${maxz.toFixed(1)})`);
}

// -- composition -----------------------------------------------------------------------------
function figCompose(comp) {
  if (!comp) return;
  const kinds = { "cnot-z": ["CNOT (Z inputs)", "var(--s1)"], "cnot-x": ["CNOT (X inputs)", "var(--s6)"], repeated_zz_k3: ["three Z⊗Z in a row", "var(--s2)"], line_n3: ["three-patch merge", "var(--s3)"] };
  const series = Object.entries(kinds).map(([k, [label, color]]) => ({ label, color, line: false, r: 3, points: [] }));
  const ratios = [];
  for (const r of comp.results) {
    Object.entries(kinds).forEach(([k], i) => {
      const e = r.experiments[k];
      if (!e) return;
      for (const [lab, v] of Object.entries(e)) {
        if (lab === "any") { ratios.push(v.predicted / v.measured); continue; }
        if (v.measured > 0 && v.predicted > 0) series[i].points.push({ x: v.measured, y: v.predicted, title: `d = ${r.d}, p = ${r.p}: ${lab}` });
      }
    });
  }
  const diag = [1e-4, 1].map((x) => ({ x, y: x }));
  series.push({ color: "var(--ink)", marker: false, width: 1, points: diag });
  series.push({ color: "var(--ink-3)", marker: false, dash: true, width: 1, points: diag.map((q) => ({ x: q.x, y: q.y * 1.25 })) });
  series.push({ color: "var(--ink-3)", marker: false, dash: true, width: 1, points: diag.map((q) => ({ x: q.x, y: q.y / 1.25 })) });
  plot($("fig-compose"), { height: 380, x: { label: "measured (circuit level)", log: true, min: 1e-4, max: 1 }, y: { label: "predicted from parts (logical level)", log: true, min: 1e-4, max: 1 }, series });
  ratios.sort((a, b) => a - b);
  setK("comp.range", `${ratios[0].toFixed(2)}–${ratios[ratios.length - 1].toFixed(2)}×`);
  setK("comp.median", `${ratios[Math.floor(ratios.length / 2)].toFixed(2)}×`);
}

// -- Shor ------------------------------------------------------------------------------------
function shorFigures(shor) {
  const base = Object.fromEntries(shor.baselines.map((b) => [b.N, b]));
  const Ns = Object.keys(base).map(Number).sort((a, b) => a - b);
  const runs = (N, f, p) => shor.runs.filter((r) => r.N === N && r.factory === f && r.p === p).sort((a, b) => a.d - b.d);
  const color = (N) => N_COLORS[Ns.indexOf(N) % N_COLORS.length];
  const pts = (rs, i) => rs.map((r) => ({ x: r.d, y: r.scores[i], lo: r.scores[i] - 2 * r.score_sigmas[i], hi: r.scores[i] + 2 * r.score_sigmas[i], dashed: r.extrapolated,
    title: `d = ${r.d}: ${r.scores[i].toFixed(3)} ± ${(2 * r.score_sigmas[i]).toFixed(3)}; ${r.expected_faults.toPrecision(3)} expected faults` }));
  const peak = [], fac = [], hp = [], hf = [];
  for (const N of Ns) {
    const rs = runs(N, "cultivation", 0.001);
    if (!rs.length) continue;
    peak.push({ label: `N = ${N}`, color: color(N), points: pts(rs, 0) });
    fac.push({ label: `N = ${N}`, color: color(N), points: pts(rs, 1) });
    hp.push({ y: base[N].noiseless[0], color: color(N), dash: "2 3" });
    hf.push({ y: base[N].random[1], color: color(N), dash: "2 3" });
  }
  const dt = [5, 7, 9, 11, 13, 15, 17, 19, 21];
  plot($("fig-shor-peak"), { width: 420, height: 320, x: { label: "code distance d", min: 5, max: 21, ticks: dt }, y: { label: "peak probability", min: 0, max: 1.02, ticks: [0, 0.25, 0.5, 0.75, 1] }, series: peak, hlines: hp });
  plot($("fig-shor-fac"), { width: 420, height: 320, x: { label: "code distance d", min: 5, max: 21, ticks: dt }, y: { label: "factors found", min: 0, max: 1.02, ticks: [0, 0.25, 0.5, 0.75, 1] }, series: fac, hlines: hf, legend: false });
  if (base[15]) setK("shor.rand15", base[15].random[1].toFixed(2));
  if (base[21]) setK("shor.rand21", base[21].random[1].toFixed(3));

  // Table: smallest d reaching 90% of the noiseless peak probability.
  const rows = [];
  for (const N of Ns) {
    const b = base[N];
    const rs = runs(N, "cultivation", 0.001);
    const hit = rs.find((r) => r.scores[0] >= 0.9 * b.noiseless[0]);
    rows.push([`${N}`, b.logical_qubits, fmtInt(b.toffoli), fmtInt(b.t), b.noiseless[0].toFixed(3),
      hit ? hit.d : "—", hit ? fmtQ(hit.physical_qubits.total) : "—", hit ? fmtTime(hit.seconds) : "—", hit ? hit.expected_faults.toPrecision(2) : "—"]);
  }
  table($("tab-shor"), ["N", "logical qubits", "Toffolis", "T gates", "noiseless peak", "d for 90%", "physical qubits", "run time", "E[faults]"], rows);

  // Factories for one N.
  const fN = [35, 21, 15].find((N) => ["cultivation", "15to1", "injected"].every((f) => runs(N, f, 0.001).length)) || Ns[0];
  setK("shor.factN", `N = ${fN}`);
  const fs = [];
  for (const f of ["cultivation", "15to1", "injected"]) for (const p of [0.001, 0.002]) {
    const rs = runs(fN, f, p);
    if (rs.length) fs.push({ label: `${F_LABEL[f]}, p = ${(p * 100).toFixed(1)}%`, color: F_COLORS[f], dash: p === 0.002, points: pts(rs, 0) });
  }
  const inj = runs(fN, "injected", 0.001);
  if (inj.length) {
    const best = inj.reduce((a, b) => (b.scores[0] > a.scores[0] ? b : a));
    setK("shor.inj", `${best.scores[0].toFixed(2)} (noiseless ${base[fN].noiseless[0].toFixed(2)}), ${best.expected_faults.toPrecision(2)} expected faults even at d = ${best.d}`);
  }
  plot($("fig-shor-fact"), { width: 420, height: 320, x: { label: "code distance d", min: 5, max: 31 }, y: { label: "peak probability", min: 0, max: 1.02, ticks: [0, 0.25, 0.5, 0.75, 1] }, series: fs });
  const brs = runs(fN, "cultivation", 0.001).map((r) => ({ label: `d = ${r.d}`, parts: Object.fromEntries(Object.keys(BUDGET_COLORS).map((k) => [k, r.budget[k] || 0])), note: `${r.expected_faults.toPrecision(2)} faults` }));
  stacked($("fig-shor-budget"), brs, BUDGET_COLORS, { width: 420, labelW: 52, noteW: 92 });
}

function shorModernTable(shor, shorMod) {
  if (!shorMod || !shorMod.baselines) return;
  const baseTb = Object.fromEntries(shor.baselines.map((b) => [b.N, b]));
  const baseMod = Object.fromEntries(shorMod.baselines.map((b) => [b.N, b]));
  const rows = [];
  for (const N of [15, 21, 35]) {
    if (!baseTb[N] || !baseMod[N]) continue;
    const bTb = baseTb[N], bMod = baseMod[N];
    const rsTb = shor.runs.filter((r) => r.N === N && r.factory === "cultivation" && r.p === 0.001);
    const rsMod = shorMod.runs.filter((r) => r.N === N && r.factory === "cultivation" && r.p === 0.001);
    const hitTb = rsTb.find((r) => r.scores[0] >= 0.9 * bTb.noiseless[0]);
    const hitMod = rsMod.find((r) => r.scores[0] >= 0.9 * bMod.noiseless[0]);
    const d11Tb = rsTb.find((r) => r.d === 11);
    const d11Mod = rsMod.find((r) => r.d === 11);
    rows.push([`${N} (Textbook Cuccaro)`, bTb.logical_qubits, fmtInt(bTb.toffoli), d11Tb ? fmtInt(d11Tb.rounds) : "—",
      d11Tb ? d11Tb.scores[0].toFixed(3) : "—", hitTb ? hitTb.d : "—", hitTb ? fmtQ(hitTb.physical_qubits.total) : "—", hitTb ? fmtTime(hitTb.seconds) : "—"]);
    rows.push([`${N} (Modern windowed + MBU)`, bMod.logical_qubits, fmtInt(bMod.toffoli), d11Mod ? fmtInt(d11Mod.rounds) : "—",
      d11Mod ? d11Mod.scores[0].toFixed(3) : "—", hitMod ? hitMod.d : "—", hitMod ? fmtQ(hitMod.physical_qubits.total) : "—", hitMod ? fmtTime(hitMod.seconds) : "—"]);
  }
  table($("tab-shor-mod"), ["N / Arithmetic", "logical qubits", "Toffolis", "rounds (d = 11)", "peak (d = 11)", "d for 90%", "physical qubits", "run time"], rows, (r, i) => i % 2 === 1);
}

// -- chemistry -------------------------------------------------------------------------------
function qpeFigures(q) {
  if (!q || !q.runs.length) return;
  const base = Object.fromEntries(q.baselines.map((b) => [b.molecule, b]));
  const mols = Object.keys(base);
  const series = [], hl = [];
  const mc = { H2: "var(--s1)", "HeH+": "var(--s2)" };
  for (const m of mols) {
    for (const f of ["cultivation", "15to1"]) {
      const rs = q.runs.filter((r) => r.molecule === m && r.factory === f && r.p === 0.001).sort((a, b) => a.d - b.d);
      if (!rs.length) continue;
      series.push({ label: `${m.replace("2", "₂").replace("+", "⁺")}, ${F_LABEL[f]}`, color: mc[m] || "var(--s3)", dash: f !== "cultivation",
        points: rs.map((r) => ({ x: r.d, y: r.scores[0], lo: r.scores[0] - 2 * r.score_sigmas[0], hi: r.scores[0] + 2 * r.score_sigmas[0], dashed: r.extrapolated, title: `d = ${r.d}: ${r.scores[0].toFixed(3)}` })) });
    }
    hl.push({ y: base[m].noiseless[0], color: mc[m], dash: "2 3" });
  }
  plot($("fig-qpe"), { x: { label: "code distance d", min: 9, max: 29 }, y: { label: "P(chemical accuracy)", min: 0, max: 1.02, ticks: [0, 0.25, 0.5, 0.75, 1] }, series, hlines: hl });
  const rows = mols.map((m) => {
    const b = base[m];
    const rs = q.runs.filter((r) => r.molecule === m && r.factory === "cultivation" && r.p === 0.001).sort((a, b) => a.d - b.d);
    const hit = rs.find((r) => r.scores[0] >= 0.9 * b.noiseless[0]);
    return [m, b.qubits + 1, b.terms, fmtSci(b.t_count), b.e_fci.toFixed(5), b.noiseless[0].toFixed(2), hit ? hit.d : "—", hit ? fmtTime(hit.seconds) : "—", hit ? fmtQ(hit.physical_qubits.total) : "—"];
  });
  table($("tab-qpe"), ["molecule", "logical qubits", "Pauli terms", "T gates", "E_FCI (Ha)", "noiseless", "d for 90%", "run time", "physical qubits"], rows);
  if (base.H2) setK("qpe.H2.t", fmtSci(base.H2.t_count));
}

function qpeHist(h) {
  const box = $("fig-qpe-hist");
  if (!h) { box.closest("figure").style.display = "none"; return; }
  const rows = [...h.runs.map((r) => ({ label: `d = ${r.d}`, e: r.energies })), { label: "noiseless", e: h.noiseless }];
  const W = 640, rowH = 46, L = 78, R = 12, T = 6, B = 30;
  const H = T + rows.length * rowH + B;
  const lim = 12; // mHa either side
  const bw = h.resolution * 1000;
  const nb = Math.ceil((2 * lim) / bw);
  const sx = (v) => L + ((v + lim) / (2 * lim)) * (W - L - R);
  const ca = 1.6;
  let s = `<svg viewBox="0 0 ${W} ${H}" style="width:100%;display:block">`;
  s += `<rect x="${sx(-ca)}" y="${T}" width="${sx(ca) - sx(-ca)}" height="${rows.length * rowH}" fill="${cssVar("var(--ok-soft)")}"/>`;
  rows.forEach((r, i) => {
    const y0 = T + i * rowH;
    const errs = r.e.map((e) => (e - h.e_fci) * 1000);
    const counts = new Array(nb).fill(0);
    let out = 0, ok = 0;
    for (const x of errs) {
      if (Math.abs(x) < ca) ok++;
      const k = Math.floor((x + lim) / bw);
      if (k >= 0 && k < nb) counts[k]++; else out++;
    }
    const mx = Math.max(...counts, 1);
    counts.forEach((c, k) => {
      if (!c) return;
      const hgt = (c / mx) * (rowH - 12);
      s += `<rect x="${sx(-lim + k * bw)}" y="${y0 + rowH - 4 - hgt}" width="${Math.max(sx(-lim + (k + 1) * bw) - sx(-lim + k * bw) - 0.6, 1)}" height="${hgt}" fill="${cssVar(i === rows.length - 1 ? "var(--ink-3)" : "var(--s1)")}"><title>${c} runs</title></rect>`;
    });
    s += `<line x1="${L}" x2="${W - R}" y1="${y0 + rowH - 4}" y2="${y0 + rowH - 4}" stroke="${cssVar("var(--rule-soft)")}"/>`;
    s += `<text x="${L - 8}" y="${y0 + rowH / 2}" text-anchor="end" style="font:500 11px var(--mono);fill:${cssVar("var(--ink-2)")}">${r.label}</text>`;
    s += `<text x="${L - 8}" y="${y0 + rowH / 2 + 13}" text-anchor="end" style="font:400 10px var(--mono);fill:${cssVar("var(--ink-3)")}">${Math.round((100 * ok) / errs.length)}% ok${out ? `, ${out} off-scale` : ""}</text>`;
  });
  for (const t of [-12, -8, -4, 0, 4, 8, 12]) s += `<text x="${sx(t)}" y="${H - 14}" text-anchor="middle" style="font:400 10.5px var(--mono);fill:${cssVar("var(--ink-3)")}">${t}</text>`;
  s += `<text x="${(L + W - R) / 2}" y="${H - 1}" text-anchor="middle" style="font:500 12px var(--sans);fill:${cssVar("var(--ink-2)")}">E − E_FCI (mHa)</text></svg>`;
  box.innerHTML = s;
}

// -- scaling ---------------------------------------------------------------------------------
function scaleFigures(sc) {
  if (!sc) return;
  const sq = [], st = [];
  for (const f of ["cultivation", "15to1"]) {
    const rows = sc.textbook.filter((r) => r.factory === f && r.d);
    sq.push({ label: `textbook, ${F_LABEL[f]}`, color: F_COLORS[f], points: rows.map((r) => ({ x: r.n, y: r.physical_qubits, title: `n = ${r.n}: d = ${r.d}` })) });
    st.push({ label: `textbook, ${F_LABEL[f]}`, color: F_COLORS[f], points: rows.map((r) => ({ x: r.n, y: r.seconds, title: `n = ${r.n}: d = ${r.d}` })) });
    if (sc.modern) {
      const mrows = sc.modern.filter((r) => r.factory === f && r.d);
      sq.push({ label: `modern, ${F_LABEL[f]}`, color: F_COLORS[f], dash: true, points: mrows.map((r) => ({ x: r.n, y: r.physical_qubits, title: `n = ${r.n}: d = ${r.d}` })) });
      st.push({ label: `modern, ${F_LABEL[f]}`, color: F_COLORS[f], dash: true, points: mrows.map((r) => ({ x: r.n, y: r.seconds, title: `n = ${r.n}: d = ${r.d}` })) });
    }
    const lastOk = Math.max(...rows.map((r) => r.n));
    if (f === "15to1") setK("scale.distill_max", `n = ${lastOk}`);
  }
  const g = sc.published.find((p) => p.name.startsWith("RSA"));
  if (g) {
    sq.push({ label: "Gidney 2025 (published)", color: "var(--s4)", line: false, points: [{ x: 2048, y: 897864 }] });
    sq.push({ label: "Gidney 2025 counts, our error model", color: "var(--s4)", line: false, hollow: true, points: [{ x: 2048, y: g.estimate.physical_qubits }] });
    st.push({ label: "Gidney 2025", color: "var(--s4)", line: false, points: [{ x: 2048, y: 4.96 * 86400 }] });
  }
  const xt = [4, 16, 64, 256, 1024, 2048];
  plot($("fig-scale-q"), { width: 420, height: 320, x: { label: "modulus size n (bits)", log: true, min: 4, max: 2048, ticks: xt, format: String }, y: { label: "physical qubits", log: true, min: 1e4, max: 1e8 }, series: sq });
  plot($("fig-scale-t"), { width: 420, height: 320, x: { label: "modulus size n (bits)", log: true, min: 4, max: 2048, ticks: xt, format: String }, y: { label: "run time (s)", log: true, min: 0.1, max: 1e9 }, series: st, legend: false });
  const tb = sc.textbook.find((r) => r.n === 2048 && r.factory === "cultivation");
  const rows = [];
  const n = 2048;
  if (tb) rows.push(["RSA-2048, textbook arithmetic (this work)", fmtInt(3 * n + 6), fmtSci(2 * n * (2 * n * (10 * n + 12) + n)), "—", "—", tb.d, fmtQ(tb.physical_qubits), fmtTime(tb.seconds)]);
  if (sc.modern) {
    const mod = sc.modern.find((r) => r.n === 2048 && r.factory === "cultivation");
    if (mod) rows.push(["RSA-2048, modern arithmetic (k=2, MBU, this work)", fmtInt(2 * n + 2 * 2 + 20), fmtSci(86107386880.0), "—", "—", mod.d, fmtQ(mod.physical_qubits), fmtTime(mod.seconds)]);
  }
  const pubs = { "RSA-2048, Gidney 2025": { d: 25, q: "898 k", t: "4.96 days" }, "FeMoco (THC), Lee et al. 2021": { d: 31, q: "≈ 4 M", t: "< 4 days" } };
  for (const p of sc.published) {
    const pub = pubs[p.name] || {};
    rows.push([p.name, fmtInt(p.logical_qubits), fmtSci(p.toffolis), pub.d ?? "—", pub.q ?? "—", p.estimate.d, fmtQ(p.estimate.physical_qubits), fmtTime(p.estimate.seconds)]);
  }
  table($("tab-scale"), ["workload", "logical qubits", "Toffolis", "published d", "published qubits", "our d", "our qubits", "run time"], rows, (r) => r[0].includes("modern"));
}

function sensitivityTable(sens) {
  if (!sens || !sens.experiments) return;
  const rows = sens.experiments.map((e) => [
    e.category,
    e.name,
    e.d,
    fmtQ(e.physical_qubits),
    `${e.delta_d >= 0 ? "+" : ""}${e.delta_d}`
  ]);
  table($("tab-tornado"), ["category", "scenario", "required d", "physical qubits", "Δd vs baseline"], rows, (r) => r[0] === "Assumed Model");
}

function opsTable() {
  const rows = [
    ["X, Y, Z", "0", "none (Pauli frame)"],
    ["H", "0", "none: transversal, orientation tracked (fast block) †"],
    ["S, S†", "d", "Z⊗Z merge channel with a |Y⟩ resource; wrong outcome → Z"],
    ["T, T†", "d", "merge channel + Z with ε_T † + wrong outcome → Z w.p. ½ (twirled S)"],
    ["CX, CZ", "2d", "measured surgery CNOT channel (idling included)"],
    ["CCX (Toffoli)", "2d", "3 merges + Z-type with ε_CCZ † + wrong outcomes → twirled CZ; d rounds idle"],
    ["M, R", "1", "one round of idling"],
    ["idle", "per round", "measured (p_X, p_Y, p_Z) of a patch"],
    ["feed-forward", "10", "decoder reaction time before a TABLE †"],
  ];
  table($("tab-ops"), ["operation", "rounds", "logical channel"], rows);
}

async function main() {
  const data = await (await fetch("data/site.json")).json();
  const cal = data.calibration;
  setK("cal.experiments", fmtInt(data.calibration_meta.experiments));
  setK("cal.shots", `${(data.calibration_meta.shots / 1e9).toFixed(1)} billion`);
  const lam = Math.exp(-cal.idle_xy.fits["0.001"].b);
  setK("cal.lambda_idle", `Λ ≈ ${lam.toFixed(1)}`);
  const y3 = cal.idle_y.points.find((q) => q.p === 0.001 && q.d === 3), x3 = cal.idle_xy.points.find((q) => q.p === 0.001 && q.d === 3);
  if (y3 && x3) setK("cal.y_ratio", `about ${Math.round((100 * y3.v) / x3.v)}%`);
  setK("demo.qubits", data.demo.logical_qubits);
  setK("demo.toffoli", fmtInt(data.demo.toffoli));
  setK("demo.t", fmtInt(data.demo.t));
  figIdle(cal);
  figOps(cal);
  xcheckTable(data.xcheck);
  figCompose(data.composition);
  shorFigures(data.shor);
  shorModernTable(data.shor, data.shor_modern);
  qpeFigures(data.qpe);
  qpeHist(data.qpe_hist);
  scaleFigures(data.scaling);
  sensitivityTable(data.sensitivity);
  opsTable();
  initDemo(data.demo, data.demo_modern);
}

main();
