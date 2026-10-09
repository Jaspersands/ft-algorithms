// The browser demo: Shor's algorithm for N = 15 (Textbook and Modern arithmetic)
// and Molecular Chemistry (H₂ Ground State QPE) on the simulated machine.

import { cssVar } from "./plot.js";

const SHOTS = 500;
const F_LABEL = { cultivation: "cultivated", "15to1": "15-to-1 distilled", injected: "injected" };

function gcd(a, b) {
  while (b) [a, b] = [b, a % b];
  return a;
}

function powmod(a, e, n) {
  let r = 1;
  a %= n;
  while (e > 0) {
    if (e & 1) r = (r * a) % n;
    a = (a * a) % n;
    e >>= 1;
  }
  return r;
}

/** The classical post-processing of one run, mirroring ftalgo.shor.factors_from_outcome. */
export function explainShor(y, m, N, a) {
  const M = 1 << m;
  if (y === 0) return { text: `y = 0 → no information`, ok: false };
  let num = y, den = M;
  let h0 = 0, h1 = 1, k0 = 1, k1 = 0;
  const qs = [];
  while (den) {
    const q = Math.floor(num / den);
    [h0, h1] = [h1, q * h1 + h0];
    [k0, k1] = [k1, q * k1 + k0];
    if (k1 >= N) break;
    qs.push([h1, k1]);
    [num, den] = [den, num - q * den];
  }
  for (const [p, q] of qs) {
    for (let t = 1; q * t < N; t++) {
      const r = q * t;
      if (powmod(a, r, N) === 1) {
        if (r % 2 === 0) {
          const h = powmod(a, r / 2, N);
          if (h !== N - 1) {
            const f = gcd(h - 1, N);
            if (f > 1 && f < N) return { text: `y = ${y} → ${y}/${M} ≈ ${p}/${q} → r = ${r} → gcd(${a}^${r / 2} − 1, ${N}) = ${f} → ${Math.min(f, N / f)} × ${Math.max(f, N / f)}`, ok: true };
          }
        }
        return { text: `y = ${y} → r = ${r} gives no factor`, ok: false };
      }
    }
  }
  return { text: `y = ${y} → no order found`, ok: false };
}

export function explainChem(y, demoChem) {
  const E = demoChem.energies[y];
  const diffMHa = (E - demoChem.e_fci) * 1000;
  const isAcc = demoChem.chemical_accuracy[y];
  const sign = diffMHa >= 0 ? "+" : "";
  return {
    text: `y = ${y} → E = ${E.toFixed(4)} Ha (diff: ${sign}${diffMHa.toFixed(2)} mHa) → ${isAcc ? "✓ Chemical Accuracy (±1.6 mHa)" : "✗ Outside Accuracy Window"}`,
    ok: Boolean(isAcc),
  };
}

export function initDemo(demo_tb, demo_mod, demo_chem) {
  const $ = (id) => document.getElementById(id);
  const state = {
    workload: "shor_tb", // 'shor_tb', 'shor_mod', 'chem_h2'
    d: 13,
    p: 0.001,
    f: "cultivation",
  };

  const curDemo = () => {
    if (state.workload === "chem_h2" && demo_chem) return demo_chem;
    if (state.workload === "shor_mod" && demo_mod) return demo_mod;
    return demo_tb;
  };

  const isChem = () => state.workload === "chem_h2";

  const cfgOf = () => {
    const cd = curDemo();
    return cd.configs.find((c) => c.d === state.d && c.p === state.p && c.factory === state.f) || cd.configs[0];
  };

  function updateHeaderCounts() {
    const cd = curDemo();
    document.querySelectorAll('[data-k="demo.qubits"]').forEach((e) => (e.textContent = cd.logical_qubits));
    document.querySelectorAll('[data-k="demo.toffoli"]').forEach((e) => (e.textContent = (cd.toffoli || 0).toLocaleString("en-US")));
    document.querySelectorAll('[data-k="demo.t"]').forEach((e) => (e.textContent = (cd.t || 0).toLocaleString("en-US")));
  }

  const workSel = $("demo-workload") || $("demo-arith");
  if (workSel) {
    workSel.value = state.workload === "shor_mod" ? "modern" : (state.workload === "chem_h2" ? "chem_h2" : "textbook");
    workSel.onchange = () => {
      const v = workSel.value;
      if (v === "modern" || v === "shor_mod") state.workload = "shor_mod";
      else if (v === "chem" || v === "chem_h2") state.workload = "chem_h2";
      else state.workload = "shor_tb";
      updateHeaderCounts();
      refresh();
    };
  }

  const seg = (container, values, label, key) => {
    container.innerHTML = "";
    for (const v of values) {
      const b = document.createElement("button");
      b.textContent = label(v);
      b.onclick = () => { state[key] = v; refresh(); };
      b.dataset.v = v;
      container.appendChild(b);
    }
  };
  const ds = [...new Set(demo_tb.configs.map((c) => c.d))].sort((a, b) => a - b);
  const ps = [...new Set(demo_tb.configs.map((c) => c.p))].sort();
  seg($("demo-d"), ds, String, "d");
  seg($("demo-p"), ps, (p) => `${(p * 100).toFixed(1)}%`, "p");
  const sel = $("demo-f");
  for (const f of ["cultivation", "15to1", "injected"]) {
    const o = document.createElement("option");
    o.value = f;
    o.textContent = F_LABEL[f];
    sel.appendChild(o);
  }
  sel.value = state.f;
  sel.onchange = () => { state.f = sel.value; refresh(); };

  const worker = new Worker(new URL("./demo-worker.js", import.meta.url), { type: "module" });
  const runBtn = $("demo-run");
  let ready = false, loadedKey = null, busy = false;
  let counts = new Array(1 << curDemo().m).fill(0), total = 0, nextShot = 0;
  const log = $("demo-log");

  function key() { return cfgOf().file; }

  function getFolder() {
    if (state.workload === "chem_h2") return "demo_chem";
    if (state.workload === "shor_mod") return "demo_modern";
    return "demo";
  }

  function drawHist() {
    const cd = curDemo();
    const M = 1 << cd.m;
    const W = 640, H = 220, L = 46, B = 26, T = 8;
    const iw = W - L - 8, ih = H - B - T;

    if (isChem()) {
      // Chemistry Histogram: Energy Error in mHa around E_FCI
      const lim = 15; // mHa either side
      const bw = cd.resolution * 1000;
      const nb = Math.ceil((2 * lim) / bw);
      const sx = (v) => L + ((v + lim) / (2 * lim)) * iw;
      const ca = 1.6;

      const maxc = Math.max(1, ...counts) / Math.max(total, 1);
      const ymax = Math.max(0.2, Math.ceil(maxc * 10) / 10);

      let s = `<svg viewBox="0 0 ${W} ${H}" style="width:100%;display:block"><g>`;
      // Green chemical accuracy window
      s += `<rect x="${sx(-ca)}" y="${T}" width="${sx(ca) - sx(-ca)}" height="${ih}" fill="${cssVar("var(--ok-soft)")}"/>`;

      for (const t of [0, ymax / 2, ymax]) {
        const y = T + ih - (t / ymax) * ih;
        s += `<line x1="${L}" x2="${L + iw}" y1="${y}" y2="${y}" stroke="${cssVar("var(--rule-soft)")}"/>`;
        s += `<text x="${L - 6}" y="${y + 3.5}" text-anchor="end" style="font:400 10.5px var(--mono);fill:${cssVar("var(--ink-3)")}">${t.toFixed(2)}</text>`;
      }

      for (let y = 0; y < M; y++) {
        if (!counts[y]) continue;
        const diffMHa = (cd.energies[y] - cd.e_fci) * 1000;
        const v = counts[y] / total;
        const hh = (v / ymax) * ih;
        const isAcc = cd.chemical_accuracy[y];
        const barX = sx(diffMHa);
        const barW = Math.max((bw / (2 * lim)) * iw - 1, 2);
        s += `<rect x="${barX}" y="${T + ih - hh}" width="${barW}" height="${hh}" fill="${cssVar(isAcc ? "var(--accent)" : "var(--ink-2)")}"><title>y = ${y}: ${diffMHa.toFixed(2)} mHa (${counts[y]} shots)</title></rect>`;
      }

      s += `<line x1="${L}" x2="${L + iw}" y1="${T + ih}" y2="${T + ih}" stroke="${cssVar("var(--ink)")}"/>`;
      for (const t of [-12, -8, -4, 0, 4, 8, 12]) {
        s += `<text x="${sx(t)}" y="${H - 9}" text-anchor="middle" style="font:400 10.5px var(--mono);fill:${cssVar("var(--ink-3)")}">${t}</text>`;
      }
      s += `</g></svg><div class="legend"><span style="color:${cssVar("var(--ok)")}"><i class="dot"></i><span>Chemical Accuracy Window (±1.6 mHa)</span></span><span style="color:${cssVar("var(--accent)")}"><i class="dot"></i><span>accurate shots</span></span><span style="color:${cssVar("var(--ink-2)")}"><i class="dot"></i><span>outside accuracy</span></span></div>`;
      $("demo-hist").innerHTML = s;
      return;
    }

    // Shor histogram
    const maxc = Math.max(1, ...counts) / Math.max(total, 1);
    const ymax = Math.max(0.3, Math.ceil(maxc * 10) / 10);
    const bw = iw / M;
    let s = `<svg viewBox="0 0 ${W} ${H}" style="width:100%;display:block"><g>`;
    for (const t of [0, ymax / 2, ymax]) {
      const y = T + ih - (t / ymax) * ih;
      s += `<line x1="${L}" x2="${L + iw}" y1="${y}" y2="${y}" stroke="${cssVar("var(--rule-soft)")}"/>`;
      s += `<text x="${L - 6}" y="${y + 3.5}" text-anchor="end" style="font:400 10.5px var(--mono);fill:${cssVar("var(--ink-3)")}">${t.toFixed(2)}</text>`;
    }
    for (let y = 0; y < M; y++) {
      const ideal = cd.ideal ? cd.ideal[y] : 0;
      if (ideal > 0.002) {
        const hh = (ideal / ymax) * ih;
        s += `<rect x="${L + y * bw - 1}" y="${T + ih - hh}" width="${bw + 2}" height="${hh}" fill="none" stroke="${cssVar("var(--ink-3)")}" stroke-dasharray="2 2"/>`;
      }
      if (!counts[y]) continue;
      const v = counts[y] / total;
      const hh = (v / ymax) * ih;
      s += `<rect x="${L + y * bw}" y="${T + ih - hh}" width="${Math.max(bw, 1.2)}" height="${hh}" fill="${cssVar(cd.peaks && cd.peaks[y] ? "var(--accent)" : "var(--ink-2)")}"><title>y = ${y}: ${counts[y]} shots</title></rect>`;
    }
    s += `<line x1="${L}" x2="${L + iw}" y1="${T + ih}" y2="${T + ih}" stroke="${cssVar("var(--ink)")}"/>`;
    for (const y of [0, 64, 128, 192, 255]) s += `<text x="${L + (y + 0.5) * bw}" y="${H - 9}" text-anchor="middle" style="font:400 10.5px var(--mono);fill:${cssVar("var(--ink-3)")}">${y}</text>`;
    s += `</g></svg><div class="legend"><span style="color:${cssVar("var(--accent)")}"><i class="dot"></i><span style="color:var(--ink-2)">outcomes on a peak</span></span><span style="color:${cssVar("var(--ink-2)")}"><i class="dot"></i><span>elsewhere</span></span><span style="color:${cssVar("var(--ink-3)")}"><i class="dash"></i><span style="color:var(--ink-2)">ideal distribution</span></span></div>`;
    $("demo-hist").innerHTML = s;
  }

  function drawStats(lastMs) {
    const cd = curDemo();
    const c = cfgOf();
    const fmt = (x) => (Number.isFinite(x) ? `${(100 * x).toFixed(1)}%` : "—");
    const q = c.qubits.total;

    if (isChem()) {
      const accCount = total ? counts.reduce((a, n, y) => a + (cd.chemical_accuracy[y] ? n : 0), 0) / total : NaN;
      const w4Count = total ? counts.reduce((a, n, y) => a + (cd.within_4x[y] ? n : 0), 0) / total : NaN;
      $("demo-stats").innerHTML = `
        <dt>workload</dt><dd>H₂ Ground State QPE (10-bit)</dd>
        <dt>chemical accuracy (&lt; 1.6 mHa)</dt><dd>${fmt(accCount)} <small>target: &gt; 90%</small></dd>
        <dt>within 6.4 mHa</dt><dd>${fmt(w4Count)}</dd>
        <dt>expected faults per run</dt><dd>${c.expected_faults < 0.01 ? c.expected_faults.toExponential(1) : c.expected_faults.toPrecision(3)}</dd>
        <dt>fault-free runs</dt><dd>${(100 * c.p0).toFixed(c.p0 > 0.99 ? 2 : 1)}%</dd>
        <dt>machine</dt><dd>${(c.rounds / 1000).toFixed(0)} ms <small>· ${(q / 1000).toFixed(0)}k physical qubits</small></dd>
        <dt>shots</dt><dd>${total}${lastMs ? ` <small>· ${lastMs.toFixed(0)} ms in your browser</small>` : ""}</dd>
        ${c.extrapolated ? `<dt>note</dt><dd><small>distance beyond directly measured range</small></dd>` : ""}`;
      return;
    }

    const peak = total ? counts.reduce((a, n, y) => a + (cd.peaks[y] ? n : 0), 0) / total : NaN;
    const fac = total ? counts.reduce((a, n, y) => a + (cd.factors[y] ? n : 0), 0) / total : NaN;
    const randPeak = cd.peaks.reduce((a, b) => a + b, 0) / cd.peaks.length;
    const randFac = cd.factors.reduce((a, b) => a + b, 0) / cd.factors.length;
    $("demo-stats").innerHTML = `
      <dt>arithmetic</dt><dd>${state.workload === "shor_mod" ? "Modern (Windowed + MBU)" : "Textbook (Cuccaro)"}</dd>
      <dt>on a peak</dt><dd>${fmt(peak)} <small>ideal 100%, random ${(100 * randPeak).toFixed(1)}%</small></dd>
      <dt>factors found</dt><dd>${fmt(fac)} <small>random outcomes: ${(100 * randFac).toFixed(0)}%</small></dd>
      <dt>expected faults per run</dt><dd>${c.expected_faults < 0.01 ? c.expected_faults.toExponential(1) : c.expected_faults.toPrecision(3)}</dd>
      <dt>fault-free runs</dt><dd>${(100 * c.p0).toFixed(c.p0 > 0.99 ? 2 : 1)}%</dd>
      <dt>machine</dt><dd>${(c.rounds / 1000).toFixed(0)} ms <small>· ${(q / 1000).toFixed(0)}k physical qubits</small></dd>
      <dt>shots</dt><dd>${total}${lastMs ? ` <small>· ${lastMs.toFixed(0)} ms in your browser</small>` : ""}</dd>
      ${c.extrapolated ? `<dt>note</dt><dd><small>distance beyond the directly measured range: channels extrapolated</small></dd>` : ""}`;
  }

  function refresh() {
    document.querySelectorAll("#demo-d button").forEach((b) => b.classList.toggle("on", +b.dataset.v === state.d));
    document.querySelectorAll("#demo-p button").forEach((b) => b.classList.toggle("on", +b.dataset.v === state.p));
    counts = new Array(1 << curDemo().m).fill(0);
    total = 0;
    nextShot = 0;
    drawHist();
    drawStats();
    if (ready) load();
  }

  function load() {
    busy = true;
    runBtn.disabled = true;
    const cd = curDemo();
    const folder = getFolder();
    log.textContent = `Fetching and parsing the noisy program (${state.workload}, d = ${state.d}, p = ${(state.p * 100).toFixed(1)}%, ${F_LABEL[state.f]} magic states)…`;
    worker.postMessage({
      type: "load",
      url: new URL(`../data/${folder}/${key()}`, import.meta.url).href,
      key: `${state.workload}:${key()}`,
      records: cd.records,
    });
  }

  worker.onmessage = (ev) => {
    const m = ev.data;
    if (m.type === "ready") {
      ready = true;
      load();
    } else if (m.type === "loaded") {
      loadedKey = m.key;
      busy = false;
      runBtn.disabled = false;
      log.textContent = `Program ready (${state.workload} parsed in ${m.parseMs.toFixed(0)} ms). Press run.`;
    } else if (m.type === "result") {
      if (m.key !== `${state.workload}:${key()}`) return;
      for (const y of m.ys) counts[y]++;
      total += m.ys.length;
      nextShot += m.ys.length;
      drawHist();
      drawStats(m.ms);
      const cd = curDemo();
      let lines = [];
      if (isChem()) {
        lines = m.ys.slice(0, 14).map((y, i) => {
          const e = explainChem(y, cd);
          return `<div class="${e.ok ? "ok" : "bad"}">shot ${m.firstShot + i}: ${e.text}</div>`;
        });
      } else {
        lines = m.ys.slice(0, 14).map((y, i) => {
          const e = explainShor(y, cd.m, cd.N, cd.a);
          return `<div class="${e.ok ? "ok" : "bad"}">shot ${m.firstShot + i}: ${e.text}</div>`;
        });
      }
      log.innerHTML = lines.join("");
      busy = false;
      runBtn.disabled = false;
    } else if (m.type === "error") {
      log.textContent = `Error: ${m.message}`;
      busy = false;
      runBtn.disabled = false;
    }
  };

  runBtn.onclick = () => {
    if (busy || loadedKey !== `${state.workload}:${key()}`) return;
    busy = true;
    runBtn.disabled = true;
    worker.postMessage({ type: "run", shots: SHOTS, seed: 2026, firstShot: nextShot, key: `${state.workload}:${key()}` });
  };

  updateHeaderCounts();
  worker.postMessage({ type: "init", wasm: new URL("../wasm/ftsim.wasm", import.meta.url).href, records: demo_tb.records });
  refresh();
}
