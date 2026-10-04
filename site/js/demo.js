// The browser demo: Shor's algorithm for N = 15 on the simulated machine at a chosen distance,
// physical error rate and magic-state source.
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
export function explain(y, m, N, a) {
  const M = 1 << m;
  if (y === 0) return { text: `y = 0 → no information`, ok: false };
  // Continued-fraction convergents of y / M.
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

export function initDemo(demo) {
  const $ = (id) => document.getElementById(id);
  const state = { d: 13, p: 0.001, f: "cultivation" };
  const cfgOf = () => demo.configs.find((c) => c.d === state.d && c.p === state.p && c.factory === state.f);
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
  const ds = [...new Set(demo.configs.map((c) => c.d))].sort((a, b) => a - b);
  const ps = [...new Set(demo.configs.map((c) => c.p))].sort();
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
  let counts = new Array(1 << demo.m).fill(0), total = 0, nextShot = 0;
  const log = $("demo-log");

  function key() { return cfgOf().file; }

  function drawHist() {
    const M = 1 << demo.m;
    const W = 640, H = 220, L = 40, B = 26, T = 8;
    const iw = W - L - 8, ih = H - B - T;
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
      const ideal = demo.ideal[y];
      if (ideal > 0.002) {
        const hh = (ideal / ymax) * ih;
        s += `<rect x="${L + y * bw - 1}" y="${T + ih - hh}" width="${bw + 2}" height="${hh}" fill="none" stroke="${cssVar("var(--ink-3)")}" stroke-dasharray="2 2"/>`;
      }
      if (!counts[y]) continue;
      const v = counts[y] / total;
      const hh = (v / ymax) * ih;
      s += `<rect x="${L + y * bw}" y="${T + ih - hh}" width="${Math.max(bw, 1.2)}" height="${hh}" fill="${cssVar(demo.peaks[y] ? "var(--accent)" : "var(--ink-2)")}"><title>y = ${y}: ${counts[y]} shots</title></rect>`;
    }
    s += `<line x1="${L}" x2="${L + iw}" y1="${T + ih}" y2="${T + ih}" stroke="${cssVar("var(--ink)")}"/>`;
    for (const y of [0, 64, 128, 192, 255]) s += `<text x="${L + (y + 0.5) * bw}" y="${H - 9}" text-anchor="middle" style="font:400 10.5px var(--mono);fill:${cssVar("var(--ink-3)")}">${y}</text>`;
    s += `</g></svg><div class="legend"><span style="color:${cssVar("var(--accent)")}"><i class="dot"></i><span style="color:var(--ink-2)">outcomes on a peak</span></span><span style="color:${cssVar("var(--ink-2)")}"><i class="dot"></i><span>elsewhere</span></span><span style="color:${cssVar("var(--ink-3)")}"><i class="dash"></i><span style="color:var(--ink-2)">ideal distribution</span></span></div>`;
    $("demo-hist").innerHTML = s;
  }

  function drawStats(lastMs) {
    const c = cfgOf();
    const peak = total ? counts.reduce((a, n, y) => a + (demo.peaks[y] ? n : 0), 0) / total : NaN;
    const fac = total ? counts.reduce((a, n, y) => a + (demo.factors[y] ? n : 0), 0) / total : NaN;
    const fmt = (x) => (Number.isFinite(x) ? `${(100 * x).toFixed(1)}%` : "—");
    const q = c.qubits.total;
    const randPeak = demo.peaks.reduce((a, b) => a + b, 0) / demo.peaks.length;
    const randFac = demo.factors.reduce((a, b) => a + b, 0) / demo.factors.length;
    $("demo-stats").innerHTML = `
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
    counts = new Array(1 << demo.m).fill(0);
    total = 0;
    nextShot = 0;
    drawHist();
    drawStats();
    if (ready) load();
  }

  function load() {
    busy = true;
    runBtn.disabled = true;
    log.textContent = `Fetching and parsing the noisy program (d = ${state.d}, p = ${(state.p * 100).toFixed(1)}%, ${F_LABEL[state.f]} magic states)…`;
    worker.postMessage({ type: "load", url: new URL(`../data/demo/${key()}`, import.meta.url).href, key: key() });
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
      log.textContent = `Program ready (parsed in ${m.parseMs.toFixed(0)} ms). Press run.`;
    } else if (m.type === "result") {
      if (m.key !== key()) return;
      for (const y of m.ys) counts[y]++;
      total += m.ys.length;
      nextShot += m.ys.length;
      drawHist();
      drawStats(m.ms);
      const lines = m.ys.slice(0, 14).map((y, i) => {
        const e = explain(y, demo.m, demo.N, demo.a);
        return `<div class="${e.ok ? "ok" : "bad"}">shot ${m.firstShot + i}: ${e.text}</div>`;
      });
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
    if (busy || loadedKey !== key()) return;
    busy = true;
    runBtn.disabled = true;
    worker.postMessage({ type: "run", shots: SHOTS, seed: 2026, firstShot: nextShot, key: key() });
  };

  worker.postMessage({ type: "init", wasm: new URL("../wasm/ftsim.wasm", import.meta.url).href, records: demo.records });
  refresh();
}
