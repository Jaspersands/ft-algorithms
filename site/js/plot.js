// A small SVG plotting primitive: axes (linear or log), gridlines, series of points with error
// bars and lines (dashed where marked), horizontal reference lines and a legend. Every chart on
// the page goes through it, so they cannot drift apart in style.

const NS = "http://www.w3.org/2000/svg";

function el(name, attrs = {}, parent) {
  const e = document.createElementNS(NS, name);
  for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
  if (parent) parent.appendChild(e);
  return e;
}

export function cssVar(v) {
  const m = /^var\((--[\w-]+)\)$/.exec(String(v));
  return m ? getComputedStyle(document.documentElement).getPropertyValue(m[1]).trim() : v;
}

function niceLinear(min, max, n = 5) {
  const span = max - min || 1;
  const step0 = span / n;
  const mag = 10 ** Math.floor(Math.log10(step0));
  const step = [1, 2, 2.5, 5, 10].map((s) => s * mag).find((s) => span / s <= n + 0.5) || 10 * mag;
  const ticks = [];
  for (let t = Math.ceil(min / step - 1e-9) * step; t <= max + 1e-9 * step; t += step) ticks.push(+t.toFixed(10));
  return ticks;
}

function logTicks(min, max) {
  const out = [];
  for (let e = Math.floor(Math.log10(min)); e <= Math.ceil(Math.log10(max)); e++) out.push(10 ** e);
  return out.filter((t) => t >= min * 0.999 && t <= max * 1.001);
}

export function fmtSci(v) {
  if (v === 0) return "0";
  const e = Math.floor(Math.log10(Math.abs(v)));
  if (e >= -2 && e <= 3) return String(+v.toPrecision(3));
  const m = v / 10 ** e;
  const sup = String(e).replace(/-/g, "−").split("").map((c) => "⁰¹²³⁴⁵⁶⁷⁸⁹"["0123456789".indexOf(c)] || (c === "−" ? "⁻" : c)).join("");
  return (Math.abs(m - 1) < 1e-9 ? "" : `${+m.toPrecision(2)}×`) + "10" + sup;
}

/**
 * plot(container, spec)
 * spec: { width, height, x: {label, log, min, max, ticks, format}, y: {…},
 *         series: [{label, color, points: [{x, y, lo, hi, dashed, title}], line, marker, dash}],
 *         hlines: [{y, color, label, dash}], legend: true }
 */
export function plot(container, spec) {
  const W = spec.width || 640, H = spec.height || 360;
  const M = { top: 14, right: 16, bottom: 44, left: 62, ...(spec.margin || {}) };
  const iw = W - M.left - M.right, ih = H - M.top - M.bottom;
  const all = spec.series.flatMap((s) => s.points);
  const xs = all.map((p) => p.x), ys = all.flatMap((p) => [p.y, p.lo ?? p.y, p.hi ?? p.y]).filter((v) => Number.isFinite(v) && (!spec.y.log || v > 0));
  (spec.hlines || []).forEach((h) => ys.push(h.y));
  const xmin = spec.x.min ?? Math.min(...xs), xmax = spec.x.max ?? Math.max(...xs);
  let ymin = spec.y.min ?? Math.min(...ys), ymax = spec.y.max ?? Math.max(...ys);
  if (spec.y.log) { ymin = spec.y.min ?? 10 ** Math.floor(Math.log10(ymin)); ymax = spec.y.max ?? 10 ** Math.ceil(Math.log10(ymax)); }
  const sx = spec.x.log ? (v) => M.left + (Math.log10(v) - Math.log10(xmin)) / (Math.log10(xmax) - Math.log10(xmin)) * iw : (v) => M.left + (v - xmin) / (xmax - xmin) * iw;
  const sy = spec.y.log ? (v) => M.top + ih - (Math.log10(v) - Math.log10(ymin)) / (Math.log10(ymax) - Math.log10(ymin)) * ih : (v) => M.top + ih - (v - ymin) / (ymax - ymin) * ih;
  const clampY = (v) => Math.max(M.top, Math.min(M.top + ih, v));

  container.innerHTML = "";
  const wrap = document.createElement("div");
  wrap.className = "chart";
  container.appendChild(wrap);
  const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": spec.aria || spec.y.label || "chart" }, wrap);
  const grid = el("g", { class: "grid" }, svg);
  const axis = el("g", { class: "axis" }, svg);

  const xt = spec.x.ticks || (spec.x.log ? logTicks(xmin, xmax) : niceLinear(xmin, xmax, spec.x.n || 6));
  const yt = spec.y.ticks || (spec.y.log ? logTicks(ymin, ymax) : niceLinear(ymin, ymax, spec.y.n || 5));
  const xf = spec.x.format || ((v) => (spec.x.log ? fmtSci(v) : String(v)));
  const yf = spec.y.format || ((v) => (spec.y.log ? fmtSci(v) : String(+v.toPrecision(3))));
  for (const t of yt) {
    el("line", { x1: M.left, x2: M.left + iw, y1: sy(t), y2: sy(t) }, grid);
    const g = el("g", { class: "tick" }, axis);
    el("text", { x: M.left - 8, y: sy(t) + 3.5, "text-anchor": "end" }, g).textContent = yf(t);
  }
  for (const t of xt) {
    el("line", { x1: sx(t), x2: sx(t), y1: M.top, y2: M.top + ih }, grid);
    const g = el("g", { class: "tick" }, axis);
    el("text", { x: sx(t), y: M.top + ih + 16, "text-anchor": "middle" }, g).textContent = xf(t);
  }
  el("path", { d: `M${M.left},${M.top}V${M.top + ih}H${M.left + iw}`, fill: "none" }, axis);
  el("text", { class: "label", x: M.left + iw / 2, y: H - 6, "text-anchor": "middle" }, svg).textContent = spec.x.label || "";
  el("text", { class: "label", transform: `translate(14 ${M.top + ih / 2}) rotate(-90)`, "text-anchor": "middle" }, svg).textContent = spec.y.label || "";

  for (const h of spec.hlines || []) {
    const c = cssVar(h.color || "var(--ink-3)");
    el("line", { x1: M.left, x2: M.left + iw, y1: sy(h.y), y2: sy(h.y), stroke: c, "stroke-width": 1.2, "stroke-dasharray": h.dash || "4 4" }, svg);
    if (h.label) el("text", { x: M.left + iw - 4, y: sy(h.y) - 4, "text-anchor": "end", style: `font: 500 10.5px var(--sans); fill: ${c}` }, svg).textContent = h.label;
  }

  for (const s of spec.series) {
    const c = cssVar(s.color || "var(--ink)");
    const pts = s.points.filter((p) => Number.isFinite(p.y) && (!spec.y.log || p.y > 0)).sort((a, b) => a.x - b.x);
    if (s.line !== false && pts.length > 1) {
      for (let i = 0; i + 1 < pts.length; i++) {
        const a = pts[i], b = pts[i + 1];
        el("line", { x1: sx(a.x), y1: clampY(sy(a.y)), x2: sx(b.x), y2: clampY(sy(b.y)), stroke: c, "stroke-width": s.width || 1.8,
          "stroke-dasharray": b.dashed || s.dash ? "5 4" : "none", "stroke-linecap": "round" }, svg);
      }
    }
    if (s.marker !== false) {
      for (const p of pts) {
        if (p.lo !== undefined && p.hi !== undefined) {
          const lo = spec.y.log ? Math.max(p.lo, ymin) : p.lo;
          el("line", { x1: sx(p.x), x2: sx(p.x), y1: clampY(sy(lo)), y2: clampY(sy(p.hi)), stroke: c, "stroke-width": 1.2 }, svg);
        }
        const circ = el("circle", { cx: sx(p.x), cy: clampY(sy(p.y)), r: s.r || 3.4, fill: p.dashed || s.hollow ? cssVar("var(--surface)") : c, stroke: c, "stroke-width": 1.4 }, svg);
        if (p.title) el("title", {}, circ).textContent = p.title;
      }
    }
  }

  if (spec.legend !== false) {
    const lg = document.createElement("div");
    lg.className = "legend";
    for (const s of spec.series.filter((s) => s.label)) {
      const sp = document.createElement("span");
      sp.style.color = cssVar(s.color || "var(--ink)");
      sp.innerHTML = `<i class="${s.dash ? "dash" : s.line === false ? "dot" : ""}"></i><span style="color: var(--ink-2)"></span>`;
      sp.lastChild.textContent = s.label;
      lg.appendChild(sp);
    }
    for (const h of (spec.hlines || []).filter((h) => h.legend)) {
      const sp = document.createElement("span");
      sp.style.color = cssVar(h.color || "var(--ink-3)");
      sp.innerHTML = `<i class="dash"></i><span style="color: var(--ink-2)"></span>`;
      sp.lastChild.textContent = h.legend;
      lg.appendChild(sp);
    }
    container.appendChild(lg);
  }
  return { svg, sx, sy };
}

/** Horizontal 100% stacked bars: rows [{label, parts: {name: value}, note}] with colors {name: color}. */
export function stacked(container, rows, colors, { width = 640, rowH = 26, labelW = 70, noteW = 120 } = {}) {
  const H = rows.length * rowH + 8;
  container.innerHTML = "";
  const wrap = document.createElement("div");
  wrap.className = "chart";
  container.appendChild(wrap);
  const svg = el("svg", { viewBox: `0 0 ${width} ${H}` }, wrap);
  const iw = width - labelW - noteW;
  rows.forEach((r, i) => {
    const y = 4 + i * rowH;
    el("text", { x: labelW - 8, y: y + rowH / 2 + 3, "text-anchor": "end", class: "tick" }, svg).textContent = r.label;
    const total = Object.values(r.parts).reduce((a, b) => a + b, 0) || 1;
    let x = labelW;
    for (const [name, v] of Object.entries(r.parts)) {
      const w = (v / total) * iw;
      const rect = el("rect", { x, y: y + 3, width: Math.max(w, 0), height: rowH - 8, fill: cssVar(colors[name] || "var(--ink-3)") }, svg);
      el("title", {}, rect).textContent = `${name}: ${(100 * v / total).toFixed(1)}%`;
      x += w;
    }
    el("text", { x: labelW + iw + 8, y: y + rowH / 2 + 3, class: "tick" }, svg).textContent = r.note || "";
  });
  const lg = document.createElement("div");
  lg.className = "legend";
  for (const [name, c] of Object.entries(colors)) {
    const sp = document.createElement("span");
    sp.style.color = cssVar(c);
    sp.innerHTML = `<i class="dot"></i><span style="color: var(--ink-2)"></span>`;
    sp.lastChild.textContent = name;
    lg.appendChild(sp);
  }
  container.appendChild(lg);
}
