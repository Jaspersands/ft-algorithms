// Resource calculator. Every row comes from the estimator of tools/scale.py, precomputed by
// tools/site_data.py over the measured noise strengths, both factories and two failure budgets;
// the page looks rows up and scales the run time with the cycle time. Nothing is modelled here.

const FAC_LABEL = { cultivation: "magic-state cultivation", "15to1": "15-to-1 distillation" };
const BUDGET_COLORS = { idle: "var(--s7)", CCZ: "var(--s2)", CNOT: "var(--s1)", T: "var(--s3)", S: "var(--s5)" };

function fmtQ(x) {
  return x >= 1e6 ? `${(x / 1e6).toPrecision(3)} M` : `${(x / 1e3).toPrecision(3)} k`;
}

function fmtTime(s) {
  if (s < 120) return `${s.toPrecision(3)} s`;
  if (s < 7200) return `${(s / 60).toPrecision(3)} min`;
  if (s < 2 * 86400) return `${(s / 3600).toPrecision(3)} h`;
  if (s < 400 * 86400) return `${(s / 86400).toPrecision(3)} days`;
  return `${(s / 86400 / 365.25).toPrecision(3)} years`;
}

function fmtSci(x) {
  const e = Math.floor(Math.log10(x));
  return `${(x / 10 ** e).toFixed(2)}×10<sup>${e}</sup>`;
}

export function initCalculator(grid) {
  const out = document.getElementById("calc-output");
  if (!out || !grid) return;
  const state = { w: "textbook", cal: "separate", p: 0.001, factory: "cultivation", target: 0.1, tau: 1.0 };

  const sel = (id, options, key, parse = (v) => v) => {
    const el = document.getElementById(id);
    if (!el) return;
    el.innerHTML = options.map(([v, label]) => `<option value="${v}">${label}</option>`).join("");
    el.value = String(state[key]);
    el.onchange = () => {
      state[key] = parse(el.value);
      render();
    };
  };
  sel("calc-preset", Object.entries(grid.workloads).map(([k, w]) => [k, w.name]), "w");
  const CAL_LABEL = { separate: "separate measure and reset layers (as calibrated)", mr: "one measure-and-reset layer (recalibrated)" };
  sel("calc-cal", (grid.calibrations || ["separate"]).map((c) => [c, CAL_LABEL[c] || c]), "cal");
  const calRow = document.getElementById("calc-cal-row");
  if (calRow && !(grid.calibrations && grid.calibrations.length > 1)) calRow.style.display = "none";
  sel("calc-p", [[0.001, "0.1%"], [0.002, "0.2%"]], "p", Number);
  sel("calc-fac", Object.entries(FAC_LABEL), "factory");
  sel("calc-target", [[0.1, "90% (E[faults] ≤ 0.1)"], [0.01, "99% (E[faults] ≤ 0.01)"]], "target", Number);
  const tau = document.getElementById("calc-tau"), tauVal = document.getElementById("calc-tau-val");
  if (tau) {
    tau.oninput = () => {
      state.tau = Number(tau.value);
      if (tauVal) tauVal.textContent = `${state.tau.toFixed(1)} µs`;
      render();
    };
  }

  function render() {
    const w = grid.workloads[state.w];
    const r = grid.rows.find((x) => x.w === state.w && (x.cal || "separate") === state.cal && x.p === state.p && x.factory === state.factory && x.target === state.target);
    if (!r) return;
    const seconds = r.rounds * state.tau * 1e-6;
    const card = (label, val, sub, hl = false) =>
      `<div class="calc-card${hl ? " highlight" : ""}"><div class="card-label">${label}</div><div class="card-val">${val}</div><div class="card-sub">${sub}</div></div>`;
    const pubNote = w.published_d
      ? `published: d = ${w.published_d} at ${(100 * w.published_success).toFixed(1)}% success`
      : "this work's circuit";
    let html = `<div class="calc-results-grid">`;
    if (r.d === null) {
      html += card("Code distance", "none", r.floor > state.target
        ? `magic-state errors alone give ${r.floor.toPrecision(2)} expected faults`
        : "the target needs d > 61", true);
    } else {
      html += card("Code distance", `d = ${r.d}`, pubNote, true);
      html += card("Physical qubits", r.physical_qubits ? fmtQ(r.physical_qubits) : "—",
        r.physical_qubits ? `${fmtQ(r.factory_qubits)} in factories` : "authors' layout not modelled");
      html += card("Run time", fmtTime(seconds), `${fmtSci(r.rounds)} rounds at ${state.tau.toFixed(1)} µs`);
      html += card("Expected faults", r.expected_faults.toPrecision(2), r.extrapolated ? "channels extrapolated beyond the measured distances" : "within the measured range");
    }
    html += `</div>`;
    if (r.d !== null) {
      const total = Object.values(r.breakdown).reduce((a, b) => a + b, 0) || 1;
      const parts = Object.entries(r.breakdown).filter(([, v]) => v / total > 0.005);
      html += `<div class="calc-breakdown"><h4>Where the expected faults come from</h4>
        <div style="display:flex;height:22px;border-radius:4px;overflow:hidden;margin:8px 0;background:var(--rule-soft)">${parts
          .map(([k, v]) => `<div title="${k}: ${((100 * v) / total).toFixed(1)}%" style="width:${(100 * v) / total}%;background:${BUDGET_COLORS[k] || "var(--ink-3)"}"></div>`)
          .join("")}</div>
        <div style="display:flex;gap:16px;flex-wrap:wrap;font:400 11px var(--mono);color:var(--ink-2)">${parts
          .map(([k, v]) => `<span><span style="color:${BUDGET_COLORS[k] || "var(--ink-3)"}">■</span> ${k} ${((100 * v) / total).toFixed(1)}%</span>`)
          .join("")}</div></div>`;
    }
    out.innerHTML = html;
  }

  render();
}
