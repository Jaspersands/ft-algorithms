// Interactive Resource Trade-off Calculator
// Real-time Pareto exploration across workloads, hardware error rates,
// cycle times, surface code architectures (standard vs XZZX), and distillation factories.

import { cssVar } from "./plot.js";

export function initCalculator() {
  const container = document.getElementById("calc-container");
  if (!container) return;

  const PRESETS = {
    rsa2048_mod: {
      name: "RSA-2048 (Modern Windowed + MBU)",
      qubits: 8199,
      toffolis: 3.5e9,
      roundsPerToffoli: 2.1,
      desc: "Modern windowed arithmetic (k=2) with measurement-based carry uncomputation.",
    },
    rsa2048_tb: {
      name: "RSA-2048 (Textbook Cuccaro)",
      qubits: 6150,
      toffolis: 8.8e9,
      roundsPerToffoli: 4.8,
      desc: "Textbook reversible Cuccaro adders and sequential modular exponentiation.",
    },
    femoco: {
      name: "FeMoco Nitrogenase (Lee et al. 2021)",
      qubits: 2142,
      toffolis: 5.3e9,
      roundsPerToffoli: 1.0,
      desc: "Tensor hypercontraction (THC) active space simulation of nitrogenase cofactor.",
    },
    lih_qpe: {
      name: "LiH Chemistry (10-bit QPE)",
      qubits: 3,
      toffolis: 0,
      t_gates: 1.8e7,
      roundsPerToffoli: 1.0,
      fixedRounds: 2.5e7,
      desc: "Active space ground state phase estimation with 10-bit phase resolution.",
    },
    h2_qpe: {
      name: "H₂ Chemistry (10-bit QPE)",
      qubits: 3,
      toffolis: 0,
      t_gates: 1.8e7,
      roundsPerToffoli: 1.0,
      fixedRounds: 1.5e7,
      desc: "Full Trotterized iterative phase estimation targeting 1.6 mHa chemical accuracy.",
    },
  };

  const state = {
    preset: "rsa2048_mod",
    p: 0.001,
    tau_us: 1.0,
    arch: "standard", // 'standard', 'xzzx'
    eta: 100.0,
    factory: "cultivation", // 'cultivation', '15to1'
    targetFaults: 0.01,
  };

  // Calibration constants
  const P_STAR = 0.008525875;
  const A_PREFACTOR = 0.06165385;

  function calculate() {
    const w = PRESETS[state.preset];
    const n_logical = w.qubits;

    // Fast-block tiles = 2n + ceil(sqrt(8n)) + 1
    const tiles = 2 * n_logical + Math.ceil(Math.sqrt(8 * n_logical)) + 1;

    // Approximate total rounds
    const rounds = w.fixedRounds
      ? w.fixedRounds
      : (w.toffolis * w.roundsPerToffoli * 25); // estimate rounds scaling with distance

    const spaceTimeRounds = n_logical * rounds;
    const target_p_per_round = state.targetFaults / Math.max(spaceTimeRounds, 1);

    // Compute required distances
    let dx = 3, dz = 3;
    const p_tot = state.p;

    if (state.arch === "xzzx") {
      const pz = (state.eta * p_tot) / (1.0 + state.eta);
      const px = p_tot / (2.0 * (1.0 + state.eta));

      // dx for Z errors
      const rz = pz / P_STAR;
      if (rz < 1.0) {
        const req_hx = Math.log((target_p_per_round / 2.0) / A_PREFACTOR) / Math.log(rz);
        dx = Math.max(3, Math.ceil(2.0 * req_hx - 1.0));
        if (dx % 2 === 0) dx++;
      } else {
        dx = 55;
      }

      // dz for X errors
      const rx = px / P_STAR;
      const req_hz = Math.log((target_p_per_round / 2.0) / A_PREFACTOR) / Math.log(rx);
      dz = Math.max(3, Math.ceil(2.0 * req_hz - 1.0));
      if (dz % 2 === 0) dz++;

      dx = Math.min(dx, 55);
      dz = Math.min(dz, 55);
    } else {
      // Standard symmetric code: d_sym protects against worst-case error
      const r = p_tot / P_STAR;
      if (r < 1.0) {
        const req_h = Math.log(target_p_per_round / A_PREFACTOR) / Math.log(r);
        dx = Math.max(3, Math.ceil(2.0 * req_h - 1.0));
        if (dx % 2 === 0) dx++;
      } else {
        dx = 55;
      }
      dz = dx;
    }

    // Footprints
    const tile_qubits = 2 * (dx + 1) * (dz + 1);
    const data_qubits = tiles * tile_qubits;

    // Factory footprint
    let fac_qubits = 0;
    if (state.factory === "cultivation") {
      // Gidney 2025: ~6 factory blocks of 1,352 qubits = 8,112 qubits base
      fac_qubits = 8112 * Math.ceil((dx + 1) / 28);
    } else {
      // 15-to-1 two-level: ~47,000 qubits
      fac_qubits = 47000 * Math.ceil((dx + 1) / 28);
    }

    const total_phys = data_qubits + fac_qubits;

    // Runtime calculation
    const total_seconds = rounds * (state.tau_us * 1e-6);

    // Baseline symmetric comparison (eta = 0.5 / symmetric standard code)
    const sym_d = Math.max(dx, dz);
    const sym_tile = 2 * (sym_d + 1) * (sym_d + 1);
    const sym_total = tiles * sym_tile + fac_qubits;
    const savingsRatio = sym_total / total_phys;

    return {
      dx,
      dz,
      tiles,
      tile_qubits,
      data_qubits,
      fac_qubits,
      total_phys,
      total_seconds,
      sym_total,
      savingsRatio,
      rounds,
    };
  }

  function fmtTime(sec) {
    if (sec < 60) return `${sec.toFixed(1)} s`;
    if (sec < 3600) return `${(sec / 60).toFixed(1)} min`;
    if (sec < 86400) return `${(sec / 3600).toFixed(1)} hours`;
    if (sec < 86400 * 365) return `${(sec / 86400).toFixed(1)} days`;
    return `${(sec / (86400 * 365)).toFixed(1)} years`;
  }

  function render() {
    const res = calculate();
    const w = PRESETS[state.preset];

    const outEl = document.getElementById("calc-output");
    if (!outEl) return;

    outEl.innerHTML = `
      <div class="calc-results-grid">
        <div class="calc-card highlight">
          <div class="card-label">Total Physical Qubits</div>
          <div class="card-val">${(res.total_phys / 1e6).toFixed(2)}M</div>
          <div class="card-sub">${res.total_phys.toLocaleString()} qubits</div>
        </div>
        <div class="calc-card">
          <div class="card-label">Execution Runtime</div>
          <div class="card-val">${fmtTime(res.total_seconds)}</div>
          <div class="card-sub">@ ${state.tau_us} µs/cycle</div>
        </div>
        <div class="calc-card">
          <div class="card-label">Code Distance</div>
          <div class="card-val">${state.arch === "xzzx" ? `d_X=${res.dx}, d_Z=${res.dz}` : `d = ${res.dx}`}</div>
          <div class="card-sub">${res.tile_qubits} qubits/tile</div>
        </div>
        <div class="calc-card">
          <div class="card-label">Footprint Savings</div>
          <div class="card-val">${res.savingsRatio > 1.05 ? `${res.savingsRatio.toFixed(2)}×` : "1.0×"}</div>
          <div class="card-sub">${state.arch === "xzzx" ? `vs symmetric d=${Math.max(res.dx, res.dz)}` : "symmetric baseline"}</div>
        </div>
      </div>

      <div class="calc-breakdown">
        <h4>Physical Resource Allocation</h4>
        <div class="bar-container" style="display:flex;height:24px;border-radius:4px;overflow:hidden;margin:8px 0;background:var(--rule-soft)">
          <div style="width:${(res.data_qubits / res.total_phys) * 100}%;background:var(--accent);display:flex;align-items:center;justify-content:center;color:#fff;font:600 11px var(--mono)">Data (${((res.data_qubits / res.total_phys) * 100).toFixed(0)}%)</div>
          <div style="width:${(res.fac_qubits / res.total_phys) * 100}%;background:var(--s4);display:flex;align-items:center;justify-content:center;color:#fff;font:600 11px var(--mono)">Factory (${((res.fac_qubits / res.total_phys) * 100).toFixed(0)}%)</div>
        </div>
        <div style="display:flex;justify-content:space-between;font:400 11px var(--mono);color:var(--ink-2)">
          <span>Data Patches + Corridors: ${res.data_qubits.toLocaleString()} (${res.tiles} tiles)</span>
          <span>Distillation Factory: ${res.fac_qubits.toLocaleString()}</span>
        </div>
      </div>
    `;
  }

  // Bind controls
  const presetSel = document.getElementById("calc-preset");
  if (presetSel) {
    presetSel.innerHTML = Object.entries(PRESETS)
      .map(([k, v]) => `<option value="${k}">${v.name}</option>`)
      .join("");
    presetSel.value = state.preset;
    presetSel.onchange = () => {
      state.preset = presetSel.value;
      render();
    };
  }

  const pSlider = document.getElementById("calc-p");
  const pVal = document.getElementById("calc-p-val");
  if (pSlider) {
    pSlider.oninput = () => {
      state.p = parseFloat(pSlider.value);
      if (pVal) pVal.textContent = `${(state.p * 100).toFixed(2)}%`;
      render();
    };
  }

  const tauSlider = document.getElementById("calc-tau");
  const tauVal = document.getElementById("calc-tau-val");
  if (tauSlider) {
    tauSlider.oninput = () => {
      state.tau_us = parseFloat(tauSlider.value);
      if (tauVal) tauVal.textContent = `${state.tau_us.toFixed(1)} µs`;
      render();
    };
  }

  const archSel = document.getElementById("calc-arch");
  const etaRow = document.getElementById("calc-eta-row");
  if (archSel) {
    archSel.value = state.arch;
    archSel.onchange = () => {
      state.arch = archSel.value;
      if (etaRow) etaRow.style.display = state.arch === "xzzx" ? "flex" : "none";
      render();
    };
  }

  const etaSlider = document.getElementById("calc-eta");
  const etaVal = document.getElementById("calc-eta-val");
  if (etaSlider) {
    etaSlider.oninput = () => {
      state.eta = parseFloat(etaSlider.value);
      if (etaVal) etaVal.textContent = `${state.eta.toFixed(0)}`;
      render();
    };
  }

  const facSel = document.getElementById("calc-fac");
  if (facSel) {
    facSel.value = state.factory;
    facSel.onchange = () => {
      state.factory = facSel.value;
      render();
    };
  }

  render();
}
