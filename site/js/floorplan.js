// 2D Lattice Surgery Floorplan Visualizer
// Interactive visualizer for Litinski's Fast Block 2D tile layout,
// multi-patch boundary merges, routing corridors, and distillation factories.

import { cssVar } from "./plot.js";

export function initFloorplan() {
  const container = document.getElementById("floorplan-canvas");
  if (!container) return;

  const state = {
    mode: "cnot", // 'cnot', 'mbu', 'distill'
    step: 0,      // 0 to 4
    playing: false,
    interval: null,
    d: 13,
  };

  const MODES = {
    cnot: {
      name: "Lattice Surgery CNOT (2d rounds)",
      desc: "Transversal CNOT between Q₀ (Control) and Q₁ (Target) mediated by routing corridor patch A₀ via sequential Z⊗Z and X⊗X boundary merges.",
      steps: [
        { label: "Step 0: Initial State", sub: "Control Q₀ and Target Q₁ separated by idle routing corridor A₀." },
        { label: "Step 1: Z⊗Z Merge (d rounds)", sub: "Control Z-boundary and Ancilla merge along horizontal corridor; joint parity extracted." },
        { label: "Step 2: Z-Split & Intermediate State", sub: "Corridor splits; measurement outcome recorded for Pauli-frame classical correction." },
        { label: "Step 3: X⊗X Merge (d rounds)", sub: "Target X-boundary merges with Ancilla across vertical corridor; transversal parity check." },
        { label: "Step 4: Final Split & Disentanglement", sub: "Lattice surgery CNOT complete in 2d rounds. Both patches restored to independent routing access." },
      ],
    },
    mbu: {
      name: "Measurement-Based Uncomputation (MBU Carry)",
      desc: "Gidney carry uncomputation: eliminates reverse adder Toffolis using mid-circuit X-basis measurement and dynamic classical CZ feed-forward.",
      steps: [
        { label: "Step 0: Carry Generation", sub: "Carry bit generated into temporary ancilla tile C₀ during forward addition." },
        { label: "Step 1: Computation Complete", sub: "Data sum bits computed; carry tile C₀ holds entangled carry register." },
        { label: "Step 2: Mid-Circuit X Measurement", sub: "Transversal Hadamard + X readout collapses C₀, extracting parity in 1 round." },
        { label: "Step 3: Feed-Forward CZ Fixup", sub: "Classical decoder triggers Pauli-frame CZ correction on operands if measurement outcome is -1." },
        { label: "Step 4: Immediate Patch Reset", sub: "Tile C₀ immediately resets to |0⟩ for next digit, saving 2d idle rounds per adder." },
      ],
    },
    distill: {
      name: "Magic State Cultivation & CCZ Injection",
      desc: "Factory block produces high-fidelity |CCZ⟩ states and routes them via fast-block corridors to data register for parallel multi-qubit Toffoli merges.",
      steps: [
        { label: "Step 0: Factory Seed Injection", sub: "Injection patch receives noisy physical magic state; 1st level distillation initiated." },
        { label: "Step 1: Cultivation Rounds", sub: "Error syndrome checks cycle on 3×4 hot patches; verified |CCZ⟩ state prepared." },
        { label: "Step 2: Corridor Routing", sub: "Corridor tiles A₀, A₁ clear routing path from factory buffer to target patches Q₀, Q₁, Q₂." },
        { label: "Step 3: 3-Way Z⊗Z Merge", sub: "Simultaneous 3-body parity merge between data patches and |CCZ⟩ resource patch." },
        { label: "Step 4: State Consumption & Cleared Buffer", sub: "Magic state consumed; Toffoli phase kicked back; factory buffer resets for next cycle." },
      ],
    },
  };

  function render() {
    const W = 720;
    const H = 340;
    const cur = MODES[state.mode];
    const s = state.step;

    // Grid layout: 4 columns x 3 rows of tiles + Factory Block on the right
    // Tile size: 84 x 74 with 14px gaps
    const tw = 84, th = 74, gap = 16, ox = 32, oy = 40;

    // Tiles definition
    // (col, row, type, label, role)
    const tiles = [
      // Row 0
      { c: 0, r: 0, type: "data", id: "Q0", name: "Q₀ (Data/Ctrl)", active: (state.mode === "cnot" && (s === 1 || s === 2)) || (state.mode === "distill" && s === 3) },
      { c: 1, r: 0, type: "corridor", id: "A0", name: "A₀ (Routing)", active: (state.mode === "cnot" && (s === 1 || s === 3)) || (state.mode === "mbu" && s === 2) || (state.mode === "distill" && (s === 2 || s === 3)) },
      { c: 2, r: 0, type: "data", id: "Q1", name: "Q₁ (Data/Tgt)", active: (state.mode === "cnot" && (s === 3 || s === 4)) || (state.mode === "distill" && s === 3) },
      { c: 3, r: 0, type: "corridor", id: "A1", name: "A₁ (Routing)", active: state.mode === "distill" && s === 2 },

      // Row 1
      { c: 0, r: 1, type: "corridor", id: "A2", name: "A₂ (Routing)", active: state.mode === "cnot" && s === 3 },
      { c: 1, r: 1, type: "ancilla", id: "C0", name: "C₀ (MBU Carry)", active: state.mode === "mbu" && (s === 1 || s === 2 || s === 3) },
      { c: 2, r: 1, type: "corridor", id: "A3", name: "A₃ (Routing)", active: false },
      { c: 3, r: 1, type: "data", id: "Q2", name: "Q₂ (Data)", active: state.mode === "distill" && s === 3 },

      // Row 2
      { c: 0, r: 2, type: "data", id: "Q3", name: "Q₃ (Data)", active: false },
      { c: 1, r: 2, type: "corridor", id: "A4", name: "A₄ (Routing)", active: false },
      { c: 2, r: 2, type: "data", id: "Q4", name: "Q₄ (Data)", active: false },
      { c: 3, r: 2, type: "corridor", id: "A5", name: "A₅ (Routing)", active: false },
    ];

    // Factory Block occupies col 4, spans rows 0..2
    const factoryActive = state.mode === "distill" && (s === 0 || s === 1 || s === 2);
    const factoryX = ox + 4 * (tw + gap) + 12;
    const factoryY = oy;
    const factoryW = 160;
    const factoryH = 3 * th + 2 * gap;

    let svg = `<svg viewBox="0 0 ${W} ${H}" style="width:100%;display:block;border-radius:6px;background:var(--bg-panel,#fff);border:1px solid var(--rule-soft,#e5e7eb);">`;

    // Background & grid guides
    svg += `<defs>
      <pattern id="fp-dots" width="16" height="16" patternUnits="userSpaceOnUse">
        <circle cx="2" cy="2" r="0.8" fill="${cssVar("var(--rule-soft)")}"/>
      </pattern>
    </defs>
    <rect width="${W}" height="${H}" fill="url(#fp-dots)"/>`;

    // Draw connection channels / boundary merges
    if (state.mode === "cnot") {
      if (s === 1) {
        // Horizontal merge between Q0 and A0
        const x1 = ox + tw, y1 = oy + th / 2;
        const x2 = ox + tw + gap, y2 = y1;
        svg += `<line x1="${x1 - 4}" y1="${y1}" x2="${x2 + 4}" y2="${y2}" stroke="${cssVar("var(--accent)")}" stroke-width="8" stroke-linecap="round"/>`;
        svg += `<text x="${(x1 + x2) / 2}" y="${y1 - 10}" text-anchor="middle" style="font:600 11px var(--mono);fill:${cssVar("var(--accent)")}">Z ⊗ Z MERGE</text>`;
      } else if (s === 3) {
        // Merge between A0 and Q1
        const x1 = ox + tw + gap + tw, y1 = oy + th / 2;
        const x2 = x1 + gap, y2 = y1;
        svg += `<line x1="${x1 - 4}" y1="${y1}" x2="${x2 + 4}" y2="${y2}" stroke="${cssVar("var(--s1)")}" stroke-width="8" stroke-linecap="round"/>`;
        svg += `<text x="${(x1 + x2) / 2}" y="${y1 - 10}" text-anchor="middle" style="font:600 11px var(--mono);fill:${cssVar("var(--s1)")}">X ⊗ X MERGE</text>`;
      }
    } else if (state.mode === "mbu") {
      if (s === 2) {
        // MBU mid-circuit measurement beam on C0
        const cx = ox + 1 * (tw + gap) + tw / 2;
        const cy = oy + 1 * (th + gap) + th / 2;
        svg += `<circle cx="${cx}" cy="${cy}" r="32" fill="none" stroke="${cssVar("var(--accent)")}" stroke-width="3" stroke-dasharray="4 3"/>`;
        svg += `<text x="${cx}" y="${cy - 38}" text-anchor="middle" style="font:600 11px var(--mono);fill:${cssVar("var(--accent)")}">M_X READOUT</text>`;
      } else if (s === 3) {
        // Feed-forward arrow to Q0 and Q1
        const x1 = ox + 1 * (tw + gap) + tw / 2;
        const y1 = oy + 1 * (th + gap);
        const x2 = ox + tw / 2;
        const y2 = oy + th;
        svg += `<path d="M ${x1} ${y1} Q ${x1 - 20} ${y1 - 20} ${x2 + 20} ${y2}" fill="none" stroke="${cssVar("var(--ok)")}" stroke-width="2.5" stroke-dasharray="4 2"/>`;
        svg += `<text x="${(x1 + x2) / 2}" y="${(y1 + y2) / 2 - 8}" text-anchor="middle" style="font:600 10.5px var(--mono);fill:${cssVar("var(--ok)")}">FEED-FORWARD CZ</text>`;
      }
    } else if (state.mode === "distill") {
      if (s === 2 || s === 3) {
        // Magic state transport channel
        const fx = factoryX;
        const fy = factoryY + th / 2;
        const tx = ox + 3 * (tw + gap) + tw;
        svg += `<line x1="${fx}" y1="${fy}" x2="${tx}" y2="${fy}" stroke="${cssVar("var(--s4)")}" stroke-width="6" stroke-dasharray="6 3"/>`;
        svg += `<text x="${(fx + tx) / 2}" y="${fy - 10}" text-anchor="middle" style="font:600 11px var(--mono);fill:${cssVar("var(--s4)")}">|CCZ⟩ ROUTE</text>`;
      }
    }

    // Render individual tiles
    tiles.forEach((t) => {
      const x = ox + t.c * (tw + gap);
      const y = oy + t.r * (th + gap);
      const isData = t.type === "data";
      const isCarry = t.type === "ancilla";

      let fill = isData ? "var(--bg,#f8fafc)" : (isCarry ? "var(--bg-soft,#f1f5f9)" : "rgba(226,232,240,0.4)");
      let stroke = isData ? "var(--rule-dark,#94a3b8)" : "var(--rule-soft,#cbd5e1)";
      let strokeWidth = 1.2;

      if (t.active) {
        fill = isData ? "rgba(239,68,68,0.12)" : (isCarry ? "rgba(16,185,129,0.15)" : "rgba(59,130,246,0.15)");
        stroke = isData ? "var(--accent)" : (isCarry ? "var(--ok)" : "var(--s1)");
        strokeWidth = 2.4;
      }

      svg += `<g class="fp-tile" data-id="${t.id}" data-name="${t.name}" style="cursor:pointer">`;
      svg += `<rect x="${x}" y="${y}" width="${tw}" height="${th}" rx="5" fill="${cssVar(fill)}" stroke="${cssVar(stroke)}" stroke-width="${strokeWidth}"/>`;

      // Draw boundary polarity indicators for fast block data patches
      if (isData) {
        // Top boundary: Z (blue)
        svg += `<line x1="${x + 4}" y1="${y + 2}" x2="${x + tw - 4}" y2="${y + 2}" stroke="${cssVar("var(--s1)")}" stroke-width="3.5" stroke-linecap="round"/>`;
        // Right boundary: X (coral/red)
        svg += `<line x1="${x + tw - 2}" y1="${y + 4}" x2="${x + tw - 2}" y2="${y + th - 4}" stroke="${cssVar("var(--accent)")}" stroke-width="3.5" stroke-linecap="round"/>`;
      }

      // Tile Labels
      svg += `<text x="${x + tw / 2}" y="${y + th / 2 - 4}" text-anchor="middle" style="font:600 12px var(--mono);fill:${cssVar(t.active ? "var(--ink)" : "var(--ink-2)")}">${t.id}</text>`;
      svg += `<text x="${x + tw / 2}" y="${y + th / 2 + 12}" text-anchor="middle" style="font:400 9.5px var(--sans);fill:${cssVar("var(--ink-3)")}">${t.type}</text>`;
      svg += `</g>`;
    });

    // Render Factory Block
    let fFill = factoryActive ? "rgba(245,158,11,0.12)" : "var(--bg-soft,#f8fafc)";
    let fStroke = factoryActive ? "var(--s4)" : "var(--rule-soft,#cbd5e1)";
    svg += `<g class="fp-factory" style="cursor:pointer">`;
    svg += `<rect x="${factoryX}" y="${factoryY}" width="${factoryW}" height="${factoryH}" rx="6" fill="${cssVar(fFill)}" stroke="${cssVar(fStroke)}" stroke-width="${factoryActive ? 2.2 : 1.2}" stroke-dasharray="${factoryActive ? 'none' : '4 2'}"/>`;
    svg += `<text x="${factoryX + factoryW / 2}" y="${factoryY + 28}" text-anchor="middle" style="font:600 13px var(--sans);fill:${cssVar("var(--ink)")}">Distillation Factory</text>`;
    svg += `<text x="${factoryX + factoryW / 2}" y="${factoryY + 46}" text-anchor="middle" style="font:400 10.5px var(--mono);fill:${cssVar("var(--ink-2)")}">3×4 Hot Patches</text>`;

    // Internal factory stage status
    const facStatus = state.mode === "distill"
      ? (s === 0 ? "Injection Phase" : (s === 1 ? "Distilling (|CCZ⟩)" : (s === 2 ? "Buffering & Output" : "Re-initializing")))
      : "Idle (Standby)";
    svg += `<rect x="${factoryX + 16}" y="${factoryY + 68}" width="${factoryW - 32}" height="32" rx="4" fill="${cssVar("var(--bg-panel)")}" stroke="${cssVar("var(--rule-soft)")}"/>`;
    svg += `<text x="${factoryX + factoryW / 2}" y="${factoryY + 88}" text-anchor="middle" style="font:500 10px var(--mono);fill:${cssVar(factoryActive ? "var(--s4)" : "var(--ink-3)")}">${facStatus}</text>`;

    svg += `<text x="${factoryX + factoryW / 2}" y="${factoryY + factoryH - 18}" text-anchor="middle" style="font:400 9.5px var(--sans);fill:${cssVar("var(--ink-3)")}">Footprint: ~47k qubits</text>`;
    svg += `</g>`;

    // Legend at bottom
    const legY = H - 16;
    svg += `<g style="font:400 10.5px var(--sans);fill:${cssVar("var(--ink-2)")}">`;
    svg += `<line x1="${ox}" y1="${legY - 4}" x2="${ox + 16}" y2="${legY - 4}" stroke="${cssVar("var(--s1)")}" stroke-width="3"/>`;
    svg += `<text x="${ox + 22}" y="${legY}">Z-Boundary</text>`;
    svg += `<line x1="${ox + 100}" y1="${legY - 4}" x2="${ox + 116}" y2="${legY - 4}" stroke="${cssVar("var(--accent)")}" stroke-width="3"/>`;
    svg += `<text x="${ox + 122}" y="${legY}">X-Boundary</text>`;
    svg += `<rect x="${ox + 200}" y="${legY - 9}" width="10" height="10" fill="none" stroke="${cssVar("var(--rule-soft)")}"/>`;
    svg += `<text x="${ox + 216}" y="${legY}">Corridor Ancilla</text>`;
    svg += `<rect x="${ox + 320}" y="${legY - 9}" width="10" height="10" fill="${cssVar("var(--s4)")}" opacity="0.4"/>`;
    svg += `<text x="${ox + 336}" y="${legY}">Factory</text>`;
    svg += `</g>`;

    svg += `</svg>`;
    container.innerHTML = svg;

    // Update explanation texts
    const titleEl = document.getElementById("floorplan-step-title");
    const subEl = document.getElementById("floorplan-step-desc");
    if (titleEl && cur.steps[s]) titleEl.textContent = cur.steps[s].label;
    if (subEl && cur.steps[s]) subEl.textContent = cur.steps[s].sub;

    const scrub = document.getElementById("floorplan-scrub");
    if (scrub) scrub.value = s;
  }

  // Attach controls
  const modeSel = document.getElementById("floorplan-mode");
  if (modeSel) {
    modeSel.value = state.mode;
    modeSel.onchange = () => {
      state.mode = modeSel.value;
      state.step = 0;
      render();
    };
  }

  const scrub = document.getElementById("floorplan-scrub");
  if (scrub) {
    scrub.oninput = () => {
      state.step = parseInt(scrub.value, 10);
      render();
    };
  }

  const playBtn = document.getElementById("floorplan-play");
  if (playBtn) {
    playBtn.onclick = () => {
      state.playing = !state.playing;
      playBtn.textContent = state.playing ? "Pause" : "Play";
      if (state.playing) {
        state.interval = setInterval(() => {
          state.step = (state.step + 1) % 5;
          render();
        }, 1200);
      } else {
        clearInterval(state.interval);
      }
    };
  }

  // Initial draw
  render();
}
