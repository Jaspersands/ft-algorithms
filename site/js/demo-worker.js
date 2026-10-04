// Runs the demo's noisy programs off the main thread: fetches a gzipped program, parses it in
// the WebAssembly engine and samples shots, returning each shot's phase-estimate outcome y.
import { loadEngine, bit } from "./engine.js";

let engine = null;
let prog = null;
let records = [];

self.onmessage = async (ev) => {
  const msg = ev.data;
  try {
    if (msg.type === "init") {
      engine = await loadEngine(msg.wasm);
      records = msg.records;
      self.postMessage({ type: "ready" });
    } else if (msg.type === "load") {
      const res = await fetch(msg.url);
      const stream = res.body.pipeThrough(new DecompressionStream("gzip"));
      const text = await new Response(stream).text();
      if (prog) prog.free();
      const t0 = performance.now();
      prog = engine.program(text);
      self.postMessage({ type: "loaded", key: msg.key, parseMs: performance.now() - t0, expectedFaults: prog.expectedFaults });
    } else if (msg.type === "run") {
      const t0 = performance.now();
      const rows = prog.sample(msg.shots, { seed: msg.seed, firstShot: msg.firstShot });
      const ys = rows.map((r) => records.reduce((y, rec, j) => y | (bit(r, rec) << j), 0));
      self.postMessage({ type: "result", key: msg.key, ys, ms: performance.now() - t0, firstShot: msg.firstShot });
    }
  } catch (e) {
    self.postMessage({ type: "error", message: String(e && e.message ? e.message : e) });
  }
};
