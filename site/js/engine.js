// ftsim in the browser: the Rust engine compiled to WebAssembly (no bindings generator; C-ABI
// exports from src/wasm.rs). Works in browsers, workers and Node.

export async function loadEngine(source) {
  let bytes;
  if (source instanceof ArrayBuffer || ArrayBuffer.isView(source)) bytes = source;
  else bytes = await (await fetch(source)).arrayBuffer();
  const { instance } = await WebAssembly.instantiate(bytes, {});
  return new Engine(instance.exports);
}

export class Engine {
  constructor(x) {
    this.x = x;
  }

  error() {
    const x = this.x;
    const p = x.ft_error_ptr(), n = x.ft_error_len();
    return new TextDecoder().decode(new Uint8Array(x.memory.buffer, p, n));
  }

  /** Parses program text; throws with the engine's message on bad input. */
  program(text) {
    const x = this.x;
    const bytes = new TextEncoder().encode(text);
    const ptr = x.ft_alloc(bytes.length);
    new Uint8Array(x.memory.buffer, ptr, bytes.length).set(bytes);
    const h = x.ft_program_new(ptr, bytes.length);
    x.ft_free(ptr, bytes.length);
    if (!h) throw new Error(this.error());
    return new WasmProgram(this, h);
  }
}

export class WasmProgram {
  constructor(engine, h) {
    this.e = engine;
    this.h = h;
    this.numMeasurements = engine.x.ft_num_measurements(h);
    this.expectedFaults = engine.x.ft_expected_faults(h);
  }

  faultDistribution(kmax) {
    const x = this.e.x;
    const p = x.ft_fault_distribution(this.h, kmax);
    return Array.from(new Float64Array(x.memory.buffer, p, kmax + 2));
  }

  /** Records of shots firstShot…firstShot+shots−1: an array of Uint8Array rows (bit-packed).
   *  faults: "plain" | "none" | k. */
  sample(shots, { seed, firstShot = 0, faults = "plain" }) {
    const x = this.e.x;
    const mode = faults === "plain" ? -1 : faults === "none" ? -2 : faults;
    const p = x.ft_sample(this.h, shots, seed, firstShot, mode);
    if (!p) throw new Error(this.e.error());
    const bps = Math.ceil(this.numMeasurements / 8);
    const all = new Uint8Array(x.memory.buffer, p, shots * bps).slice();
    const rows = [];
    for (let s = 0; s < shots; s++) rows.push(all.subarray(s * bps, (s + 1) * bps));
    return rows;
  }

  free() {
    this.e.x.ft_program_free(this.h);
  }
}

export function bit(row, j) {
  return (row[j >> 3] >> (j & 7)) & 1;
}
