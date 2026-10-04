// The WebAssembly engine gives the same records as the native one (same program, seed, shots).
import { readFileSync } from "node:fs";
import { execFileSync } from "node:child_process";
import { loadEngine, bit } from "../site/js/engine.js";

const text = "H 0 1 2\nREPEAT 30 {\nCCX 0 1 2\nT 0\nH 1\nDEPOLARIZE1(0.05) 0 1 2\nM 2\nTABLE rec[-1] {S 0} {H 0}\n}\nM 0 1\n";
const engine = await loadEngine(readFileSync(new URL("../site/wasm/ftsim.wasm", import.meta.url)));
const prog = engine.program(text);
const rows = prog.sample(200, { seed: 3 });
const ours = rows.map((r) => Array.from({ length: prog.numMeasurements }, (_, j) => bit(r, j)).join("")).join("\n");
const py = execFileSync(".venv/bin/python", ["-c", `
import ftalgo, sys
p = ftalgo.Program(sys.stdin.read())
r = p.sample(200, seed=3, backend="sparse", threads=1)
print("\\n".join("".join(str(int(b)) for b in row) for row in r))
`], { input: text }).toString().trim();
if (ours !== py) { console.error("MISMATCH"); process.exit(1); }
let bad = false;
try { engine.program("CX 0 0"); bad = true; } catch (e) { if (!String(e.message).includes("line 1")) bad = true; }
if (bad) { console.error("error path broken"); process.exit(1); }
console.log(`wasm == native on 200 shots × ${prog.numMeasurements} records; E[faults] ${prog.expectedFaults.toFixed(3)}`);
