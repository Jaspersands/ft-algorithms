"""Circuit-level noise models applied layer by layer to a noiseless Stim circuit, and the
CNOT-to-CZ rewrite needed to run the surface code under SI1000.

Two models, both with a single strength p:

- ``sd6``: standard depolarizing, as stabilizer-qec and Stim's generated circuits use it:
  DEPOLARIZE2(p) after every two-qubit gate, DEPOLARIZE1(p) after every single-qubit gate and on
  every idle qubit of a layer, a flip with probability p before every measurement and after every
  reset. Applied to stabilizer-qec's noiseless memory circuit it reproduces that circuit's own
  error model exactly (tests/test_circuit_noise.py), which is what makes the comparison below fair.
- ``si1000``: Gidney, Newman and McEwen's superconducting-inspired model, as tabulated in Gidney,
  Newman, Brooks and Jones, "Yoked surface codes" (arXiv:2312.04522, Table 3): on circuits of
  single-qubit unitaries, CZ, Z-basis reset and Z-basis measurement, DEPOLARIZE2(p) after CZ,
  DEPOLARIZE1(p/10) after any single-qubit unitary and on idle qubits, X_ERROR(2p) after reset,
  a 5p flip of the result and DEPOLARIZE1(p) after measurement, and DEPOLARIZE1(2p) on qubits
  idle during a layer of resets or measurements.

A third, ``uniform``, is Stim's own generated rotated memory with every parameter at p
(``stim.Circuit.generated``: DEPOLARIZE1(p) on data qubits once per round, DEPOLARIZE2(p) after every
CX and DEPOLARIZE1(p) after every H, flips on reset and measurement), which unlike sd6 leaves qubits
that idle during a gate layer noiseless.

``to_cz`` turns a CX circuit into a CZ circuit (H on each CX target before and after its CZ, X-basis
resets and measurements into Z-basis ones with an H) and merges adjacent single-qubit layers.
"""

from __future__ import annotations

import stim

ANNOTATIONS = {"DETECTOR", "OBSERVABLE_INCLUDE", "QUBIT_COORDS", "SHIFT_COORDS"}
SINGLE = {"H", "I"}
TWO = {"CX", "CZ"}
MEASURE = {"M", "MX"}
RESET = {"R", "RX"}


def _moments(circuit: stim.Circuit) -> list[list[stim.CircuitInstruction]]:
    """The noiseless, flattened circuit as a list of moments (instructions between TICKs)."""
    out: list[list[stim.CircuitInstruction]] = [[]]
    for inst in circuit.without_noise().flattened():
        if inst.name == "TICK":
            out.append([])
        else:
            out[-1].append(inst)
    return [m for m in out if m]


def _qubits(inst: stim.CircuitInstruction) -> list[int]:
    return [t.value for t in inst.targets_copy() if t.is_qubit_target]


def _gates(moment) -> list[stim.CircuitInstruction]:
    return [i for i in moment if i.name not in ANNOTATIONS]


def _is_single(moment) -> bool:
    g = _gates(moment)
    return bool(g) and all(i.name in SINGLE for i in g) and len(g) == len(moment)


def to_cz(circuit: stim.Circuit) -> list[list[stim.CircuitInstruction]]:
    """Moments of the same computation over {H, CZ, R, M}."""
    moments: list[list[stim.CircuitInstruction]] = []
    for m in _moments(circuit):
        names = {i.name for i in _gates(m)}
        if names & TWO:
            if not names <= TWO:
                raise ValueError(f"two-qubit gates share a layer with {names - TWO}")
            targets = [q for i in m if i.name == "CX" for q in _qubits(i)[1::2]]
            pairs = [q for i in _gates(m) for q in _qubits(i)]
            if targets:
                moments.append([stim.CircuitInstruction("H", targets)])
            moments.append([stim.CircuitInstruction("CZ", pairs)] + [i for i in m if i.name in ANNOTATIONS])
            if targets:
                moments.append([stim.CircuitInstruction("H", targets)])
        elif names & {"RX"}:
            xq = [q for i in m if i.name == "RX" for q in _qubits(i)]
            zq = [q for i in m if i.name == "R" for q in _qubits(i)]
            moments.append([stim.CircuitInstruction("R", zq + xq)] + [i for i in m if i.name in ANNOTATIONS])
            moments.append([stim.CircuitInstruction("H", xq)])
        elif names & {"MX"}:
            xq = [q for i in m if i.name == "MX" for q in _qubits(i)]
            if names != {"MX"}:
                raise ValueError("MX shares a layer with other gates")
            moments.append([stim.CircuitInstruction("H", xq)])
            moments.append([stim.CircuitInstruction("M", xq)] + [i for i in m if i.name in ANNOTATIONS])
        else:
            moments.append(m)
    # Merge adjacent single-qubit layers (H·H = I) and drop the ones that become empty.
    merged: list[list[stim.CircuitInstruction]] = []
    for m in moments:
        if _is_single(m) and merged and _is_single(merged[-1]):
            odd: dict[int, None] = {}
            for i in merged[-1] + m:
                for q in _qubits(i):
                    if q in odd:
                        del odd[q]
                    else:
                        odd[q] = None
            merged[-1] = [stim.CircuitInstruction("H", list(odd))] if odd else []
            if not odd:
                merged.pop()
        else:
            merged.append(m)
    return merged


def merge_measure_reset(moments):
    """Each measurement layer followed by a reset of the same qubits becomes one measure-and-reset
    layer, as hardware and Gidney's generated circuits run it: qubits waiting for it idle once a
    round instead of twice."""
    out = []
    for m in moments:
        g = _gates(m)
        if out and g and all(i.name == "R" for i in g) and len(g) == len(m):
            prev = _gates(out[-1])
            if prev and all(i.name == "M" for i in prev) and \
                    sorted(q for i in prev for q in _qubits(i)) == sorted(q for i in g for q in _qubits(i)):
                ann = [i for i in out[-1] if i.name in ANNOTATIONS]
                out[-1] = prev + g + ann
                continue
        out.append(m)
    return out


def apply_noise(moments, model: str, p: float, num_qubits: int) -> stim.Circuit:
    """The moments with ``model``'s noise at strength p; every qubit 0..num_qubits-1 is live."""
    if model not in ("sd6", "si1000"):
        raise ValueError(model)
    si = model == "si1000"
    c = stim.Circuit()
    for k, m in enumerate(moments):
        gates = _gates(m)
        names = {i.name for i in gates}
        busy = {q for i in gates for q in _qubits(i)}
        idle = [q for q in range(num_qubits) if q not in busy]
        pre, post = stim.Circuit(), stim.Circuit()
        reset_here = {q for i in gates if i.name in RESET for q in _qubits(i)}
        for i in gates:
            qs = _qubits(i)
            if i.name in MEASURE:
                pre.append("Z_ERROR" if i.name == "MX" else "X_ERROR", qs, 5 * p if si else p)
                kept = [q for q in qs if q not in reset_here]  # a reset in the same layer supersedes it
                if si and kept:
                    post.append("DEPOLARIZE1", kept, p)
            elif i.name in RESET:
                post.append("Z_ERROR" if i.name == "RX" else "X_ERROR", qs, 2 * p if si else p)
            elif i.name in TWO:
                post.append("DEPOLARIZE2", qs, p)
            elif i.name in SINGLE:
                post.append("DEPOLARIZE1", qs, p / 10 if si else p)
            else:
                raise ValueError(f"no noise rule for {i.name}")
        if idle:
            if names & (MEASURE | RESET):
                q = 2 * p if si else p
            else:
                q = p / 10 if si else p
            # stabilizer-qec leaves qubits idle during the initial all-qubit reset noiseless; so do we.
            if not (k == 0 and names <= RESET):
                post.append("DEPOLARIZE1", idle, q)
        c += pre
        for i in gates:
            c.append(i)
        c += post
        for i in m:
            if i.name in ANNOTATIONS:
                c.append(i)
        if k + 1 < len(moments):
            c.append("TICK")
    return c


def memory(d: int, rounds: int, basis: str, model: str, p: float, merge_mr: bool | None = None) -> stim.Circuit:
    """A rotated memory with ``model``'s noise: stabilizer-qec's CX circuit for sd6 (its separate
    measurement and reset layers, as calibrated) and sd6_mr (one measure-and-reset layer), its CZ
    rewrite with one measure-and-reset layer for si1000, Stim's generated circuit for uniform."""
    if model == "sd6_mr":
        model, merge_mr = "sd6", True
    import stabilizer_qec as sq

    if model == "uniform":
        return stim.Circuit.generated(f"surface_code:rotated_memory_{basis}", distance=d, rounds=rounds,
                                      after_clifford_depolarization=p, before_round_data_depolarization=p,
                                      before_measure_flip_probability=p, after_reset_flip_probability=p)

    base = stim.Circuit(str(sq.memory_circuit(distance=d, rounds=rounds, p=p, basis=basis)))
    moments = to_cz(base) if model == "si1000" else _moments(base)
    if merge_mr if merge_mr is not None else model == "si1000":
        moments = merge_measure_reset(moments)
    return apply_noise(moments, model, p, base.num_qubits)


NOISE = {"X_ERROR", "Z_ERROR", "DEPOLARIZE1", "DEPOLARIZE2"}


def merge_measure_reset_layers(circuit: stim.Circuit) -> stim.Circuit:
    """A noisy circuit with each reset layer that follows a measurement layer folded into it: the
    resets and their flips join the measurement layer, and the reset layer's idle depolarizing is
    dropped, so qubits waiting for the two wait once. Every other noise channel is kept as it is,
    so this applies to any of stabilizer-qec's circuits (surgery and reference-qubit memories
    included, whose noiseless segments stay noiseless). On a memory it gives exactly sd6_mr."""
    moments: list[list[stim.CircuitInstruction]] = [[]]
    for inst in circuit.flattened():
        if inst.name == "TICK":
            moments.append([])
        else:
            moments[-1].append(inst)
    out: list[list[stim.CircuitInstruction]] = []
    for m in moments:
        gates = [i for i in m if i.name not in NOISE and i.name not in ANNOTATIONS]
        if out and gates and all(i.name in RESET for i in gates) and not any(i.name in ANNOTATIONS for i in m):
            prev = [i for i in out[-1] if i.name not in NOISE and i.name not in ANNOTATIONS]
            if prev and all(i.name in MEASURE for i in prev):
                out[-1] = out[-1] + [i for i in m if i.name != "DEPOLARIZE1"]
                continue
        out.append(m)
    c = stim.Circuit()
    for k, m in enumerate(out):
        for i in m:
            c.append(i)
        if k + 1 < len(out):
            c.append("TICK")
    return c
