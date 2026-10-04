"""A reference-qubit (Choi) memory experiment: the full logical Pauli channel of an idle patch.

Starting from stabilizer-qec's SD6 X-basis memory circuit (same noise and schedule as every
other calibration), the patch is projected noiselessly before any noise (MPP of every Z check,
which round 1's Z checks are then compared with) and a noiseless reference qubit R in |+⟩ is
entangled with it by a noiseless MPP Z_R·Z_L. Without the projection, an X error before the
first Z-check readout would flip Z_R Z_L unseen, since first-round Z outcomes are otherwise
random. The transversal X readout at the end is replaced by a noiseless perfect round (MPP of
every check, compared with the last noisy round) and noiseless MPPs of X_R·X_L and Z_R·Z_L.
Observable 0 (X_R X_L) flips on logical Z or Y, observable 1 (Z_R Z_L, with the initial MPP's
outcome) on logical X or Y: their joint distribution is the patch's channel, including the Y
component that separate X- and Z-basis memories cannot see.
"""

from __future__ import annotations

import re


def choi_memory(d: int, rounds: int, p: float) -> str:
    import stabilizer_qec as sq

    text = str(sq.memory_circuit(distance=d, rounds=rounds, p=p, basis="x"))
    lines = text.splitlines()
    coords: dict[int, tuple[int, int]] = {}
    for l in lines:
        m = re.match(r"QUBIT_COORDS\((\d+), (\d+)\) (\d+)", l)
        if m:
            coords[int(m.group(3))] = (int(m.group(1)), int(m.group(2)))
    final_mx = max(i for i, l in enumerate(lines) if l.startswith("MX "))
    data = [int(t) for t in lines[final_mx].split()[1:]]
    first_h = next(l for l in lines if l.startswith("H "))
    xchecks = [int(t) for t in first_h.split()[1:]]
    last_m = max(i for i, l in enumerate(lines) if l.startswith("M "))
    measured = [int(t) for t in lines[last_m].split()[1:]]
    pos = {coords[q]: q for q in data}

    def support(chk: int) -> list[int]:
        x, y = coords[chk]
        return [pos[(x + dx, y + dy)] for dx in (-1, 1) for dy in (-1, 1) if (x + dx, y + dy) in pos]

    xs = sorted({c[0] for c in (coords[q] for q in data)})
    ys = sorted({c[1] for c in (coords[q] for q in data)})
    x_logical = [pos[(x, ys[0])] for x in xs]  # a row: the circuit's own X observable
    z_logical = [pos[(xs[0], y)] for y in ys]  # a column: anticommutes with the row once
    ref = max(coords) + 1

    K = len(measured)
    zchecks = [c for c in measured if c not in xchecks]
    Kz = len(zchecks)
    out: list[str] = []
    rx_done = False
    first_m_done = False
    for i, l in enumerate(lines[:final_mx]):
        if l.startswith("Z_ERROR") and i == final_mx - 1:
            continue  # the readout's flip noise goes with the readout
        out.append(l)
        if not rx_done and l.startswith("RX "):
            out.insert(len(out) - 1, f"QUBIT_COORDS(-2, -2) {ref}")
            out.append(f"RX {ref}")
            for chk in zchecks:
                out.append("MPP " + "*".join(f"Z{q}" for q in support(chk)))
            out.append("MPP " + "*".join([f"Z{ref}"] + [f"Z{q}" for q in z_logical]))
            rx_done = True
        elif rx_done and not first_m_done and l.startswith("M "):
            first_m_done = True
            now = Kz + 1 + K
            for i_m, chk in enumerate(measured):
                if chk in xchecks:
                    continue
                x, y = coords[chk]
                zi = zchecks.index(chk)
                out.append(f"DETECTOR({x}, {y}, 0) rec[-{K - i_m}] rec[-{now - zi}]")
    # Perfect final round, in the order of the last noisy round, compared with it.
    out.append("TICK")
    for chk in measured:
        b = "X" if chk in xchecks else "Z"
        out.append("MPP " + "*".join(f"{b}{q}" for q in support(chk)))
    for i, chk in enumerate(measured):
        x, y = coords[chk]
        out.append(f"DETECTOR({x}, {y}, {rounds + 1}) rec[-{K - i}] rec[-{2 * K - i}]")
    out.append("MPP " + "*".join([f"X{ref}"] + [f"X{q}" for q in x_logical]))
    out.append("MPP " + "*".join([f"Z{ref}"] + [f"Z{q}" for q in z_logical]))
    # The projection, the initial Z_R Z_L, the rounds, the perfect round, the two readouts.
    total = Kz + 1 + rounds * K + K + 2
    out.append("OBSERVABLE_INCLUDE(0) rec[-2]")
    out.append(f"OBSERVABLE_INCLUDE(1) rec[-1] rec[-{total - Kz}]")
    return "\n".join(out) + "\n"
