"""Translation of straight-line ``ftsim`` programs to Qiskit, for cross-checks.

Supported: gates, Pauli channels (as Aer ``pauli_error`` instructions) and measurements. Qiskit
orders qubits little-endian like ``ftsim`` (qubit q is bit q of a basis index); Qiskit Pauli
labels are written right to left, so a two-qubit label "BA" puts A on the first target.
"""

from __future__ import annotations

import re

_GATES = {
    "I": "id", "X": "x", "Y": "y", "Z": "z", "H": "h", "S": "s", "S_DAG": "sdg", "T": "t",
    "T_DAG": "tdg", "SQRT_X": "sx", "SQRT_X_DAG": "sxdg", "CX": "cx", "CY": "cy", "CZ": "cz",
    "SWAP": "swap", "CCX": "ccx", "CCZ": "ccz",
}
_ARITY = {"cx": 2, "cy": 2, "cz": 2, "swap": 2, "ccx": 3, "ccz": 3}
_P = "IXYZ"


def _channel(name: str, args: list[float]):
    """(arity, [(qiskit label, probability)]) of a noise instruction."""
    if name == "PAULI_CHANNEL_1":
        terms = [(_P[j + 1], args[j]) for j in range(3)]
        ar = 1
    elif name == "X_ERROR":
        terms, ar = [("X", args[0])], 1
    elif name == "Y_ERROR":
        terms, ar = [("Y", args[0])], 1
    elif name == "Z_ERROR":
        terms, ar = [("Z", args[0])], 1
    elif name == "DEPOLARIZE1":
        terms, ar = [(c, args[0] / 3) for c in "XYZ"], 1
    elif name in ("PAULI_CHANNEL_2", "DEPOLARIZE2"):
        probs = args if name == "PAULI_CHANNEL_2" else [args[0] / 15] * 15
        terms = []
        for j in range(15):
            k = j + 1
            a, b = _P[k // 4], _P[k % 4]
            terms.append((b + a, probs[j]))
        ar = 2
    elif name == "Z_CHANNEL_3":
        terms = []
        for m in range(1, 8):
            label = "".join("Z" if (m >> q) & 1 else "I" for q in (2, 1, 0))
            terms.append((label, args[m - 1]))
        ar = 3
    else:
        raise ValueError(f"unsupported noise {name}")
    terms = [(lab, p) for lab, p in terms if p > 0]
    total = sum(p for _, p in terms)
    return ar, terms + [("I" * ar, 1 - total)]


def to_qiskit(text: str, *, noise: bool = True, measure: bool = True):
    """(QuantumCircuit, measured qubits in record order)."""
    from qiskit import QuantumCircuit
    from qiskit_aer.noise import pauli_error

    lines = []
    nq = 0
    for raw in text.splitlines():
        line = raw.split("#")[0].strip()
        if not line:
            continue
        m = re.match(r"^([A-Z_0-9]+)(?:\(([^)]*)\))?\s*(.*)$", line)
        name, args, targets = m.group(1), m.group(2), m.group(3)
        qs = [int(t) for t in targets.split()]
        if qs:
            nq = max(nq, max(qs) + 1)
        lines.append((name, [float(a) for a in args.split(",")] if args else [], qs))
    qc = QuantumCircuit(nq)
    measured = []
    for name, args, qs in lines:
        if name in ("TICK",) or name.startswith("MARK"):
            continue
        if name in _GATES:
            g = _GATES[name]
            ar = _ARITY.get(g, 1)
            for i in range(0, len(qs), ar):
                t = qs[i : i + ar]
                if g == "ccz":  # Aer's density-matrix method lacks ccz: H·CCX·H on the target.
                    qc.h(t[2])
                    qc.ccx(*t)
                    qc.h(t[2])
                else:
                    getattr(qc, g)(*t)
        elif name == "M":
            measured.extend(qs)
        elif name in ("REPEAT", "TABLE", "MX", "MY", "R", "RX"):
            raise ValueError(f"{name} is not supported by the Qiskit bridge")
        else:
            ar, terms = _channel(name, args)
            if not noise:
                continue
            err = pauli_error(terms)
            for i in range(0, len(qs), ar):
                qc.append(err.to_instruction(), qs[i : i + ar])
    return qc, measured


def exact_distribution(text: str) -> dict[int, float]:
    """Exact outcome probabilities of the program's measurements, keyed by the record bits as
    an integer (record j = bit j), from Aer's density-matrix simulator."""
    from qiskit_aer import AerSimulator

    qc, measured = to_qiskit(text)
    qc.save_probabilities_dict(qubits=measured, label="p")
    sim = AerSimulator(method="density_matrix")
    res = sim.run(qc, shots=1).result()
    probs = res.data(0)["p"]
    return {int(k, 16) if isinstance(k, str) else int(k): float(v) for k, v in probs.items()}
