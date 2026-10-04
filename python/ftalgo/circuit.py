"""A small IR for logical circuits that keeps REPEAT and TABLE structure.

``Circuit`` is what the algorithm builders produce, what the scheduler (``ftalgo.schedule``)
decorates with noise, and what ``to_text`` turns into the ``ftsim`` program format.

Ops (tuples):

- ``("gate", name, targets)``
- ``("measure", basis, q, rid)`` where ``rid`` is the record id within this circuit's scope
- ``("reset", basis, q)``
- ``("repeat", n, body)`` with ``body`` a ``Circuit`` (its records are local to the body)
- ``("table", rids, cases)`` with ``cases`` a list of 2^len(rids) ``Circuit``s
- ``("mark", label)``
- ``("noise", name, args, targets)`` raw noise instructions (used by the noise compiler)
"""

from __future__ import annotations

from collections import Counter

ARITY = {
    "I": 1, "X": 1, "Y": 1, "Z": 1, "H": 1, "S": 1, "S_DAG": 1, "T": 1, "T_DAG": 1,
    "SQRT_X": 1, "SQRT_X_DAG": 1, "CX": 2, "CY": 2, "CZ": 2, "SWAP": 2, "CCX": 3, "CCZ": 3,
}
DAGGER = {
    "S": "S_DAG", "S_DAG": "S", "T": "T_DAG", "T_DAG": "T", "SQRT_X": "SQRT_X_DAG",
    "SQRT_X_DAG": "SQRT_X",
}


class Circuit:
    def __init__(self, num_qubits: int = 0):
        self.ops: list[tuple] = []
        self.num_qubits = num_qubits
        self.num_records = 0
        self._names: dict[str, list[int]] = {}

    # -- registers -----------------------------------------------------------------------------
    def alloc(self, n: int, name: str) -> list[int]:
        """n fresh qubits (in |0⟩), named for reports."""
        qs = list(range(self.num_qubits, self.num_qubits + n))
        self.num_qubits += n
        self._names[name] = qs
        return qs

    @property
    def registers(self) -> dict[str, list[int]]:
        return dict(self._names)

    def _touch(self, *qs: int) -> None:
        for q in qs:
            if q < 0:
                raise ValueError(f"negative qubit {q}")
            if q >= self.num_qubits:
                self.num_qubits = q + 1

    # -- instructions --------------------------------------------------------------------------
    def gate(self, name: str, *targets: int) -> None:
        if name not in ARITY:
            raise ValueError(f"unknown gate {name}")
        if len(targets) != ARITY[name]:
            raise ValueError(f"{name} takes {ARITY[name]} targets, got {len(targets)}")
        if len(set(targets)) != len(targets):
            raise ValueError(f"{name} on repeated qubits {targets}")
        self._touch(*targets)
        self.ops.append(("gate", name, tuple(targets)))

    def x(self, q): self.gate("X", q)
    def y(self, q): self.gate("Y", q)
    def z(self, q): self.gate("Z", q)
    def h(self, q): self.gate("H", q)
    def s(self, q): self.gate("S", q)
    def sdg(self, q): self.gate("S_DAG", q)
    def t(self, q): self.gate("T", q)
    def tdg(self, q): self.gate("T_DAG", q)
    def cx(self, c, t): self.gate("CX", c, t)
    def cz(self, a, b): self.gate("CZ", a, b)
    def swap(self, a, b): self.gate("SWAP", a, b)
    def ccx(self, a, b, t): self.gate("CCX", a, b, t)
    def ccz(self, a, b, c): self.gate("CCZ", a, b, c)

    def cswap(self, c: int, a: int, b: int) -> None:
        """Fredkin: CX(b,a) CCX(c,a,b) CX(b,a)."""
        self.cx(b, a)
        self.ccx(c, a, b)
        self.cx(b, a)

    def gates(self, seq: list[str], q: int) -> None:
        """A sequence of single-qubit gates on q, in time order."""
        for g in seq:
            self.gate(g, q)

    def measure(self, q: int, basis: str = "Z") -> int:
        self._touch(q)
        rid = self.num_records
        self.ops.append(("measure", basis, q, rid))
        self.num_records += 1
        return rid

    def reset(self, q: int, basis: str = "Z") -> None:
        self._touch(q)
        self.ops.append(("reset", basis, q))

    def mark(self, label: str) -> None:
        self.ops.append(("mark", label))

    def noise(self, name: str, args: tuple, *targets: int) -> None:
        self._touch(*targets)
        self.ops.append(("noise", name, tuple(args), tuple(targets)))

    def repeat(self, n: int, body: "Circuit") -> None:
        if n < 1:
            raise ValueError("repeat count must be ≥ 1")
        if not body.ops:
            return
        self._touch(*range(body.num_qubits))
        self.ops.append(("repeat", n, body))
        self.num_records += n * body.num_records

    def table(self, rids: list[int], cases: list["Circuit"]) -> None:
        if len(cases) != 1 << len(rids):
            raise ValueError(f"TABLE on {len(rids)} records needs {1 << len(rids)} cases")
        for r in rids:
            if not 0 <= r < self.num_records:
                raise ValueError(f"record {r} does not exist yet")
        m = {c.num_records for c in cases}
        if len(m) != 1:
            raise ValueError("TABLE cases must have equal measurement counts")
        for c in cases:
            self._touch(*range(c.num_qubits))
        self.ops.append(("table", tuple(rids), cases))
        self.num_records += m.pop()

    def extend(self, other: "Circuit") -> None:
        """Appends other's ops; other's record ids are shifted into this scope."""
        shift = self.num_records
        for op in other.ops:
            if op[0] == "measure":
                self.ops.append(("measure", op[1], op[2], op[3] + shift))
            elif op[0] == "table":
                self.ops.append(("table", tuple(r + shift for r in op[1]), op[2]))
            else:
                self.ops.append(op)
        self._touch(*range(other.num_qubits))
        self.num_records += other.num_records

    def inverse(self) -> "Circuit":
        """The inverse of a gate-only circuit (MARKs kept in place)."""
        inv = Circuit(self.num_qubits)
        for op in reversed(self.ops):
            kind = op[0]
            if kind == "gate":
                inv.ops.append(("gate", DAGGER.get(op[1], op[1]), op[2]))
            elif kind == "repeat":
                inv.ops.append(("repeat", op[1], op[2].inverse()))
            elif kind == "mark":
                inv.ops.append(op)
            else:
                raise ValueError(f"cannot invert a circuit containing {kind}")
        return inv

    # -- output --------------------------------------------------------------------------------
    def to_text(self) -> str:
        lines: list[str] = []
        self._emit(lines, "")
        return "\n".join(lines) + "\n"

    def _emit(self, lines: list[str], ind: str) -> None:
        count = 0
        pending: list[str] | None = None  # merge consecutive identical gate names on one line

        def flush():
            nonlocal pending
            if pending:
                lines.append(ind + " ".join(pending))
            pending = None

        for op in self.ops:
            kind = op[0]
            if kind == "gate":
                name, ts = op[1], op[2]
                if pending is not None and pending[0] == name and len(pending) < 400:
                    pending.extend(map(str, ts))
                else:
                    flush()
                    pending = [name, *map(str, ts)]
                continue
            flush()
            if kind == "measure":
                lines.append(f"{ind}M{'' if op[1] == 'Z' else op[1]} {op[2]}")
                count += 1
            elif kind == "reset":
                lines.append(f"{ind}R{'' if op[1] == 'Z' else op[1]} {op[2]}")
            elif kind == "mark":
                lines.append(f"{ind}MARK({op[1]})")
            elif kind == "noise":
                args = ",".join(repr(float(a)) for a in op[2])
                lines.append(f"{ind}{op[1]}({args}) " + " ".join(map(str, op[3])))
            elif kind == "repeat":
                lines.append(f"{ind}REPEAT {op[1]} {{")
                op[2]._emit(lines, ind + "  ")
                lines.append(f"{ind}}}")
                count += op[1] * op[2].num_records
            elif kind == "table":
                refs = " ".join(f"rec[-{count - r}]" for r in op[1])
                lines.append(f"{ind}TABLE {refs} {{")
                for i, case in enumerate(op[2]):
                    if i:
                        lines.append(f"{ind}}} {{")
                    case._emit(lines, ind + "  ")
                lines.append(f"{ind}}}")
                count += op[2][0].num_records
        flush()

    def counts(self) -> Counter:
        """Gate counts (plus "M", "R"); REPEAT multiplies, a TABLE counts its largest case per kind."""
        c: Counter = Counter()
        for op in self.ops:
            kind = op[0]
            if kind == "gate":
                c[op[1]] += 1
            elif kind == "measure":
                c["M"] += 1
            elif kind == "reset":
                c["R"] += 1
            elif kind == "repeat":
                for k, v in op[2].counts().items():
                    c[k] += op[1] * v
            elif kind == "table":
                best: Counter = Counter()
                for case in op[2]:
                    for k, v in case.counts().items():
                        best[k] = max(best[k], v)
                c.update(best)
        return c

    def toffoli_count(self) -> int:
        cnt = self.counts()
        return cnt["CCX"] + cnt["CCZ"]

    def t_count(self) -> int:
        cnt = self.counts()
        return cnt["T"] + cnt["T_DAG"]
