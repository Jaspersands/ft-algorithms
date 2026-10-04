"""The noise compiler: an ideal logical ``Circuit`` → a scheduled noisy ``Circuit`` for ftsim.

Each qubit has a clock (in rounds). An operation starts when all its qubits are free (and, for a
TABLE, when the records it reads are available plus the reaction time); every qubit that waited
receives the composed idle channel for the gap; the operation's own channel follows it, marked
with its class (``MARK``) for fault attribution. A qubit starts its life at its first operation.

REPEAT is a barrier: all live qubits are synchronized before it, the body is scheduled once with
every live qubit idling to the body's end, and the clocks advance by n times its duration. A
TABLE pads every case with idling to the latest case's end, so the schedule does not depend on
measurement outcomes.

The expected number of faults per class (Σ of site probabilities, REPEATs multiplied) is
accumulated as the compiler emits noise: the run's error budget before any sampling.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field

from .arch import Architecture, compose_p1
from .circuit import Circuit

FREE = {"I", "X", "Y", "Z", "H"}
MERGE_ONE = {"S", "S_DAG", "T", "T_DAG"}
TWO = {"CX", "CZ", "CY"}
THREE = {"CCX", "CCZ"}


@dataclass
class NoisyProgram:
    circuit: Circuit
    rounds: int
    budget: Counter
    counts: Counter
    num_qubits: int
    arch: Architecture
    t_count: int = 0
    ccz_count: int = 0
    idle_qubit_rounds: float = 0.0

    @property
    def text(self) -> str:
        return self.circuit.to_text()

    @property
    def seconds(self) -> float:
        return self.rounds * 1e-6

    def physical_qubits(self) -> dict:
        data = self.arch.data_qubits(self.num_qubits)
        fac = self.arch.factory_qubits(self.t_count, self.ccz_count, self.rounds)
        return {"data": data, "factories": fac, "total": data + fac}


def _p1(q) -> tuple:
    return tuple(float(max(v, 0.0)) for v in q)


@dataclass
class _Scope:
    clock: dict[int, int]
    rec_time: list[int] = field(default_factory=list)


class _Compiler:
    def __init__(self, arch: Architecture):
        self.arch = arch
        self.budget: Counter = Counter()
        self.counts: Counter = Counter()
        self.t_count = 0
        self.ccz_count = 0
        self.idle_qr = 0.0

    # -- emission helpers ------------------------------------------------------------------------
    def noise1(self, out: Circuit, cls: str, q: int, chan, mult: int) -> None:
        chan = _p1(chan)
        tot = sum(chan)
        if tot <= 0:
            return
        out.mark(cls)
        out.noise("PAULI_CHANNEL_1", chan, q)
        self.budget[cls] += tot * mult

    def idle_to(self, out: Circuit, sc: _Scope, q: int, t: int, mult: int) -> None:
        cur = sc.clock.get(q)
        if cur is None:
            sc.clock[q] = t
            return
        gap = t - cur
        if gap > 0:
            self.noise1(out, "idle", q, self.arch.idle(gap), mult)
            self.idle_qr += gap * mult
            sc.clock[q] = t

    def start(self, out: Circuit, sc: _Scope, qs, extra: int = 0, mult: int = 1) -> int:
        s = max([sc.clock.get(q, 0) or 0 for q in qs] + [extra])
        for q in qs:
            self.idle_to(out, sc, q, s, mult)
        return s

    # -- the walk --------------------------------------------------------------------------------
    def block(self, circ: Circuit, out: Circuit, sc: _Scope, mult: int) -> None:
        a = self.arch
        d = a.d
        for op in circ.ops:
            kind = op[0]
            if kind == "gate":
                name, ts = op[1], op[2]
                self.counts[name] += mult
                if name in FREE:
                    out.gate(name, *ts)
                elif name in MERGE_ONE:
                    s = self.start(out, sc, ts, mult=mult)
                    out.gate(name, *ts)
                    if name.startswith("T"):
                        self.noise1(out, "T", ts[0], a.t_channel(), mult)
                        self.t_count += mult
                    else:
                        self.noise1(out, "S", ts[0], a.s_channel(), mult)
                    sc.clock[ts[0]] = s + d
                elif name in TWO or name == "SWAP":
                    s = self.start(out, sc, ts, mult=mult)
                    reps = 3 if name == "SWAP" else 1
                    out.gate(name, *ts)
                    out.mark("CNOT")
                    for r in range(reps):
                        pair = ts if r % 2 == 0 else ts[::-1]
                        out.noise("PAULI_CHANNEL_2", tuple(a.cnot), *pair)
                        self.budget["CNOT"] += sum(a.cnot) * mult
                    for q in ts:
                        sc.clock[q] = s + 2 * d * reps
                elif name in THREE:
                    s = self.start(out, sc, ts, mult=mult)
                    out.gate(name, *ts)
                    z3, qx = a.ccz_channel()
                    out.mark("CCZ")
                    if name == "CCX":  # CCX = H_t CCZ H_t: the target's Z-type errors are X-type
                        out.h(ts[2])
                        out.noise("Z_CHANNEL_3", tuple(z3), *ts)
                        out.h(ts[2])
                    else:
                        out.noise("Z_CHANNEL_3", tuple(z3), *ts)
                    self.budget["CCZ"] += sum(z3) * mult
                    for q in ts:
                        self.noise1(out, "CCZ", q, (qx, 0.0, 0.0), mult)
                        self.noise1(out, "idle", q, a.idle(d), mult)
                        sc.clock[q] = s + 2 * d
                    self.ccz_count += mult
                else:
                    raise ValueError(f"no operation model for {name}")
            elif kind == "measure":
                basis, q = op[1], op[2]
                s = self.start(out, sc, [q], mult=mult)
                self.noise1(out, "measure", q, a.idle(1), mult)
                out.measure(q, basis)
                sc.clock[q] = s + 1
                sc.rec_time.append(s + 1)
                self.counts["M"] += mult
            elif kind == "reset":
                basis, q = op[1], op[2]
                s = self.start(out, sc, [q], mult=mult)
                out.reset(q, basis)
                self.noise1(out, "measure", q, a.idle(1), mult)
                sc.clock[q] = s + 1
                self.counts["R"] += mult
            elif kind == "mark":
                pass  # algorithm-region marks are not noise classes
            elif kind == "noise":
                out.ops.append(op)
            elif kind == "table":
                self.table(op, out, sc, mult)
            elif kind == "repeat":
                self.repeat(op, out, sc, mult)
            else:
                raise ValueError(kind)

    def table(self, op, out: Circuit, sc: _Scope, mult: int) -> None:
        rids, cases = op[1], op[2]
        qs = sorted({q for c in cases for q in _qubits(c)})
        ready = max(sc.rec_time[r] for r in rids) + self.arch.reaction
        s = self.start(out, sc, qs, extra=ready, mult=mult)
        outs, ends = [], []
        # Dynamic noise inside the cases does not enter the static budget; count the largest case.
        saved = (self.budget.copy(), self.counts.copy(), self.t_count, self.ccz_count)
        best = None
        for c in cases:
            self.budget, self.counts, self.t_count, self.ccz_count = saved[0].copy(), saved[1].copy(), saved[2], saved[3]
            o = Circuit(out.num_qubits)
            csc = _Scope({q: s for q in qs}, [])
            self.block(c, o, csc, mult)
            outs.append((o, csc))
            ends.append(max(csc.clock.values(), default=s))
            if best is None or sum(self.budget.values()) > sum(best[0].values()):
                best = (self.budget.copy(), self.counts.copy(), self.t_count, self.ccz_count)
        self.budget, self.counts, self.t_count, self.ccz_count = best
        end = max(ends)
        padded = []
        for o, csc in outs:
            for q in qs:
                gap = end - csc.clock.get(q, s)
                if gap > 0:
                    o.mark("idle")
                    o.noise("PAULI_CHANNEL_1", _p1(self.arch.idle(gap)), q)
                    if o is outs[0][0]:
                        self.idle_qr += gap * mult
            padded.append(o)
        out.table(list(rids), padded)
        for q in qs:
            sc.clock[q] = end
        m = cases[0].num_records
        sc.rec_time.extend([end] * m)

    def repeat(self, op, out: Circuit, sc: _Scope, mult: int) -> None:
        n, body = op[1], op[2]
        live = sorted(set(sc.clock) | set(_qubits(body)))
        s = self.start(out, sc, live, mult=mult)
        bsc = _Scope({q: 0 for q in live}, [])
        bo = Circuit(out.num_qubits)
        self.block(body, bo, bsc, mult * n)
        D = max(bsc.clock.values(), default=0)
        for q in live:
            self.idle_to(bo, bsc, q, D, mult * n)
        out.repeat(n, bo)
        for q in live:
            sc.clock[q] = s + n * D
        sc.rec_time.extend([s + n * D] * (n * body.num_records))


def _qubits(c: Circuit) -> set[int]:
    qs: set[int] = set()
    for op in c.ops:
        if op[0] == "gate":
            qs.update(op[2])
        elif op[0] in ("measure", "reset"):
            qs.add(op[2])
        elif op[0] == "noise":
            qs.update(op[3])
        elif op[0] == "repeat":
            qs |= _qubits(op[2])
        elif op[0] == "table":
            for cc in op[2]:
                qs |= _qubits(cc)
    return qs


def compile_noisy(circ: Circuit, arch: Architecture) -> NoisyProgram:
    comp = _Compiler(arch)
    out = Circuit(circ.num_qubits)
    sc = _Scope({}, [])
    comp.block(circ, out, sc, 1)
    rounds = max(sc.clock.values(), default=0)
    return NoisyProgram(out, rounds, comp.budget, comp.counts, circ.num_qubits, arch, comp.t_count, comp.ccz_count, comp.idle_qr)


__all__ = ["compile_noisy", "NoisyProgram", "compose_p1"]
