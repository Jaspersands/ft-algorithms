"""A bit-level simulator for reversible classical circuits (X, CX, CCX, SWAP), independent of the
state-vector engine: arithmetic blocks are verified on every input with it.

Phase-only gates (Z, S, T, CZ, CCZ and daggers) do not change a basis state's label; they are
accepted only with ``ignore_phases=True``.
"""

from __future__ import annotations

from .circuit import Circuit

_PHASES = {"Z", "S", "S_DAG", "T", "T_DAG", "CZ", "CCZ", "I"}


def run_reversible(c: Circuit, state: int, *, ignore_phases: bool = False) -> int:
    for op in c.ops:
        kind = op[0]
        if kind == "gate":
            name, t = op[1], op[2]
            if name == "X":
                state ^= 1 << t[0]
            elif name == "CX":
                if state >> t[0] & 1:
                    state ^= 1 << t[1]
            elif name == "CCX":
                if state >> t[0] & 1 and state >> t[1] & 1:
                    state ^= 1 << t[2]
            elif name == "SWAP":
                a, b = state >> t[0] & 1, state >> t[1] & 1
                if a != b:
                    state ^= (1 << t[0]) | (1 << t[1])
            elif name in _PHASES and (ignore_phases or name == "I"):
                pass
            else:
                raise ValueError(f"run_reversible: {name} is not a classical reversible gate")
        elif kind == "repeat":
            for _ in range(op[1]):
                state = run_reversible(op[2], state, ignore_phases=ignore_phases)
        elif kind == "mark":
            pass
        else:
            raise ValueError(f"run_reversible: {kind} not supported")
    return state


def get(state: int, reg: list[int]) -> int:
    """The integer held by a little-endian register."""
    v = 0
    for i, q in enumerate(reg):
        v |= (state >> q & 1) << i
    return v


def put(state: int, reg: list[int], value: int) -> int:
    for i, q in enumerate(reg):
        bit = value >> i & 1
        state = (state & ~(1 << q)) | (bit << q)
    return state
