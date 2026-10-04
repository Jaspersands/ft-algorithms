"""Single-qubit Z rotations as Clifford + T sequences (Ross–Selinger gridsynth).

``rz(theta, eps)`` returns gates, in time order, whose product is within ``eps`` (operator norm,
minimized over global phase) of diag(1, e^{iθ}); every sequence is checked as a matrix before
it is returned. Multiples of π/4 are exact and need no search. Results are cached on disk
(``data/synth/cache.json``), keyed by the exact float repr of θ and ε.
"""

from __future__ import annotations

import json
import math
import os
import threading

import numpy as np

_H = np.array([[1, 1], [1, -1]], dtype=complex) / math.sqrt(2)
MATS = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.diag([1, -1]).astype(complex),
    "H": _H,
    "S": np.diag([1, 1j]),
    "S_DAG": np.diag([1, -1j]),
    "T": np.diag([1, np.exp(1j * math.pi / 4)]),
    "T_DAG": np.diag([1, np.exp(-1j * math.pi / 4)]),
}
_EXACT = {0: [], 1: ["T"], 2: ["S"], 3: ["S", "T"], 4: ["Z"], 5: ["Z", "T"], 6: ["S_DAG"], 7: ["T_DAG"]}

_CACHE_PATH = os.environ.get(
    "FTALGO_SYNTH_CACHE",
    os.path.join(os.path.dirname(__file__), "..", "..", "data", "synth", "cache.json"),
)
_cache: dict[str, str] | None = None
_dirty = False
_lock = threading.Lock()


def matrix(seq: list[str]) -> np.ndarray:
    """The unitary of a gate sequence applied in time order."""
    u = np.eye(2, dtype=complex)
    for g in seq:
        u = MATS[g] @ u
    return u


def distance_up_to_phase(u: np.ndarray, v: np.ndarray) -> float:
    """min over φ of ‖u − e^{iφ} v‖ (operator norm), for 2×2 unitaries."""
    tr = np.trace(v.conj().T @ u)
    ph = tr / abs(tr) if abs(tr) > 1e-300 else 1.0
    return float(np.linalg.norm(u - ph * v, 2))


def target(theta: float) -> np.ndarray:
    return np.diag([1, np.exp(1j * theta)])


def _load() -> dict[str, str]:
    global _cache
    if _cache is None:
        try:
            with open(_CACHE_PATH) as f:
                _cache = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError):
            _cache = {}
    return _cache


def save_cache() -> None:
    global _dirty
    with _lock:
        if not _dirty or _cache is None:
            return
        os.makedirs(os.path.dirname(_CACHE_PATH), exist_ok=True)
        tmp = _CACHE_PATH + ".tmp"
        with open(tmp, "w") as f:
            json.dump(_cache, f, sort_keys=True, separators=(",", ":"))
        os.replace(tmp, _CACHE_PATH)
        _dirty = False


def _gridsynth(theta: float, eps: float) -> list[str]:
    import mpmath
    from pygridsynth.gridsynth import gridsynth_gates

    mpmath.mp.dps = max(30, int(-math.log10(eps) * 3) + 20)
    # gridsynth approximates Rz(θ) = e^{−iθZ/2} = e^{−iθ/2}·diag(1, e^{iθ}); "W" is a global phase.
    s = gridsynth_gates(theta=mpmath.mpf(repr(theta)), epsilon=mpmath.mpf(repr(eps)))
    # The string is a matrix product (leftmost applied last), so time order is reversed.
    return [c for c in reversed(s) if c != "W"]


def rz(theta: float, eps: float) -> list[str]:
    """Clifford + T gates (time order) ≈ diag(1, e^{iθ}) within eps up to global phase."""
    global _dirty
    theta = float(theta)
    eps = float(eps)
    if not (0 < eps < 0.5):
        raise ValueError("eps must be in (0, 0.5)")
    k = theta / (math.pi / 4)
    kr = round(k)
    if abs(k - kr) < 1e-12:
        return list(_EXACT[kr % 8])
    key = f"{theta!r}|{eps!r}"
    cache = _load()
    with _lock:
        hit = cache.get(key)
    if hit is not None:
        seq = hit.split(",") if hit else []
    else:
        seq = _gridsynth(theta, eps)
        with _lock:
            cache[key] = ",".join(seq)
            _dirty = True
    d = distance_up_to_phase(matrix(seq), target(theta))
    if d > eps * (1 + 1e-6) + 1e-14:
        raise RuntimeError(f"synthesized sequence for θ={theta} misses by {d} > ε={eps}")
    return seq


def t_count(seq: list[str]) -> int:
    return sum(1 for g in seq if g in ("T", "T_DAG"))
