"""Sampling and decoding surface-code circuits with stabilizer-qec, reduced to counts of joint
observable-flip patterns (bit j of a pattern = observable j predicted wrongly)."""

from __future__ import annotations

import time
from collections import Counter
from dataclasses import dataclass, field

import numpy as np


@dataclass
class Counts:
    shots: int
    num_obs: int
    patterns: dict[int, int] = field(default_factory=dict)  # non-zero patterns only
    seconds: float = 0.0

    @property
    def failures(self) -> int:
        return sum(self.patterns.values())

    def probabilities(self) -> np.ndarray:
        p = np.zeros(1 << self.num_obs)
        for k, v in self.patterns.items():
            p[k] = v / self.shots
        p[0] = 1 - p[1:].sum()
        return p

    def marginal(self, j: int) -> float:
        return sum(v for k, v in self.patterns.items() if k >> j & 1) / self.shots

    def merge(self, o: "Counts") -> "Counts":
        if o.num_obs != self.num_obs:
            raise ValueError("observable counts differ")
        pats = Counter(self.patterns)
        pats.update(o.patterns)
        return Counts(self.shots + o.shots, self.num_obs, dict(pats), self.seconds + o.seconds)

    def to_json(self) -> dict:
        return {"shots": self.shots, "num_obs": self.num_obs, "patterns": {str(k): v for k, v in sorted(self.patterns.items())}, "seconds": round(self.seconds, 3)}

    @staticmethod
    def from_json(d: dict) -> "Counts":
        return Counts(d["shots"], d["num_obs"], {int(k): v for k, v in d["patterns"].items()}, d.get("seconds", 0.0))


def _packed_ints(rows: np.ndarray) -> np.ndarray:
    """Bit-packed uint8 rows (little-endian bits) → one integer per row."""
    out = np.zeros(rows.shape[0], dtype=np.int64)
    for b in range(rows.shape[1]):
        out |= rows[:, b].astype(np.int64) << (8 * b)
    return out


def run(circuit, *, shots: int, seed: int, target: int | None = None, batch: int = 1 << 16, correlated: bool = True, threads: int = 0) -> Counts:
    """Samples `shots` shots (stopping early once `target` failures are seen, at a multiple of
    `batch`), decodes with (correlated) matching, and counts joint flip patterns."""
    import stabilizer_qec as sq

    c = circuit if isinstance(circuit, sq.Circuit) else sq.Circuit(circuit)
    if c.num_observables > 62:
        raise ValueError("at most 62 observables")
    dem = c.detector_error_model(decompose_errors=True)
    m = sq.Matching(dem, enable_correlations=correlated)
    sampler = c.compile_detector_sampler(seed=seed)
    batch = max(64, batch // 64 * 64)
    pats: Counter = Counter()
    done = 0
    t0 = time.perf_counter()
    while done < shots:
        b = min(batch, shots - done)
        dets, obs = sampler.sample(b, separate_observables=True, bit_packed=True, threads=threads)
        pred = m.decode_batch(dets, bit_packed_shots=True, bit_packed_predictions=True, threads=threads)
        flips = _packed_ints(np.asarray(pred, dtype=np.uint8) ^ np.asarray(obs, dtype=np.uint8))
        nz = flips[flips != 0]
        if nz.size:
            u, cnt = np.unique(nz, return_counts=True)
            for k, v in zip(u.tolist(), cnt.tolist()):
                pats[k] += v
        done += b
        if target is not None and sum(pats.values()) >= target:
            break
    return Counts(done, c.num_observables, dict(pats), time.perf_counter() - t0)


def run_stim(text: str, *, shots: int, seed: int, correlated: bool = True) -> Counts:
    """The same experiment sampled by Stim and decoded by PyMatching (cross-check)."""
    import pymatching
    import stim

    c = stim.Circuit(text)
    dem = c.detector_error_model(decompose_errors=True)
    m = pymatching.Matching.from_detector_error_model(dem, enable_correlations=correlated)
    sampler = c.compile_detector_sampler(seed=seed)
    t0 = time.perf_counter()
    pats: Counter = Counter()
    done = 0
    while done < shots:
        b = min(1 << 16, shots - done)
        dets, obs = sampler.sample(b, separate_observables=True, bit_packed=True)
        # PyMatching needs enable_correlations at decode time as well as at construction.
        pred = m.decode_batch(dets, bit_packed_shots=True, bit_packed_predictions=True, enable_correlations=correlated)
        flips = _packed_ints(np.asarray(pred, dtype=np.uint8) ^ np.asarray(obs, dtype=np.uint8))
        nz = flips[flips != 0]
        if nz.size:
            u, cnt = np.unique(nz, return_counts=True)
            for k, v in zip(u.tolist(), cnt.tolist()):
                pats[k] += v
        done += b
    return Counts(done, c.num_observables, dict(pats), time.perf_counter() - t0)
