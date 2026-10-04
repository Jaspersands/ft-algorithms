"""numpy front end of the Rust engine (``ftalgo._ftsim``)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Union

import numpy as np

from . import _ftsim


@dataclass
class Faults:
    """Faults that fired: one entry per fault. ``site`` is the static site index, or
    2**64 − 1 for a site inside a TABLE case; ``pauli`` packs two bits per target
    (0 = I, 1 = X, 2 = Y, 3 = Z), first target lowest; ``cls`` indexes ``Program.classes``."""

    shot: np.ndarray
    site: np.ndarray
    cls: np.ndarray
    pauli: np.ndarray

    DYNAMIC = np.iinfo(np.uint64).max


class Program:
    """A logical program in the ``ftsim`` text format."""

    def __init__(self, text: str):
        self.text = text
        self._p = _ftsim.Program(text)

    @property
    def num_qubits(self) -> int:
        return self._p.num_qubits

    @property
    def num_measurements(self) -> int:
        return self._p.num_measurements

    @property
    def num_sites(self) -> int:
        """Static noise sites (target groups of noise instructions outside TABLE cases)."""
        return self._p.num_sites

    @property
    def expected_faults(self) -> float:
        """Σ q over static sites."""
        return self._p.expected_faults

    @property
    def classes(self) -> list[str]:
        return list(self._p.classes)

    def fault_count_distribution(self, kmax: int) -> np.ndarray:
        """[P(K=0), …, P(K=kmax), P(K>kmax)] for the number K of static faults."""
        return np.array(self._p.fault_count_distribution(kmax))

    def sample(
        self,
        shots: int,
        *,
        seed: int,
        faults: Union[str, int] = "plain",
        backend: str = "auto",
        threads: int = 0,
        first_shot: int = 0,
        log_faults: bool = False,
    ):
        """Runs ``shots`` shots. ``faults``: "plain" (every site independently), "none", or an
        integer k (exactly k static faults, from their conditional law; TABLE-case noise still
        fires independently). Returns a (shots, num_measurements) bool array, and with
        ``log_faults`` also a :class:`Faults`."""
        d = self._p.sample(
            shots,
            seed=seed,
            faults=faults,
            backend=backend,
            threads=threads,
            first_shot=first_shot,
            log_faults=log_faults,
        )
        m = self.num_measurements
        bps = d["bytes_per_shot"]
        raw = np.frombuffer(d["records"], dtype=np.uint8).reshape(shots, bps)
        rec = np.unpackbits(raw, axis=1, bitorder="little")[:, :m].astype(bool)
        if not log_faults:
            return rec
        f = Faults(
            shot=np.frombuffer(d["fault_shot"], dtype=np.uint64),
            site=np.frombuffer(d["fault_site"], dtype=np.uint64),
            cls=np.frombuffer(d["fault_class"], dtype=np.uint16),
            pauli=np.frombuffer(d["fault_pauli"], dtype=np.uint32),
        )
        return rec, f

    def final_state(self, *, seed: int = 0, backend: str = "auto"):
        """One noiseless shot: ({basis index: amplitude}, records)."""
        idx, re, im, recs = self._p.final_state(seed=seed, backend=backend)
        return {int(i): complex(a, b) for i, a, b in zip(idx, re, im)}, np.array(recs, dtype=bool)

    def statevector(self, *, seed: int = 0, backend: str = "auto") -> np.ndarray:
        """The final state of one noiseless shot as a dense vector (qubit q = bit q)."""
        amps, _ = self.final_state(seed=seed, backend=backend)
        if self.num_qubits > 26:
            raise ValueError("statevector: too many qubits for a dense vector; use final_state")
        v = np.zeros(1 << max(self.num_qubits, 1), dtype=complex)
        for i, a in amps.items():
            v[i] = a
        return v

    def __repr__(self) -> str:
        return f"Program(num_qubits={self.num_qubits}, num_measurements={self.num_measurements}, num_sites={self.num_sites})"
