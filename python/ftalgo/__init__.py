"""Quantum algorithms on a simulated fault-tolerant quantum computer.

``ftalgo.Program`` runs logical Clifford + T + Toffoli programs (the ``ftsim`` text format) on a
sparse or dense state vector with Pauli-channel noise, sampled per site or stratified by the
exact number of faults.
"""

from ._ftsim import __version__ as engine_version
from .engine import Faults, Program

__all__ = ["Program", "Faults", "engine_version"]
