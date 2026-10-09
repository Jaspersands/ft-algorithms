"""Lithium Hydride (LiH) active-space, frozen-core, and Z2 tapering module.

Re-exports ftalgo.chem.lih functions for top-level convenience.
"""

from ftalgo.chem.lih import (
    LiHDissociationPoint,
    lih_dissociation_curve,
    lih_problem,
)

__all__ = [
    "LiHDissociationPoint",
    "lih_problem",
    "lih_dissociation_curve",
]
