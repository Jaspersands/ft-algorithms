"""Lithium hydride (LiH) ab initio Hamiltonian, active space reduction, and dissociation curve.

LiH (4 electrons, 6 spatial STO-3G orbitals):
- 1s² core on Li is frozen, shifting the 1-body integrals and nuclear repulsion.
- Active space: 2 electrons in 2 spatial orbitals (HOMO and LUMO: σ and σ*).
- Under Jordan-Wigner, this yields 4 qubits.
- Under Z₂ spin-parity tapering (spin-up and spin-down parity), this tapers to 2 qubits.
- Orbitals: CASSCF(2,2)-optimised by default (within ~2 mHa of STO-3G FCI along the curve);
  the active-space Hamiltonian reproduces PySCF CASSCF/CASCI energies to < 1e-10 Hartree.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

import numpy as np
from pyscf import ao2mo, fci, gto, mcscf, scf

from .qubit import Molecule, PauliSum, QubitProblem, fermion_hamiltonian, taper


@dataclass
class LiHDissociationPoint:
    r_angstrom: float
    r_bohr: float
    e_hf: float
    e_casci: float
    e_casscf: float
    e_fci: float
    e_nuc: float


def _casscf(mf, ncas: int, mo_guess=None):
    """Orbital-optimised CAS(2, ncas) on top of RHF; returns the converged CASSCF object."""
    mc = mcscf.CASSCF(mf, ncas, 2)
    mc.verbose = 0
    if mo_guess is not None:
        mo_guess = mcscf.project_init_guess(mc, mo_guess)
    mc.kernel(mo_guess)
    if not mc.converged:
        raise RuntimeError("CASSCF did not converge")
    return mc


def lih_problem(r_angstrom: float = 1.595, ncas: int = 2, orbitals: str = "casscf", mo_guess=None) -> QubitProblem:
    """Constructs the active-space tapered Hamiltonian for LiH at internuclear distance r_angstrom.

    Freezes the 1s² core orbital of Li (orbital 0), keeping ncas active valence orbitals.
    For ncas = 2 (default), active space has 2 electrons in 2 spatial orbitals (4 spin-orbitals),
    which tapers to 2 qubits with 9 Pauli terms.

    orbitals = "casscf" (default) uses orbital-optimised CASSCF(2, ncas) orbitals: with ncas = 2 this
    stays within ~2 mHa of full STO-3G FCI over R in [1, 3] Å. orbitals = "rhf" uses canonical RHF
    orbitals, whose CAS(2, 2) recovers well under 1 mHa of correlation (kept for comparison).
    The exact ground energy of the returned Hamiltonian is stored as e_fci (the active-space target).
    """
    r_bohr = r_angstrom * 1.8897261254578281
    mol = gto.M(atom=f"Li 0 0 0; H 0 0 {r_angstrom}", basis="sto-3g", verbose=0)
    mf = scf.RHF(mol).run()

    if orbitals == "casscf":
        mc = _casscf(mf, ncas, mo_guess)
        mo = mc.mo_coeff
        e_active = float(mc.e_tot)
    elif orbitals == "rhf":
        mo = mf.mo_coeff
        e_active = float(mcscf.CASCI(mf, ncas=ncas, nelecas=2).kernel()[0])
    else:
        raise ValueError(f"unknown orbitals {orbitals!r}")

    # Transform integrals to MO basis
    h_mo = mo.T @ mf.get_hcore() @ mo
    eri_mo = ao2mo.kernel(mol, mo, compact=False).reshape((mol.nao, mol.nao, mol.nao, mol.nao))

    core = [0]
    active = list(range(1, 1 + ncas))
    e_nuc = float(mol.energy_nuc())

    # Core energy shift: E_core = E_nuc + 2 Σ_i h_ii + Σ_ij (2 (ii|jj) - (ij|ji))
    e_core = e_nuc + sum(2.0 * h_mo[i, i] for i in core)
    for i in core:
        for j in core:
            e_core += 2.0 * eri_mo[i, i, j, j] - eri_mo[i, j, j, i]

    # Effective 1-body Hamiltonian: h_eff_pq = h_pq + Σ_i (2 (pq|ii) - (pi|iq))
    h_eff = np.zeros((ncas, ncas))
    for p_idx, p in enumerate(active):
        for q_idx, q in enumerate(active):
            val = float(h_mo[p, q])
            for i in core:
                val += 2.0 * eri_mo[p, q, i, i] - eri_mo[p, i, i, q]
            h_eff[p_idx, q_idx] = val

    # 2-body Hamiltonian in active space: (pq|rs)
    g_eff = eri_mo[np.ix_(active, active, active, active)]

    full = fermion_hamiltonian(h_eff, g_eff, e_core)
    hf_bits = (1 << 2) - 1  # 2 electrons in lowest active mode (modes 0 and 1)
    tapered, hf_t, syms = taper(full, hf_bits)

    m_obj = Molecule("LiH", [("Li", (0.0, 0.0, 0.0)), ("H", (0.0, 0.0, r_bohr))], 0)
    return QubitProblem(
        molecule=m_obj,
        e_hf=float(mf.e_tot),
        e_fci=e_active,
        hamiltonian=tapered,
        hf_bits=hf_t,
        full=full,
        full_hf_bits=hf_bits,
        symmetries=syms,
    )


def lih_dissociation_curve(r_values: Sequence[float] | None = None) -> list[LiHDissociationPoint]:
    """Computes the dissociation potential energy curve of LiH across a range of bond lengths.

    Each point carries RHF, CASCI(2,2) on RHF orbitals, CASSCF(2,2) and full STO-3G FCI energies.
    CASSCF orbitals are carried along the curve as the initial guess for the next geometry.
    """
    if r_values is None:
        r_values = [1.0, 1.2, 1.4, 1.595, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0]
    pts = []
    mo = None
    for r in r_values:
        mol = gto.M(atom=f"Li 0 0 0; H 0 0 {r}", basis="sto-3g", verbose=0)
        mf = scf.RHF(mol).run()
        e_cas = float(mcscf.CASCI(mf, ncas=2, nelecas=2).kernel()[0])
        mc = _casscf(mf, 2, mo)
        mo = mc.mo_coeff
        e_full_fci = float(fci.FCI(mf).kernel()[0])
        pts.append(
            LiHDissociationPoint(
                r_angstrom=float(r),
                r_bohr=float(r * 1.8897261254578281),
                e_hf=float(mf.e_tot),
                e_casci=e_cas,
                e_casscf=float(mc.e_tot),
                e_fci=e_full_fci,
                e_nuc=float(mol.energy_nuc()),
            )
        )
    return pts
