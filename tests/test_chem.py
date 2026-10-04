"""Chemistry checked against Szabo & Ostlund's published values and against PySCF."""

import numpy as np
import pytest

from ftalgo.chem.basis import ao_integrals
from ftalgo.chem.qubit import (
    PauliSum,
    jw_annihilation,
    jw_creation,
    molecule,
    pauli_mul,
    qubit_problem,
)
from ftalgo.chem.scf import rhf

MOLS = ["H2", "HeH+", "H3+", "H4"]


def test_szabo_ostlund_h2():
    # Szabo & Ostlund §3.5.2 / Table 4.x: H2, STO-3G, R = 1.4 bohr.
    prob = qubit_problem(molecule("H2", szabo=True))
    assert prob.e_hf == pytest.approx(-1.1167, abs=1e-4)
    assert prob.e_fci == pytest.approx(-1.1373, abs=1e-4)


def test_szabo_ostlund_heh_plus():
    # Szabo & Ostlund §3.5.3: HeH+ at R = 1.4632 bohr with ζ_He = 2.0925, ζ_H = 1.24:
    # electronic energy −4.227529, total −2.860662 hartree. We get −2.8606586; our integrals match
    # PySCF to 1e-8 (test_against_pyscf), so the 3e-6 residual is the textbook program's own
    # numerics, not ours; the precise check is the PySCF one.
    mol = molecule("HeH+", szabo=True)
    ints = ao_integrals(mol.atoms, mol.basis())
    scf = rhf(ints, 2)
    assert scf.energy == pytest.approx(-2.860662, abs=5e-6)


def test_pauli_algebra():
    n = 2
    for a in range(16):
        for b in range(16):
            ka = (a & 3, a >> 2)
            kb = (b & 3, b >> 2)
            ph, k = pauli_mul(ka, kb, n)
            A = PauliSum(n, {ka: 1}).matrix()
            B = PauliSum(n, {kb: 1}).matrix()
            C = PauliSum(n, {k: ph}).matrix()
            assert np.allclose(A @ B, C)


def test_jw_anticommutation():
    n = 4
    for p in range(n):
        for q in range(n):
            ap, aq = jw_annihilation(p, n).matrix(), jw_annihilation(q, n).matrix()
            adq = jw_creation(q, n).matrix()
            assert np.allclose(ap @ adq + adq @ ap, np.eye(16) * (p == q))
            assert np.allclose(ap @ aq + aq @ ap, 0)


pyscf = pytest.importorskip("pyscf")


def pyscf_mol(mol):
    from pyscf import gto

    m = gto.M(
        atom=[(e, c) for e, c in mol.atoms],
        basis="sto-3g",
        unit="Bohr",
        charge=mol.charge,
        spin=0,
        verbose=0,
    )
    return m


@pytest.mark.parametrize("name", MOLS)
def test_against_pyscf(name):
    from pyscf import fci, scf

    mol = molecule(name)
    pm = pyscf_mol(mol)
    ints = ao_integrals(mol.atoms, mol.basis())
    assert np.allclose(ints.S, pm.intor("int1e_ovlp"), atol=1e-8)
    assert np.allclose(ints.T, pm.intor("int1e_kin"), atol=1e-8)
    assert np.allclose(ints.V, pm.intor("int1e_nuc"), atol=1e-8)
    assert np.allclose(ints.eri, pm.intor("int2e"), atol=1e-8)
    assert ints.e_nuc == pytest.approx(pm.energy_nuc(), abs=1e-10)
    mf = scf.RHF(pm)
    mf.conv_tol = 1e-12
    e_hf = mf.kernel()
    e_fci = fci.FCI(mf).kernel()[0]
    prob = qubit_problem(mol)
    assert prob.e_hf == pytest.approx(e_hf, abs=1e-8)
    assert prob.e_fci == pytest.approx(e_fci, abs=1e-8)


@pytest.mark.parametrize("name", MOLS)
def test_tapering_preserves_the_spectrum(name):
    prob = qubit_problem(molecule(name))
    H = prob.hamiltonian
    assert H.n == prob.full.n - 2
    ev = np.linalg.eigvalsh(H.matrix())
    assert np.min(np.abs(ev - prob.e_fci)) < 1e-9
    assert H.expectation_basis(prob.hf_bits).real == pytest.approx(prob.e_hf, abs=1e-9)
    assert prob.full.expectation_basis(prob.full_hf_bits).real == pytest.approx(prob.e_hf, abs=1e-9)
