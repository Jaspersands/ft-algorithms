import numpy as np
import scipy.linalg
from ftalgo.chem.lih import lih_dissociation_curve, lih_problem
from ftalgo.qpe import CHEMICAL_ACCURACY, qpe_circuit


def test_lih_problem_matches_pyscf():
    prob = lih_problem(1.595, ncas=2)
    assert prob.hamiltonian.n == 2
    assert len(prob.hamiltonian.terms) > 0

    # Spectrum of full untapered Hamiltonian
    M_full = prob.full.matrix()
    evals_full = scipy.linalg.eigvalsh(M_full)
    assert abs(evals_full[0] - prob.e_fci) < 1e-9

    # Spectrum of tapered Hamiltonian
    M_tap = prob.hamiltonian.matrix()
    evals_tap = scipy.linalg.eigvalsh(M_tap)
    assert abs(evals_tap[0] - prob.e_fci) < 1e-9


def test_lih_rhf_orbitals_match_pyscf_casci():
    prob = lih_problem(1.595, ncas=2, orbitals="rhf")
    assert abs(scipy.linalg.eigvalsh(prob.hamiltonian.matrix())[0] - prob.e_fci) < 1e-11
    # Canonical-RHF CAS(2,2) is nearly uncorrelated: under 1 mHa below RHF.
    assert 0 < prob.e_hf - prob.e_fci < 1e-3


def test_lih_dissociation_curve_casscf_tracks_fci():
    pts = lih_dissociation_curve([1.2, 1.4, 1.595, 1.8, 2.2, 3.0])
    assert len(pts) == 6
    for p in pts:
        # Variational ordering: FCI <= CASSCF <= CASCI(RHF orbitals) <= RHF.
        assert p.e_fci <= p.e_casscf + 1e-9 <= p.e_casci + 2e-9 <= p.e_hf + 3e-9
        # CASSCF(2,2) stays within ~2 mHa of full STO-3G FCI along the curve.
        assert p.e_casscf - p.e_fci < 2.5e-3
    energies = [p.e_casscf for p in pts]
    min_idx = int(np.argmin(energies))
    assert min_idx in (1, 2)  # minimum near 1.6 Å



def test_lih_qpe_circuit_construction():
    prob = lih_problem(1.595, ncas=2)
    c, info = qpe_circuit(
        prob.hamiltonian,
        prob.hf_bits,
        m=10,
        tau=2 * np.pi,
        steps=4,
        e_ref=prob.e_hf,
        order=4,
    )
    assert c.num_qubits == 3  # 1 control + 2 system
    assert info.m == 10
    assert info.steps == 4
    assert info.order == 4
    assert len(info.records) == 10
    assert c.t_count() > 0
