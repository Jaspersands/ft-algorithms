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
    assert abs(evals_full[0] - prob.e_fci) < 1e-11

    # Spectrum of tapered Hamiltonian
    M_tap = prob.hamiltonian.matrix()
    evals_tap = scipy.linalg.eigvalsh(M_tap)
    assert abs(evals_tap[0] - prob.e_fci) < 1e-11


def test_lih_dissociation_curve_minimum():
    pts = lih_dissociation_curve([1.2, 1.4, 1.595, 1.8, 2.2])
    assert len(pts) == 5
    # The minimum of E_FCI and E_CASCI should be around 1.595 A
    energies = [p.e_casci for p in pts]
    min_idx = int(np.argmin(energies))
    assert min_idx in (2, 1)  # near 1.6 A
    assert pts[min_idx].e_casci < pts[0].e_casci
    assert pts[min_idx].e_casci < pts[-1].e_casci


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
