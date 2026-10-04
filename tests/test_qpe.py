import math

import numpy as np
import pytest
import scipy.linalg as sl
from scipy import stats

import ftalgo
from ftalgo import qpe
from ftalgo.chem.qubit import PauliSum, molecule, qubit_problem


@pytest.fixture(scope="module")
def h2():
    return qubit_problem(molecule("H2"))


@pytest.mark.parametrize("order", [2, 4])
def test_trotter_converges(h2, order):
    H = h2.hamiltonian
    Hm = H.matrix() - h2.e_hf * np.eye(4)
    exact = sl.expm(-1j * Hm * 0.5)
    errs = []
    for r in (4, 8, 16):
        U = np.linalg.matrix_power(qpe.trotter_unitary(H, 0.5 / r, h2.e_hf, order), r)
        errs.append(np.linalg.norm(U - exact, 2))
    # Error ∝ r^(−order): halving dt divides it by ~2^order.
    assert errs[0] / errs[1] == pytest.approx(2**order, rel=0.15)
    assert errs[1] / errs[2] == pytest.approx(2**order, rel=0.15)


def test_controlled_rotation_matrix():
    from ftalgo.circuit import Circuit

    H = PauliSum(2, {(0b01, 0b11): 1.0})  # Y on qubit 0, Z on qubit 1
    theta = 0.37
    for ctrl_val in (0, 1):
        c = Circuit(3)
        if ctrl_val:
            c.x(0)
        c.h(1)
        c.t(1)
        c.h(2)
        prep = Circuit(3)
        prep.extend(c)
        qpe.controlled_pauli_rotation(c, 0, [1, 2], (0b01, 0b11), theta, 1e-10)
        v = ftalgo.Program(c.to_text()).statevector()
        w = ftalgo.Program(prep.to_text()).statevector()
        P = PauliSum(3, {(0b010, 0b110): 1.0}).matrix()
        R = np.cos(theta) * np.eye(8) - 1j * np.sin(theta) * P
        want = R @ w if ctrl_val else w
        assert abs(abs(np.vdot(want, v)) - 1) < 1e-9


def test_noiseless_qpe_matches_ideal_distribution(h2):
    c, info = qpe.qpe_circuit(h2.hamiltonian, h2.hf_bits, m=5, tau=2 * math.pi, steps=2, e_ref=h2.e_hf, order=2, eps=1e-6)
    p = ftalgo.Program(c.to_text())
    shots = 6000
    y = qpe.estimate(p.sample(shots, seed=3, faults="none"), info.records)
    ideal = qpe.ideal_distribution(h2.hamiltonian, h2.hf_bits, info)
    h = np.bincount(y, minlength=32).astype(float)
    keep = ideal * shots >= 5
    f_exp = np.append(ideal[keep] * shots, ideal[~keep].sum() * shots)
    f_obs = np.append(h[keep], h[~keep].sum())
    if f_exp[-1] < 5:
        f_exp[-2] += f_exp[-1]
        f_obs[-2] += f_obs[-1]
        f_exp, f_obs = f_exp[:-1], f_obs[:-1]
    assert stats.chisquare(f_obs, f_exp).pvalue > 1e-4


def test_h2_chemical_accuracy_noiseless(h2):
    """The production parameters reach chemical accuracy in most noiseless runs."""
    c, info = qpe.qpe_circuit(h2.hamiltonian, h2.hf_bits, m=10, tau=2 * math.pi, steps=4, e_ref=h2.e_hf, order=4)
    ideal = qpe.ideal_distribution(h2.hamiltonian, h2.hf_bits, info)
    y = np.arange(1 << info.m)
    phi = np.where(y / 1024 >= 0.5, y / 1024 - 1, y / 1024)
    E = info.e_ref - 2 * math.pi * phi / info.tau
    p_ok = ideal[np.abs(E - h2.e_fci) < qpe.CHEMICAL_ACCURACY].sum()
    assert p_ok > 0.85
    rows = ftalgo.Program(c.to_text()).sample(40, seed=4, faults="none")
    got = np.abs(info.energies(rows) - h2.e_fci) < qpe.CHEMICAL_ACCURACY
    assert got.mean() > 0.6
