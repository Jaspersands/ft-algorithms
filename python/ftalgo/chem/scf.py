"""Restricted Hartree–Fock and molecular-orbital integrals."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .basis import AOIntegrals


@dataclass
class RHF:
    energy: float          # total, including nuclear repulsion
    e_nuc: float
    C: np.ndarray          # MO coefficients (columns)
    eps: np.ndarray        # orbital energies
    n_occ: int
    iterations: int


def rhf(ints: AOIntegrals, n_electrons: int, *, tol: float = 1e-12, max_iter: int = 500) -> RHF:
    if n_electrons % 2:
        raise ValueError("closed-shell RHF needs an even electron count")
    n_occ = n_electrons // 2
    S, H, G = ints.S, ints.H, ints.eri
    w, U = np.linalg.eigh(S)
    X = U @ np.diag(w**-0.5) @ U.T  # symmetric orthogonalization

    def solve(F):
        e, Cp = np.linalg.eigh(X.T @ F @ X)
        return e, X @ Cp

    eps, C = solve(H)
    D = C[:, :n_occ] @ C[:, :n_occ].T
    e_old = None
    focks, errs = [], []
    for it in range(1, max_iter + 1):
        J = np.einsum("ijkl,kl->ij", G, D)
        K = np.einsum("ikjl,kl->ij", G, D)
        F = H + 2 * J - K
        e_el = float(np.sum(D * (H + F)))
        # DIIS on the commutator FDS − SDF.
        err = F @ D @ S - S @ D @ F
        focks.append(F)
        errs.append(err)
        if len(focks) > 8:
            focks.pop(0)
            errs.pop(0)
        if len(focks) >= 2:
            k = len(focks)
            B = -np.ones((k + 1, k + 1))
            B[k, k] = 0
            for a in range(k):
                for b in range(k):
                    B[a, b] = np.sum(errs[a] * errs[b])
            rhs = np.zeros(k + 1)
            rhs[k] = -1
            try:
                coef = np.linalg.solve(B, rhs)[:k]
                F = sum(c * f for c, f in zip(coef, focks))
            except np.linalg.LinAlgError:
                pass
        eps, C = solve(F)
        D_new = C[:, :n_occ] @ C[:, :n_occ].T
        if e_old is not None and abs(e_el - e_old) < tol and np.max(np.abs(D_new - D)) < 1e-9:
            D = D_new
            break
        D, e_old = D_new, e_el
    else:
        raise RuntimeError("RHF did not converge")
    # Final energy and orbitals from the converged density.
    J = np.einsum("ijkl,kl->ij", G, D)
    K = np.einsum("ikjl,kl->ij", G, D)
    F = H + 2 * J - K
    e_el = float(np.sum(D * (H + F)))
    eps, C = solve(F)
    return RHF(e_el + ints.e_nuc, ints.e_nuc, C, eps, n_occ, it)


def mo_integrals(ints: AOIntegrals, C: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(h_pq, (pq|rs)) in the MO basis."""
    h = C.T @ ints.H @ C
    g = np.einsum("pi,qj,rk,sl,pqrs->ijkl", C, C, C, C, ints.eri, optimize=True)
    return h, g
