"""Qubit Hamiltonians: Pauli algebra, Jordan–Wigner, and Z2 tapering (Bravyi, Gambetta,
Mezzacapo and Temme 2017).

A Pauli string is a pair of bit masks (x, z): qubit q carries X if only x has bit q, Z if only z
does, and Y if both. ``PauliSum`` maps (x, z) → complex coefficient.

Spin orbitals are interleaved: orbital i has spin-up mode 2i and spin-down mode 2i + 1.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .basis import STO3G_COEFFS, SFunction, ao_integrals, sto3g_exponents
from .scf import mo_integrals, rhf

ANGSTROM = 1.8897261254578281  # bohr per ångström (CODATA 2018)

# Single-qubit products: index 0 I, 1 X, 2 Y, 3 Z; MUL[a][b] = (phase, c) with a·b = phase·c.
_MUL = [
    [(1, 0), (1, 1), (1, 2), (1, 3)],
    [(1, 1), (1, 0), (1j, 3), (-1j, 2)],
    [(1, 2), (-1j, 3), (1, 0), (1j, 1)],
    [(1, 3), (1j, 2), (-1j, 1), (1, 0)],
]


def _letter(x: int, z: int, q: int) -> int:
    xb, zb = x >> q & 1, z >> q & 1
    return {(0, 0): 0, (1, 0): 1, (1, 1): 2, (0, 1): 3}[(xb, zb)]


def _bits(letter: int) -> tuple[int, int]:
    return [(0, 0), (1, 0), (1, 1), (0, 1)][letter]


def pauli_mul(a: tuple[int, int], b: tuple[int, int], n: int) -> tuple[complex, tuple[int, int]]:
    ph = 1
    x = z = 0
    for q in range(n):
        p, c = _MUL[_letter(*a, q)][_letter(*b, q)]
        ph *= p
        cx, cz = _bits(c)
        x |= cx << q
        z |= cz << q
    return ph, (x, z)


def pauli_label(key: tuple[int, int], n: int) -> str:
    """Qubit 0 first: e.g. 'XZI'."""
    return "".join("IXYZ"[_letter(*key, q)] for q in range(n))


class PauliSum:
    def __init__(self, n: int, terms: dict | None = None):
        self.n = n
        self.terms: dict[tuple[int, int], complex] = dict(terms or {})

    @staticmethod
    def identity(n: int, c: complex = 1.0) -> "PauliSum":
        return PauliSum(n, {(0, 0): c})

    def __add__(self, o: "PauliSum") -> "PauliSum":
        t = dict(self.terms)
        for k, v in o.terms.items():
            t[k] = t.get(k, 0) + v
        return PauliSum(self.n, t)

    def scale(self, c: complex) -> "PauliSum":
        return PauliSum(self.n, {k: v * c for k, v in self.terms.items()})

    def __mul__(self, o: "PauliSum") -> "PauliSum":
        t: dict = {}
        for ka, va in self.terms.items():
            for kb, vb in o.terms.items():
                ph, k = pauli_mul(ka, kb, self.n)
                t[k] = t.get(k, 0) + ph * va * vb
        return PauliSum(self.n, t)

    def simplify(self, tol: float = 1e-12) -> "PauliSum":
        return PauliSum(self.n, {k: v for k, v in self.terms.items() if abs(v) > tol})

    def real(self) -> "PauliSum":
        for v in self.terms.values():
            if abs(complex(v).imag) > 1e-10:
                raise ValueError("non-Hermitian term")
        return PauliSum(self.n, {k: float(complex(v).real) for k, v in self.terms.items()})

    def matrix(self) -> np.ndarray:
        """Dense 2^n matrix; qubit q is bit q of the basis index."""
        dim = 1 << self.n
        M = np.zeros((dim, dim), dtype=complex)
        idx = np.arange(dim)
        for (x, z), c in self.terms.items():
            # P|b⟩ = i^{|x∧z|} (−1)^{|z∧b|} |b ⊕ x⟩ with Y = iXZ.
            ny = bin(x & z).count("1")
            sign = np.array([(-1) ** bin(z & int(b)).count("1") for b in idx])
            M[idx ^ x, idx] += c * (1j**ny) * sign
        return M

    def expectation_basis(self, b: int) -> complex:
        """⟨b|H|b⟩ for a computational basis state."""
        return sum(c * (-1) ** bin(z & b).count("1") for (x, z), c in self.terms.items() if x == 0)

    def __len__(self) -> int:
        return len(self.terms)


def jw_annihilation(p: int, n: int) -> PauliSum:
    """a_p = Z_0 … Z_{p−1} (X_p + iY_p)/2."""
    zs = (1 << p) - 1
    return PauliSum(n, {(1 << p, zs): 0.5, (1 << p, zs | 1 << p): 0.5j})


def jw_creation(p: int, n: int) -> PauliSum:
    zs = (1 << p) - 1
    return PauliSum(n, {(1 << p, zs): 0.5, (1 << p, zs | 1 << p): -0.5j})


def fermion_hamiltonian(h: np.ndarray, g: np.ndarray, e_nuc: float) -> PauliSum:
    """H = E_nuc + Σ h_pq a†_pσ a_qσ + ½ Σ (pq|rs) a†_pσ a†_rτ a_sτ a_qσ under Jordan–Wigner."""
    norb = h.shape[0]
    n = 2 * norb
    a = [jw_annihilation(p, n) for p in range(n)]
    ad = [jw_creation(p, n) for p in range(n)]
    H = PauliSum.identity(n, e_nuc)
    for p in range(norb):
        for q in range(norb):
            if abs(h[p, q]) < 1e-14:
                continue
            for s in range(2):
                H = H + (ad[2 * p + s] * a[2 * q + s]).scale(h[p, q])
    for p in range(norb):
        for q in range(norb):
            for r in range(norb):
                for s_ in range(norb):
                    v = g[p, q, r, s_]
                    if abs(v) < 1e-14:
                        continue
                    for sg in range(2):
                        for tu in range(2):
                            P, Q, R, S = 2 * p + sg, 2 * q + sg, 2 * r + tu, 2 * s_ + tu
                            if P == R or Q == S:
                                continue
                            H = H + (ad[P] * ad[R] * a[S] * a[Q]).scale(0.5 * v)
    return H.simplify().real()


def taper(H: PauliSum, hf_bits: int) -> tuple[PauliSum, int, list[dict]]:
    """Removes the qubits of the spin-up and spin-down parity symmetries.

    τ = Π Z over one spin's modes commutes with H. The Clifford U = (X_q + τ)/√2 (q the
    symmetry's highest mode) maps τ ↦ X_q and leaves every term with I or X on q; X_q is replaced
    by τ's eigenvalue on the Hartree–Fock state and qubit q is dropped. The HF basis state maps
    to itself with bit q removed. Returns (tapered H, tapered HF bits, the symmetries used)."""
    n = H.n
    syms = []
    for spin in (1, 0):  # the higher pivot first, so the lower index is unaffected
        tau_z = sum(1 << q for q in range(spin, n, 2))
        pivot = max(q for q in range(spin, n, 2))
        sector = (-1) ** bin(hf_bits & tau_z).count("1")
        syms.append({"spin": "down" if spin else "up", "pivot": pivot, "sector": sector})
        new: dict = {}
        for key, c in H.terms.items():
            x, z = key
            # Check commutation with τ (Z-string): anticommute iff |x ∧ τ| odd.
            if bin(x & tau_z).count("1") % 2:
                raise ValueError("term does not commute with the symmetry")
            # Anticommutes with X_pivot iff z has the pivot bit.
            if z >> pivot & 1:
                # U P U = −P · X_q · τ
                ph1, k1 = pauli_mul(key, (1 << pivot, 0), n)
                ph2, k2 = pauli_mul(k1, (0, tau_z), n)
                c = -c * ph1 * ph2
                key = k2
            x, z = key
            assert not z >> pivot & 1
            if x >> pivot & 1:
                c = c * sector
                x &= ~(1 << pivot)
            # Drop qubit `pivot`: shift higher bits down.
            lo = (1 << pivot) - 1
            x = (x & lo) | ((x >> (pivot + 1)) << pivot)
            z = (z & lo) | ((z >> (pivot + 1)) << pivot)
            new[(x, z)] = new.get((x, z), 0) + c
        lo = (1 << pivot) - 1
        hf_bits = (hf_bits & lo) | ((hf_bits >> (pivot + 1)) << pivot)
        n -= 1
        H = PauliSum(n, new).simplify().real()
    return H, hf_bits, syms


@dataclass
class Molecule:
    name: str
    atoms: list  # [(element, (x, y, z) in bohr)]
    charge: int
    zetas: dict | None = None  # element → ζ (None: Basis Set Exchange STO-3G)

    @property
    def n_electrons(self) -> int:
        from .basis import CHARGE

        return sum(CHARGE[e] for e, _ in self.atoms) - self.charge

    def basis(self) -> list[SFunction]:
        out = []
        for el, c in self.atoms:
            z = None if self.zetas is None else self.zetas.get(el)
            out.append(SFunction(tuple(c), sto3g_exponents(el, z), STO3G_COEFFS))
        return out


def molecule(name: str, R: float | None = None, *, szabo: bool = False) -> Molecule:
    """Built-in molecules (distances in bohr). Defaults: H2 at 1.4 bohr (Szabo & Ostlund's
    geometry, ≈ 0.741 Å), HeH+ at 1.4632 bohr, H3+ an equilateral triangle of side 1.65 bohr,
    H4 a linear chain with 1.8 bohr spacing. ``szabo=True`` uses Szabo & Ostlund's ζ values
    (H 1.24, He 2.0925) instead of the Basis Set Exchange STO-3G."""
    zetas = {"H": 1.24, "He": 2.0925} if szabo else None
    if name == "H2":
        R = 1.4 if R is None else R
        return Molecule(name, [("H", (0.0, 0.0, 0.0)), ("H", (0.0, 0.0, R))], 0, zetas)
    if name == "HeH+":
        R = 1.4632 if R is None else R
        return Molecule(name, [("He", (0.0, 0.0, 0.0)), ("H", (0.0, 0.0, R))], 1, zetas)
    if name == "H3+":
        R = 1.65 if R is None else R
        h = R * math.sqrt(3) / 2
        return Molecule(name, [("H", (0.0, 0.0, 0.0)), ("H", (R, 0.0, 0.0)), ("H", (R / 2, h, 0.0))], 1, zetas)
    if name == "H4":
        R = 1.8 if R is None else R
        return Molecule(name, [("H", (0.0, 0.0, i * R)) for i in range(4)], 0, zetas)
    raise ValueError(f"unknown molecule {name}")


@dataclass
class QubitProblem:
    molecule: Molecule
    e_hf: float
    e_fci: float
    hamiltonian: PauliSum   # tapered
    hf_bits: int            # tapered HF basis state
    full: PauliSum          # untapered JW
    full_hf_bits: int
    symmetries: list


def number_sector_ground(H: PauliSum, n_up: int, n_down: int) -> float:
    """Lowest eigenvalue of a JW Hamiltonian within fixed spin-up and spin-down counts."""
    n = H.n
    idx = [b for b in range(1 << n) if bin(b & sum(1 << q for q in range(0, n, 2))).count("1") == n_up
           and bin(b & sum(1 << q for q in range(1, n, 2))).count("1") == n_down]
    M = H.matrix()[np.ix_(idx, idx)]
    return float(np.linalg.eigvalsh(M)[0])


def qubit_problem(mol: Molecule) -> QubitProblem:
    ints = ao_integrals(mol.atoms, mol.basis())
    ne = mol.n_electrons
    scf = rhf(ints, ne)
    h, g = mo_integrals(ints, scf.C)
    full = fermion_hamiltonian(h, g, ints.e_nuc)
    hf = (1 << ne) - 1  # lowest ne spin orbitals (interleaved: both spins of the lowest orbitals)
    e_fci = number_sector_ground(full, ne // 2, ne // 2)
    tapered, hf_t, syms = taper(full, hf)
    return QubitProblem(mol, scf.energy, e_fci, tapered, hf_t, full, hf, syms)
