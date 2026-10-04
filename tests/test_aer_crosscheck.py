"""ftsim against Qiskit Aer: noiseless amplitudes, and noisy outcome distributions under the
same Pauli channels (Aer exact density matrices vs our Monte Carlo, plain and stratified)."""

import numpy as np
import pytest

import ftalgo

qiskit = pytest.importorskip("qiskit")
pytest.importorskip("qiskit_aer")
from qiskit.quantum_info import Statevector  # noqa: E402
from scipy import stats  # noqa: E402

from ftalgo.qiskit_bridge import exact_distribution, to_qiskit  # noqa: E402

GATES1 = ["I", "X", "Y", "Z", "H", "S", "S_DAG", "T", "T_DAG", "SQRT_X", "SQRT_X_DAG"]
GATES2 = ["CX", "CY", "CZ", "SWAP"]
GATES3 = ["CCX", "CCZ"]


def random_gates(rng, n, count):
    out = []
    for _ in range(count):
        r = rng.random()
        if r < 0.55:
            g, ar = rng.choice(GATES1), 1
        elif r < 0.85:
            g, ar = rng.choice(GATES2), 2
        else:
            g, ar = rng.choice(GATES3), 3
        qs = rng.choice(n, size=ar, replace=False)
        out.append(f"{g} " + " ".join(map(str, qs)))
    return out


def random_noise(rng, n):
    kind = rng.integers(3)
    if kind == 0:
        p = rng.uniform(0, 0.05, size=3)
        q = rng.integers(n)
        return f"PAULI_CHANNEL_1({','.join(map(str, p))}) {q}"
    if kind == 1:
        p = rng.uniform(0, 0.05 / 15 * 2, size=15)
        a, b = rng.choice(n, size=2, replace=False)
        return f"PAULI_CHANNEL_2({','.join(map(str, p))}) {a} {b}"
    p = rng.uniform(0, 0.05 / 7 * 2, size=7)
    a, b, c = rng.choice(n, size=3, replace=False)
    return f"Z_CHANNEL_3({','.join(map(str, p))}) {a} {b} {c}"


@pytest.mark.parametrize("seed", range(100))
def test_noiseless_amplitudes(seed):
    rng = np.random.default_rng(seed)
    text = "\n".join(random_gates(rng, 6, 80))
    ours = ftalgo.Program(text + "\nI 5").statevector(backend="dense")
    qc, _ = to_qiskit(text + "\nI 5", measure=False)
    ref = np.asarray(Statevector(qc).data)
    assert np.allclose(ours, ref, atol=1e-10)
    sparse = ftalgo.Program(text + "\nI 5").statevector(backend="sparse")
    assert np.allclose(sparse, ref, atol=1e-10)


def noisy_program(seed):
    rng = np.random.default_rng(1000 + seed)
    lines = []
    for _ in range(8):
        lines += random_gates(rng, 4, 4)
        lines.append(random_noise(rng, 4))
    lines.append("M 0 1 2 3")
    return "\n".join(lines)


def histogram(records):
    keys = (records.astype(np.int64) << np.arange(records.shape[1])).sum(axis=1)
    return np.bincount(keys, minlength=1 << records.shape[1])


@pytest.mark.parametrize("seed", range(10))
def test_noisy_distribution_plain(seed):
    text = noisy_program(seed)
    exact = exact_distribution(text)
    probs = np.array([exact.get(k, 0.0) for k in range(16)])
    shots = 200_000
    counts = histogram(ftalgo.Program(text).sample(shots, seed=seed))
    keep = probs * shots >= 5
    f_exp = probs[keep] * shots
    f_obs = counts[keep].astype(float)
    # Fold the rare outcomes into one bin so the test keeps its degrees of freedom honest.
    f_exp = np.append(f_exp, shots - f_exp.sum())
    f_obs = np.append(f_obs, shots - f_obs.sum())
    if f_exp[-1] < 5:
        f_exp[-2] += f_exp[-1]
        f_obs[-2] += f_obs[-1]
        f_exp, f_obs = f_exp[:-1], f_obs[:-1]
    p = stats.chisquare(f_obs, f_exp).pvalue
    assert p > 1e-4, (p, text)


@pytest.mark.parametrize("seed", range(10))
def test_noisy_distribution_stratified(seed):
    text = noisy_program(seed)
    exact = exact_distribution(text)
    probs = np.array([exact.get(k, 0.0) for k in range(16)])
    prog = ftalgo.Program(text)
    kmax = 6
    pk = prog.fault_count_distribution(kmax)
    assert pk[-1] < 1e-5
    est = np.zeros(16)
    var = np.zeros(16)
    n = 40_000
    for k in range(kmax + 1):
        if pk[k] < 1e-9 or k > prog.num_sites:
            continue
        h = histogram(prog.sample(n, seed=seed * 100 + k, faults=k)) / n
        est += pk[k] * h
        var += pk[k] ** 2 * h * (1 - h) / n
    z = np.abs(est - probs) / np.sqrt(var + 1e-12)
    assert np.all((z < 5) | (np.abs(est - probs) < 2e-5)), (z, text)
