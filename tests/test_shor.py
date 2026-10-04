import numpy as np
import pytest
from scipy import stats

import ftalgo
from ftalgo import shor


def test_ideal_distribution_n15():
    p = shor.ideal_distribution(15, 7, 8)
    assert p.sum() == pytest.approx(1.0)
    # r = 4 divides 2^8: all weight on multiples of 64.
    assert p[[0, 64, 128, 192]] == pytest.approx([0.25] * 4)
    p21 = shor.ideal_distribution(21, 2, 10)
    assert p21.sum() == pytest.approx(1.0)


def test_post_processing():
    assert shor.factors_from_outcome(64, 8, 15, 7) == (3, 5)
    assert shor.factors_from_outcome(0, 8, 15, 7) is None
    # 128/256 = 1/2 gives the candidate 2 then its multiple 4.
    assert shor.factors_from_outcome(128, 8, 15, 7) == (3, 5)
    assert shor.convergent_denominators(171, 1024, 21) == [1, 5, 6]


@pytest.mark.parametrize("N", shor.INSTANCES)
def test_bases(N):
    a = shor.choose_base(N)
    assert a == 2
    assert shor.is_good_base(a, N)


def histogram(y, m):
    return np.bincount(y.astype(np.int64), minlength=1 << m)


@pytest.mark.parametrize("N", [15, 21])
def test_exact_qft_matches_ideal_distribution(N):
    c, info = shor.order_finding(N, cutoff=None if False else 2 * N.bit_length())
    p = ftalgo.Program(c.to_text())
    shots = 4000
    y = info.outcomes(p.sample(shots, seed=1, faults="none"))
    ideal = shor.ideal_distribution(N, info.a, info.m)
    keep = ideal * shots >= 5
    f_exp = np.append(ideal[keep] * shots, shots * ideal[~keep].sum())
    h = histogram(y, info.m)
    f_obs = np.append(h[keep], h[~keep].sum()).astype(float)
    if f_exp[-1] < 5:
        f_exp[-2] += f_exp[-1]
        f_obs[-2] += f_obs[-1]
        f_exp, f_obs = f_exp[:-1], f_obs[:-1]
    assert stats.chisquare(f_obs, f_exp).pvalue > 1e-4


@pytest.mark.parametrize("N", [15, 21, 35])
def test_default_aqft_success_close_to_exact(N):
    c, info = shor.order_finding(N)
    p = ftalgo.Program(c.to_text())
    shots = 3000
    y = info.outcomes(p.sample(shots, seed=2, faults="none"))
    ok = shor.success_table(N, info.a, info.m)[y].mean()
    exact = shor.ideal_success(N, info.a, info.m)
    assert abs(ok - exact) < 4 * np.sqrt(exact * (1 - exact) / shots) + 0.01, (ok, exact)
