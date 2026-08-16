"""Sanity tests: every instrument is checked against brute force on small ranges.

Run with `python tests/test_goldbach.py` (no pytest needed) or via pytest.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np

from goldbach.counts import goldbach_counts, unordered_counts
from goldbach.hardy_littlewood import TWIN_PRIME_C2, li2, singular_series
from goldbach.minimal_prime import minimal_goldbach_prime, records
from goldbach.sieve import prime_sieve, primes_up_to
from goldbach.spectrum import prime_spectrum


def brute_force_ordered(n, primes_set):
    return sum(1 for p in primes_set if p != 2 and (n - p) in primes_set and n - p != 2)


def test_sieve():
    assert list(primes_up_to(30)) == [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    assert prime_sieve(1)[0] == False and prime_sieve(1)[1] == False


def test_counts_match_brute_force():
    n_max = 3000
    g, is_prime, residual = goldbach_counts(n_max)
    assert residual < 0.01, f"FFT residual too large: {residual}"
    primes_set = set(primes_up_to(n_max).tolist())
    for n in range(6, n_max + 1, 2):
        assert g[n] == brute_force_ordered(n, primes_set), f"mismatch at n={n}"
    # Known values: 10 = 3+7 = 5+5 = 7+3 -> ordered 3; unordered 2.
    assert g[10] == 3
    r = unordered_counts(g, is_prime)
    assert r[10] == 2
    for n in range(6, n_max + 1, 2):
        diag = 1 if (n // 2 in primes_set and n // 2 != 2) else 0
        assert r[n] == (g[n] + diag) // 2


def test_no_goldbach_failures_small():
    g, _, _ = goldbach_counts(10**5)
    evens = np.arange(6, 10**5 + 1, 2)
    assert (g[evens] > 0).all()


def test_singular_series():
    s = singular_series(10**4)
    for n in (6, 30, 128, 9240):
        expected = 1.0
        for p in primes_up_to(n):
            if p > 2 and n % p == 0:
                expected *= (p - 1) / (p - 2)
        assert abs(s[n] - expected) < 1e-12, f"S({n}) wrong"
    assert s[128] == 1.0  # power of two: no odd prime factors


def test_li2_against_fine_trapezoid():
    # Independent reference: the integrand is symmetric about n/2, so
    # Li2 = 2 * integral_2^{n/2}; trapezoid on a dense geometric grid keeps
    # the reference accurate inside the boundary layer near t = 2.
    for n in (100.0, 10_000.0, 2_000_000.0):
        t = np.geomspace(2.0, n / 2.0, 4_000_001)
        ref = 2.0 * np.trapezoid(1.0 / (np.log(t) * np.log(n - t)), t)
        val = li2(np.array([n]))[0]
        assert abs(val - ref) / ref < 1e-6, f"Li2({n}): {val} vs {ref}"


def test_twin_prime_constant():
    # Recompute C2 = prod_{2 < p <= P} (1 - 1/(p-1)^2); tail beyond P is O(1/P).
    p = primes_up_to(10**7)[1:].astype(np.float64)
    approx = np.prod(1.0 - 1.0 / (p - 1.0) ** 2)
    assert abs(approx - TWIN_PRIME_C2) < 1e-7


def test_minimal_prime():
    n_max = 10**4
    is_prime = prime_sieve(n_max)
    pmin = minimal_goldbach_prime(n_max, is_prime)
    primes = primes_up_to(n_max)
    for k, n in enumerate(range(6, n_max + 1, 2)):
        expected = next(int(p) for p in primes if p > 2 and is_prime[n - p])
        assert pmin[k] == expected, f"p_min({n})"
    recs = records(pmin)
    assert recs[0] == (6, 3)
    assert all(b[1] > a[1] and b[0] > a[0] for a, b in zip(recs, recs[1:]))


def test_spectrum_peaks():
    # |S(a/q)| should be ~ pi(N)/phi(q) at low-denominator rationals and tiny
    # at "generic" alpha; check the q=1,2,3 spikes dominate.
    alphas, mags, n_primes = prime_spectrum(10**6, grid=720720)
    def at(frac):
        return mags[int(round(frac * len(alphas)))]
    assert at(0.0) == n_primes  # S(0) = pi(N)
    assert at(1 / 2) > 0.95 * n_primes  # almost all primes are odd
    assert at(1 / 3) > 0.45 * n_primes  # phi(3) = 2
    assert at(1 / 5) < 0.30 * n_primes  # phi(5) = 4
    generic = mags[12_345]  # far from any small-denominator rational
    assert generic < 0.01 * n_primes


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_"):
            fn()
            print(f"PASS {name}")
    print("all tests passed")
