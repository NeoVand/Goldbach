"""The Hardy–Littlewood model of Goldbach counts.

The extended Goldbach conjecture (Hardy & Littlewood, 1923) predicts that
the ordered count g(n) for even n is asymptotic to

    g_HL(n) = 2 * C2 * S(n) * Li2(n),

where C2 is the twin prime constant, the singular series

    S(n) = prod over odd primes p | n of (p - 1)/(p - 2)

is the arithmetic correction for the residues n occupies, and

    Li2(n) = integral from 2 to n-2 of dt / (ln t * ln(n - t))

is the smooth density of prime pairs.  Everything here is a PREDICTION,
not a theorem — the experiments in this repo measure how well it holds.
"""

import numpy as np

from goldbach.sieve import primes_up_to

# Twin prime constant C2 = prod over odd primes p of (1 - 1/(p-1)^2).
TWIN_PRIME_C2 = 0.6601618158468695739278121100145557784326233602847334133194484233354056

# Number of Gauss–Legendre nodes for Li2; the integrand is smooth on [2, n-2].
_GL_NODES = 96


def singular_series(n_max: int) -> np.ndarray:
    """S(n) for every n in [0, n_max], by one multiplicative sweep per prime.

    S(n) = prod_{p | n, p odd prime} (p-1)/(p-2).  Convention: S(n) = 1 when n
    has no odd prime factor (powers of two), matching the HL formula for even n.
    """
    s = np.ones(n_max + 1, dtype=np.float64)
    for p in primes_up_to(n_max // 2):
        if p == 2:
            continue
        s[p::p] *= (p - 1) / (p - 2)
    # Odd primes larger than n_max // 2 divide no even n <= n_max, but do divide
    # themselves; only even-n values of S are used by the model, so leave them.
    return s


def li2(n: np.ndarray | float) -> np.ndarray:
    """Li2(n) = integral_2^{n-2} dt / (ln t * ln(n-t)), Gauss–Legendre quadrature.

    The integrand is symmetric about t = n/2, so fold to [2, n/2] and integrate
    in u = ln t, where the boundary layer near t = 2 flattens out and fixed-order
    quadrature converges far below the prime-fluctuation noise the model absorbs.
    """
    n = np.atleast_1d(np.asarray(n, dtype=np.float64))
    nodes, weights = np.polynomial.legendre.leggauss(_GL_NODES)
    lo, hi = np.log(2.0), np.log(n[:, None] / 2.0)
    half_span = (hi - lo) / 2.0
    u = lo + half_span * (nodes[None, :] + 1.0)
    t = np.exp(u)
    integrand = t / (u * np.log(n[:, None] - t))
    return 2.0 * (half_span[:, 0]) * (integrand * weights[None, :]).sum(axis=1)


def hl_prediction(even_n: np.ndarray, s_values: np.ndarray) -> np.ndarray:
    """Predicted ordered count g_HL(n) = 2 * C2 * S(n) * Li2(n) for even n.

    Li2 is evaluated on a dense logarithmic grid and interpolated in log n —
    it is smooth, so this is cheap and accurate for millions of n at once.
    """
    even_n = np.asarray(even_n, dtype=np.float64)
    grid = np.exp(np.linspace(np.log(6.0), np.log(even_n.max()), 4096))
    li2_grid = li2(grid)
    li2_n = np.interp(np.log(even_n), np.log(grid), li2_grid)
    return 2.0 * TWIN_PRIME_C2 * s_values * li2_n
