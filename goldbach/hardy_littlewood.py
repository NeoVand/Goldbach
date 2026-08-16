"""The Hardy-Littlewood prediction for Goldbach counts.

The extended Goldbach conjecture (Hardy & Littlewood 1923) predicts that the
number of ordered representations of even n as p + q is asymptotic to

    R2(n) ~ 2 C2 * prod_{p | n, p odd} (p-1)/(p-2) * J(n),

where C2 = 0.660161... is the twin-prime constant and
J(n) = integral_2^{n-2} dt / (log t * log(n-t)).

`singular_series` computes the jagged arithmetic factor for every even n at
once with a sieve; `hl_integral` computes the smooth factor J(n) by quadrature
on a log-spaced grid plus spline interpolation.
"""

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import CubicSpline

from .sieve import primes_up_to

TWIN_PRIME_C2 = 0.6601618158468695739278121100145557784326


def singular_series(limit: int) -> np.ndarray:
    """prod_{p | n, p odd} (p-1)/(p-2) for even n; array indexed so that
    entry i corresponds to n = 4 + 2i (matching goldbach_counts' evens)."""
    evens = np.arange(4, limit + 1, 2)
    s = np.ones(len(evens), dtype=np.float64)
    for p in primes_up_to(limit // 2):
        if p == 2:
            continue
        # even multiples of p among evens: n = 2*p*k -> index (2*p*k - 4)/2
        first = 2 * p
        idx = np.arange((first - 4) // 2, len(s), p)
        s[idx] *= (p - 1) / (p - 2)
    return s


def _J(n: float) -> float:
    val, _ = quad(lambda t: 1.0 / (np.log(t) * np.log(n - t)), 2.0, n - 2.0,
                  limit=200)
    return val


def hl_integral(ns: np.ndarray, grid_points: int = 200) -> np.ndarray:
    """J(n) for an array of even n, via quadrature on a log grid + spline."""
    ns = np.asarray(ns, dtype=np.float64)
    lo, hi = ns.min(), ns.max()
    if lo < 6:
        raise ValueError("J(n) needs n >= 6")
    grid = np.exp(np.linspace(np.log(lo), np.log(hi), grid_points))
    vals = np.array([_J(g) for g in grid])
    spline = CubicSpline(np.log(grid), vals)
    return spline(np.log(ns))


def hl_prediction(limit: int) -> tuple[np.ndarray, np.ndarray]:
    """Predicted ordered representation count for every even n in [6, limit].

    Returns (evens, prediction) with evens = [6, 8, ..., limit'].
    """
    s = singular_series(limit)[1:]  # drop n = 4
    evens = np.arange(6, limit + 1, 2)
    j = hl_integral(evens)
    return evens, 2.0 * TWIN_PRIME_C2 * s * j
