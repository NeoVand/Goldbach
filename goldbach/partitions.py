"""Exact Goldbach partition counts for every even number at once.

The number of ordered representations n = p + q (p, q prime) is the
self-convolution of the prime indicator, so a single FFT computes the whole
table. Values are integers of moderate size, and float64 FFT round-off is
~1e-9 here, so rounding recovers them exactly; `goldbach_count_direct`
provides an independent check.
"""

import numpy as np

from .sieve import prime_indicator


def ordered_representations(limit: int) -> np.ndarray:
    """R2[m] = #{(p, q) prime, ordered : p + q = m} for all m <= limit."""
    ind = prime_indicator(limit).astype(np.float64)
    m = 1 << int(np.ceil(np.log2(2 * limit + 1)))
    spectrum = np.fft.rfft(ind, n=m)
    conv = np.fft.irfft(spectrum * spectrum, n=m)[: limit + 1]
    return np.rint(conv).astype(np.int64)


def goldbach_counts(limit: int) -> tuple[np.ndarray, np.ndarray]:
    """Unordered Goldbach partition counts g(n) for all even n <= limit.

    Returns (evens, g) where evens = [4, 6, ..., limit'] and g[i] is the number
    of unordered prime pairs {p, q} with p + q = evens[i].
    """
    r2 = ordered_representations(limit)
    is_prime = prime_indicator(limit // 2)
    evens = np.arange(4, limit + 1, 2)
    diag = is_prime[evens // 2].astype(np.int64)  # p = q = n/2 contributes once
    g = (r2[evens] + diag) // 2
    return evens, g


def goldbach_count_direct(n: int) -> int:
    """Brute-force unordered count for a single even n (for verification)."""
    if n < 4 or n % 2:
        return 0
    ind = prime_indicator(n)
    ps = np.nonzero(ind[: n // 2 + 1])[0]
    return int(np.count_nonzero(ind[n - ps]))
