"""The comet camera: every Goldbach count at once, from a single FFT.

The number of ordered ways to write n as a sum of two odd primes is the
autocorrelation of the prime indicator:

    g(n) = sum_{i+j=n} a(i) a(j),   a = indicator of the odd primes.

Autocorrelation is multiplication in Fourier space, so one FFT of the
prime indicator "photographs" g(n) for every n <= N simultaneously in
O(N log N) — no per-n work at all.  This is also exactly the circle
method made executable: the FFT of a IS the exponential sum S(alpha)
sampled on a grid, and g(n) is the inverse transform of |S|^2.
"""

import numpy as np

from goldbach.sieve import prime_sieve


def goldbach_counts(n_max: int):
    """Ordered Goldbach counts for every n up to n_max, by FFT autocorrelation.

    Returns (g, is_prime, residual) where
      g[n]      = #{(p, q) : p + q = n, p and q odd primes}  (ordered pairs;
                  the diagonal p = q = n/2 contributes 1),
      is_prime  = boolean sieve up to n_max,
      residual  = max |raw FFT value - nearest integer| over all n, the
                  floating-point safety margin (must be far below 0.5 for
                  the rounded counts to be exact).
    """
    if n_max < 6:
        raise ValueError("n_max must be at least 6")
    is_prime = prime_sieve(n_max)
    a = is_prime.astype(np.float64)
    a[2] = 0.0  # odd primes only; 4 = 2 + 2 is handled by callers explicitly

    # Linear (non-circular) convolution needs transform length >= 2*n_max + 1.
    length = 1 << (2 * n_max + 1).bit_length()
    fa = np.fft.rfft(a, length)
    conv = np.fft.irfft(fa * fa, length)[: n_max + 1]

    g = np.rint(conv).astype(np.int64)
    residual = float(np.max(np.abs(conv - g)))
    return g, is_prime, residual


def unordered_counts(g: np.ndarray, is_prime: np.ndarray) -> np.ndarray:
    """Unordered counts r(n) = #{{p, q} : p + q = n, p <= q} from ordered g(n).

    Ordered pairs double-count {p, q} with p != q and count the diagonal once:
    r(n) = (g(n) + [n/2 is an odd prime]) / 2.
    """
    n_max = len(g) - 1
    diag = np.zeros(n_max + 1, dtype=np.int64)
    half = np.arange(0, n_max + 1, 2) // 2
    diag[::2] = is_prime[half] & (half != 2)
    return (g + diag) // 2
