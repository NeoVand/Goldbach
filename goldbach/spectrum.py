"""The spectrum of the primes: the circle method as a measurement.

The exponential sum S(alpha) = sum over primes p <= N of e^(2*pi*i*p*alpha)
is the heart of the circle method, because the Goldbach count is exactly

    g(n) = integral_0^1 |S(alpha)|^2 e^(-2*pi*i*n*alpha) d(alpha):

the comet is the Fourier transform of the primes' power spectrum.  On a
grid alpha = k/M, S is computable exactly with one length-M FFT of the
histogram of primes modulo M.  Choosing M with many small divisors puts
every rational a/q with q | M exactly on the grid, so the major-arc
spikes — height ~ pi(N) * mu(q)/phi(q) at a/q — appear at full height.
"""

import numpy as np

from goldbach.sieve import primes_up_to

# 720720 = 2^4 * 3^2 * 5 * 7 * 11 * 13 = lcm(1..16): every a/q with q <= 16
# lies on the grid, and the length is FFT-smooth.
DEFAULT_GRID = 720720


def prime_spectrum(n_max: int, grid: int = DEFAULT_GRID):
    """|S(k/grid)| for k = 0..grid-1, exactly, via one FFT of prime residue counts.

    Returns (alphas, magnitudes, prime_count).  S(k/M) depends on the primes
    only through their residues mod M, so binning then transforming is exact.
    """
    primes = primes_up_to(n_max)
    residue_counts = np.bincount(primes % grid, minlength=grid).astype(np.float64)
    spectrum = np.fft.fft(residue_counts)
    alphas = np.arange(grid) / grid
    return alphas, np.abs(spectrum), len(primes)
