"""The Lambda-weighted comet camera, for listening to the zeta zeros.

Analytic number theory weighs primes by the von Mangoldt function
Lambda(p^k) = log p, because that is the weighting the Riemann zeta function
speaks in.  The weighted Goldbach count

    G(n) = sum_{a+b=n} Lambda(a) Lambda(b)

obeys Fujii's explicit formula (1991): under RH,

    sum_{n<=X} G(n) = X^2/2 - 2 sum_rho X^(rho+1)/(rho(rho+1)) + smaller,

so the deviation of the cumulative count from X^2/2, rescaled by X^(3/2),
oscillates in u = log X at angular frequencies equal to the imaginary parts
of the zeta zeros, each with amplitude 4/|rho(rho+1)|.  The zeros are IN the
Goldbach data; this module computes the signal that carries them.
"""

import numpy as np

from goldbach.sieve import primes_up_to

# Imaginary parts of the first ten nontrivial zeta zeros (Odlyzko's tables).
ZETA_ZEROS = np.array([
    14.134725141734693, 21.022039638771555, 25.010857580145688,
    30.424876125859513, 32.935061587739190, 37.586178158825671,
    40.918719012147495, 43.327073280914999, 48.005150881167160,
    49.773832477672302,
])


def mangoldt(n_max: int) -> np.ndarray:
    """Von Mangoldt function Lambda(k) for k in [0, n_max]."""
    lam = np.zeros(n_max + 1, dtype=np.float64)
    for p in primes_up_to(n_max):
        log_p = np.log(float(p))
        pk = int(p)
        while pk <= n_max:
            lam[pk] = log_p
            pk *= int(p)
    return lam


def weighted_counts(weights: np.ndarray) -> np.ndarray:
    """G(n) = sum_{a+b=n} w(a) w(b) for all n, by FFT autocorrelation.

    Same instrument as goldbach.counts, but for real-valued weights (von
    Mangoldt, or a synthetic control set); the result is a float signal.
    """
    n_max = len(weights) - 1
    length = 1 << (2 * n_max + 1).bit_length()
    fw = np.fft.rfft(weights, length)
    return np.fft.irfft(fw * fw, length)[: n_max + 1]


def cumulative_deviation(g_weighted: np.ndarray):
    """D(X) = (sum_{n<=X} G(n) - X^2/2) / X^(3/2) for every X up to n_max.

    Under RH this equals -4 sum over zeros gamma of
    cos(gamma log X - phase) / |rho(rho+1)| up to o(1): a superposition of
    log-periodic tones, one per zeta zero.  Returns (X values, D values)
    starting at X = 2.
    """
    cum = np.cumsum(g_weighted)
    x = np.arange(2, len(g_weighted), dtype=np.float64)
    d = (cum[2:] - x * x / 2.0) / x**1.5
    return x, d


def log_periodogram(x: np.ndarray, d: np.ndarray, x_min: float, x_max: float,
                    samples: int = 4096, detrend_degree: int = 2):
    """Amplitude spectrum of D as a signal in u = log X, Hann-windowed.

    Resamples D uniformly in u on [log x_min, log x_max], removes a smooth
    polynomial trend (secular drift from lower-order terms), applies a Hann
    window, and returns (angular frequencies, amplitudes) calibrated so that
    a pure tone A*cos(gamma*u + phase) appears with height ~ A at its peak.
    """
    u = np.linspace(np.log(x_min), np.log(x_max), samples)
    d_u = d[np.minimum(np.exp(u).astype(np.int64), len(d) + 1) - 2]
    trend = np.polynomial.Polynomial.fit(u, d_u, detrend_degree)
    signal = d_u - trend(u)
    window = np.hanning(samples)
    spectrum = np.fft.rfft(signal * window)
    du = u[1] - u[0]
    freqs = 2.0 * np.pi * np.fft.rfftfreq(samples, du)
    amps = 2.0 * np.abs(spectrum) / window.sum()
    return freqs, amps


def predicted_amplitudes(zeros: np.ndarray = ZETA_ZEROS) -> np.ndarray:
    """Fujii amplitude 4/|rho(rho+1)| for each zero rho = 1/2 + i*gamma."""
    rho = 0.5 + 1j * zeros
    return 4.0 / np.abs(rho * (rho + 1.0))
