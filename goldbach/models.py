"""Random universes: Goldbach in ensembles of fake primes.

Three nested models, each keeping more of the primes' *local* (divisibility)
structure while staying random globally:

  naive   -- Cramer's model: k is "prime" with probability 1/log k.
  parity  -- odd k only, probability 2/log k (density-corrected).
  local m -- k coprime to m only, probability (m/phi(m))/log k.

The point: Goldbach-type statements hold with overwhelming probability in all
of them, and the local models even reproduce the banded structure of the real
Goldbach comet. What separates the real primes from these universes is only
the unproven assertion that primes behave "randomly enough" -- the models make
that gap quantitative.
"""

from math import gcd

import numpy as np


def survival_probabilities(limit: int, modulus: int = 1) -> np.ndarray:
    """pi[k] = model probability that k is 'prime', for k = 0..limit."""
    k = np.arange(limit + 1, dtype=np.float64)
    with np.errstate(divide="ignore"):
        pi = 1.0 / np.log(np.maximum(k, 3.0))
    if modulus > 1:
        phi = sum(1 for r in range(modulus) if gcd(r, modulus) == 1)
        coprime = np.array([gcd(int(x), modulus) == 1 for x in range(modulus)])
        pi *= modulus / phi
        pi[~coprime[np.arange(limit + 1) % modulus]] = 0.0
    pi[:3] = 0.0
    return np.clip(pi, 0.0, 1.0)


def sample_universe(limit: int, rng: np.random.Generator,
                    modulus: int = 1) -> np.ndarray:
    """One random universe: boolean 'prime' indicator under the model."""
    return rng.random(limit + 1) < survival_probabilities(limit, modulus)


def representation_counts(indicator: np.ndarray) -> np.ndarray:
    """Ordered two-'prime' representation counts of every n, via FFT."""
    n = len(indicator) - 1
    m = 1 << int(np.ceil(np.log2(2 * n + 2)))
    spec = np.fft.rfft(indicator.astype(np.float64), n=m)
    return np.rint(np.fft.irfft(spec * spec, n=m)[: n + 1]).astype(np.int64)


def failure_frequencies(limit: int, trials: int, rng: np.random.Generator,
                        modulus: int = 1, batch: int = 64) -> np.ndarray:
    """Empirical P(even n has NO two-'prime' representation), per even n >= 6.

    Returns an array over evens = 6, 8, ..., aligned with theory_curve.
    """
    m = 1 << int(np.ceil(np.log2(2 * limit + 2)))
    pi = survival_probabilities(limit, modulus)
    evens = np.arange(6, limit + 1, 2)
    fails = np.zeros(len(evens), dtype=np.int64)
    done = 0
    while done < trials:
        b = min(batch, trials - done)
        universes = rng.random((b, limit + 1)) < pi
        spec = np.fft.rfft(universes.astype(np.float64), n=m, axis=1)
        conv = np.fft.irfft(spec * spec, n=m, axis=1)[:, evens]
        fails += (conv < 0.5).sum(axis=0)
        done += b
    return fails / trials


def theory_curve(limit: int, modulus: int = 1) -> np.ndarray:
    """Exact model probability that even n has no representation, per even n.

    Independence makes this a closed product over unordered pairs:
    P(fail) = prod_{3 <= k <= n/2} (1 - pi_k * pi_{n-k}), with the k = n/2
    term being pi_{n/2} alone (a single site, not a pair).
    """
    pi = survival_probabilities(limit, modulus)
    evens = np.arange(6, limit + 1, 2)
    out = np.empty(len(evens))
    for i, n in enumerate(evens):
        ks = np.arange(3, n // 2)
        log_p = np.log1p(-pi[ks] * pi[n - ks]).sum()
        log_p += np.log1p(-pi[n // 2])
        out[i] = np.exp(log_p)
    return out
