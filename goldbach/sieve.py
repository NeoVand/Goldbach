"""Vectorized sieve of Eratosthenes."""

import numpy as np


def prime_indicator(limit: int) -> np.ndarray:
    """Boolean array `is_prime` of length limit+1; is_prime[k] iff k is prime."""
    if limit < 2:
        return np.zeros(limit + 1, dtype=bool)
    is_prime = np.ones(limit + 1, dtype=bool)
    is_prime[:2] = False
    for p in range(2, int(limit**0.5) + 1):
        if is_prime[p]:
            is_prime[p * p :: p] = False
    return is_prime


def prime_indicator_odd_only(limit: int) -> np.ndarray:
    """Indicator of odd primes (2 excluded). Used where p=2 pairs are irrelevant."""
    ind = prime_indicator(limit)
    ind[2] = False
    return ind


def primes_up_to(limit: int) -> np.ndarray:
    """Array of all primes <= limit."""
    return np.nonzero(prime_indicator(limit))[0]
