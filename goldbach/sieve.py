"""Sieve of Eratosthenes, vectorised."""

import numpy as np


def prime_sieve(n: int) -> np.ndarray:
    """Boolean array `is_prime` of length n+1: is_prime[k] iff k is prime."""
    if n < 1:
        return np.zeros(n + 1, dtype=bool)
    is_prime = np.ones(n + 1, dtype=bool)
    is_prime[:2] = False
    for p in range(2, int(n**0.5) + 1):
        if is_prime[p]:
            is_prime[p * p :: p] = False
    return is_prime


def primes_up_to(n: int) -> np.ndarray:
    """All primes <= n as an int64 array."""
    return np.nonzero(prime_sieve(n))[0].astype(np.int64)
