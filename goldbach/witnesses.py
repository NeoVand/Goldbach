"""Smallest Goldbach witness: the least prime p with n - p prime."""

import numpy as np

from .sieve import prime_indicator, primes_up_to


def first_witnesses(limit: int) -> tuple[np.ndarray, np.ndarray]:
    """For every even n in [4, limit], the least prime p with n - p prime.

    Returns (evens, witness). Raises if some even has no witness (which would
    be a counterexample to the conjecture -- report it loudly, don't crash).
    """
    is_prime = prime_indicator(limit)
    evens = np.arange(4, limit + 1, 2)
    witness = np.zeros(len(evens), dtype=np.int64)
    remaining = np.arange(len(evens))
    for p in primes_up_to(limit - 2):
        rem = evens[remaining] - p
        ok = (rem >= 2) & is_prime[rem]
        witness[remaining[ok]] = p
        remaining = remaining[~ok]
        if len(remaining) == 0:
            break
    if len(remaining):
        bad = evens[remaining]
        raise AssertionError(
            f"GOLDBACH COUNTEREXAMPLE CANDIDATE(S): {bad[:10]}..."
        )
    return evens, witness


def witness_records(evens: np.ndarray, witness: np.ndarray):
    """Running maxima of the first witness: the numbers hardest to represent."""
    idx = np.nonzero(witness > np.maximum.accumulate(
        np.concatenate([[0], witness[:-1]])))[0]
    return [(int(evens[i]), int(witness[i])) for i in idx]
