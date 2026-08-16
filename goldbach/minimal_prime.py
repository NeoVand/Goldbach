"""The stubbornness measure: the smallest prime that Goldbach-splits n.

p_min(n) is the least prime p with n - p also prime.  It measures how hard
an even number resists: almost all n split with a single-digit prime, and
the record-setters — the most stubborn n — grow only polylogarithmically.
"""

import numpy as np

from goldbach.sieve import primes_up_to


def minimal_goldbach_prime(n_max: int, is_prime: np.ndarray) -> np.ndarray:
    """p_min(n) for every even n in [6, n_max], vectorised smallest-prime passes.

    Returns an int32 array pmin indexed so pmin[k] corresponds to n = 6 + 2k.
    Sweeps primes in increasing order, resolving in bulk the n with n - p prime;
    the unresolved set shrinks geometrically, so a few hundred passes suffice.
    Raises if any n up to n_max has no representation (a Goldbach failure).
    """
    evens = np.arange(6, n_max + 1, 2, dtype=np.int64)
    pmin = np.zeros(evens.size, dtype=np.int32)
    unresolved = np.arange(evens.size)
    for p in primes_up_to(n_max // 2):
        if p == 2:
            continue
        if unresolved.size == 0:
            break
        candidates = evens[unresolved]
        ok = (candidates >= 2 * p) & is_prime[candidates - p]
        pmin[unresolved[ok]] = p
        unresolved = unresolved[~ok]
    if unresolved.size:
        raise AssertionError(
            f"Goldbach failure(s) found at n = {evens[unresolved][:5]} — "
            "if this ever triggers on verified ranges, suspect the sieve."
        )
    return pmin


def records(pmin: np.ndarray, first_n: int = 6, step: int = 2):
    """Running-maximum records of p_min: list of (n, p_min(n)) where a new max is set."""
    previous_best = np.maximum.accumulate(np.concatenate(([0], pmin[:-1])))
    record_idx = np.nonzero(pmin > previous_best)[0]
    return [(int(first_n + step * k), int(pmin[k])) for k in record_idx]
