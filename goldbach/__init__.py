"""Instruments for exploring the structure of Goldbach partitions.

This package does NOT prove the Goldbach conjecture — nobody has.
It provides honest tools: exact bulk computation of partition counts,
the Hardy–Littlewood predicted density, and the spectral view of the
primes that links the two through the circle method.
"""

from goldbach.sieve import prime_sieve, primes_up_to
from goldbach.counts import goldbach_counts, unordered_counts
from goldbach.hardy_littlewood import TWIN_PRIME_C2, singular_series, li2, hl_prediction
from goldbach.minimal_prime import minimal_goldbach_prime, records
from goldbach.spectrum import prime_spectrum

__all__ = [
    "prime_sieve",
    "primes_up_to",
    "goldbach_counts",
    "unordered_counts",
    "TWIN_PRIME_C2",
    "singular_series",
    "li2",
    "hl_prediction",
    "minimal_goldbach_prime",
    "records",
    "prime_spectrum",
]
