"""Goldbach Lab: computational instruments for the Goldbach conjecture.

This package does not prove the conjecture (nobody has). It provides exact
large-scale computation of Goldbach partition counts, the Hardy-Littlewood
deconvolution, a spectral "Farey lens" implementation of the circle method,
random-model Monte Carlo, and first-witness record hunting.
"""

from .sieve import primes_up_to, prime_indicator
from .partitions import goldbach_counts, goldbach_count_direct
from .hardy_littlewood import singular_series, hl_integral, hl_prediction
from .lens import FareyLens
from .witnesses import first_witnesses, witness_records

__all__ = [
    "primes_up_to",
    "prime_indicator",
    "goldbach_counts",
    "goldbach_count_direct",
    "singular_series",
    "hl_integral",
    "hl_prediction",
    "FareyLens",
    "first_witnesses",
    "witness_records",
]
