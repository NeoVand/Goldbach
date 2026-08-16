"""The Farey lens: the circle method as a computational instrument.

Weight each prime p <= N by log p and treat the result as a signal on the
integers. Its Fourier transform S(alpha) = sum_p log p * e(2 pi i p alpha)
concentrates at rational frequencies a/q with small q ("major arcs"); the
weighted Goldbach count R(n) = sum_{p+q=n} log p log q is the n-th Fourier
coefficient of S^2.

The lens filters the spectrum of S^2, keeping only bins within a half-width
of `bins` around every Farey fraction a/q with q <= Q, and inverse-transforms.
One pass reconstructs the major-arc part of R(n) for every n simultaneously.
What the filter discards is exactly the minor-arc contribution -- the term no
one can bound well enough to prove the conjecture. Here we can measure it.
"""

from math import gcd

import numpy as np

from .sieve import prime_indicator


def farey_fractions(q_max: int) -> list[tuple[int, int]]:
    """All reduced fractions a/q with 1 <= q <= q_max, 0 <= a/q <= 1/2.

    (Fractions in (1/2, 1) are conjugate mirrors and are handled implicitly
    by the real-signal symmetry of the rfft.)
    """
    fracs = [(0, 1)]
    for q in range(2, q_max + 1):
        for a in range(1, q // 2 + 1):
            if gcd(a, q) == 1 and 2 * a <= q:
                fracs.append((a, q))
    return fracs


class FareyLens:
    def __init__(self, n_max: int, fft_size: int | None = None):
        self.n_max = n_max
        m = fft_size or 1 << int(np.ceil(np.log2(2 * n_max + 2)))
        if m < 2 * n_max + 1:
            raise ValueError("fft_size must be >= 2*n_max + 1 (no wraparound)")
        self.m = m
        ind = prime_indicator(n_max)
        signal = np.zeros(n_max + 1)
        ks = np.nonzero(ind)[0]
        signal[ks] = np.log(ks)
        self.spectrum = np.fft.rfft(signal, n=m)  # S at alpha = j/m, j <= m/2
        self._sq = self.spectrum * self.spectrum  # spectrum of the convolution

    def magnitude(self, stride: int = 1) -> tuple[np.ndarray, np.ndarray]:
        """(alpha, |S(alpha)|) on the half-circle grid, for plotting."""
        mag = np.abs(self.spectrum[::stride])
        alpha = np.arange(0, self.m // 2 + 1, stride) / self.m
        return alpha, mag

    def full_field(self) -> np.ndarray:
        """Exact weighted count R(n) = sum_{p+q=n} log p log q for all n."""
        return np.fft.irfft(self._sq, n=self.m)[: 2 * self.n_max + 1]

    def arc_mask(self, q_max: int, bins: int) -> np.ndarray:
        """Boolean mask over rfft bins: within `bins` of some a/q, q <= q_max."""
        mask = np.zeros(self.m // 2 + 1, dtype=bool)
        half = self.m // 2
        for a, q in farey_fractions(q_max):
            center = round(a * self.m / q)
            lo, hi = max(0, center - bins), min(half, center + bins)
            mask[lo : hi + 1] = True
        return mask

    def major_field(self, q_max: int, bins: int) -> np.ndarray:
        """Major-arc reconstruction of R(n) for all n <= 2*n_max."""
        filtered = np.where(self.arc_mask(q_max, bins), self._sq, 0)
        return np.fft.irfft(filtered, n=self.m)[: 2 * self.n_max + 1]

    def minor_field(self, q_max: int, bins: int) -> np.ndarray:
        """The discarded minor-arc contribution: R - R_major."""
        return self.full_field() - self.major_field(q_max, bins)
