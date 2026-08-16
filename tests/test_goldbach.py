"""Verification: every instrument is checked against independent ground truth."""

import numpy as np
import pytest

from goldbach.hardy_littlewood import hl_prediction, singular_series
from goldbach.lens import FareyLens, farey_fractions
from goldbach.partitions import goldbach_count_direct, goldbach_counts
from goldbach.sieve import primes_up_to
from goldbach.witnesses import first_witnesses, witness_records

# OEIS A045917: Goldbach partition counts g(2n) for 2n = 4, 6, 8, ...
A045917 = [1, 1, 1, 2, 1, 2, 2, 2, 2, 3, 3, 3, 2, 3, 2, 4, 4, 2, 3, 4, 3, 4,
           5, 4, 3, 5, 3, 4, 6, 3, 5, 6, 2, 5, 6, 5, 5, 7, 4, 5, 8, 5, 4, 9,
           4, 5, 7, 3, 6, 8, 5]

# OEIS A025018/A025019-style record pairs (n, least witness p), starting at 4.
KNOWN_RECORDS = [(4, 2), (6, 3), (12, 5), (30, 7), (98, 19), (220, 23),
                 (308, 31), (556, 47), (992, 73), (2642, 103), (5372, 139),
                 (7426, 173), (43532, 211), (54244, 233), (63274, 293),
                 (113672, 313), (128168, 331), (194428, 359), (194470, 383),
                 (413572, 389), (503222, 523)]


def test_counts_match_oeis():
    evens, g = goldbach_counts(104)
    assert list(g) == A045917


def test_fft_counts_match_brute_force():
    evens, g = goldbach_counts(2000)
    for n in [4, 10, 100, 666, 1000, 1998]:
        assert g[(n - 4) // 2] == goldbach_count_direct(n)


def test_singular_series_values():
    s = singular_series(30)
    evens = np.arange(4, 31, 2)
    # n = 4, 8, 16: no odd prime factor -> 1;  n = 6: factor 3 -> 2;
    # n = 30: factors 3 and 5 -> 2 * (4/3).
    lookup = dict(zip(evens.tolist(), s.tolist()))
    assert lookup[4] == 1 and lookup[8] == 1 and lookup[16] == 1
    assert lookup[6] == pytest.approx(2.0)
    assert lookup[30] == pytest.approx(2.0 * 4.0 / 3.0)


def test_hl_prediction_tracks_truth():
    limit = 200_000
    evens, g = goldbach_counts(limit)
    ev, pred = hl_prediction(limit)
    r2 = 2 * g[1:].astype(float)  # ordered counts (diagonal negligible here)
    ratio = r2[-5000:] / pred[-5000:]
    assert abs(np.mean(ratio) - 1.0) < 0.02  # HL is good to ~1% at 2e5


def test_lens_full_field_is_exact():
    lens = FareyLens(5000)
    field = lens.full_field()
    ps = primes_up_to(5000)
    for n in [100, 1234, 4444]:
        direct = sum(np.log(p) * np.log(n - p)
                     for p in ps[ps < n] if (n - p) in set(ps.tolist()))
        assert field[n] == pytest.approx(direct, rel=1e-9, abs=1e-6)


def test_lens_major_plus_minor_is_exact():
    lens = FareyLens(5000)
    total = lens.major_field(5, 16) + lens.minor_field(5, 16)
    assert np.allclose(total, lens.full_field(), atol=1e-6)


def test_farey_fractions_are_reduced_and_half():
    fracs = farey_fractions(7)
    assert (0, 1) in fracs and (1, 2) in fracs and (3, 7) in fracs
    assert (2, 4) not in fracs and (4, 7) not in fracs


def test_witness_records_match_known():
    evens, w = first_witnesses(600_000)
    recs = witness_records(evens, w)
    assert recs[: len(KNOWN_RECORDS)] == KNOWN_RECORDS
