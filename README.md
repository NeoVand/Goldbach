# Goldbach Lab

An experimental research program for the Goldbach conjecture: every even
number greater than 2 is the sum of two primes.

**This repository does not prove the conjecture, and does not pretend to.**
The conjecture is open. What lives here instead is a set of computational
instruments — built from scratch, verified against independent ground truth —
that make the conjecture's structure *visible* and its difficulty *measurable*:

| Instrument | What it does |
|---|---|
| `goldbach/partitions.py` | Exact partition counts g(n) for **every** even n ≤ 10⁷ in one FFT |
| `goldbach/hardy_littlewood.py` | The Hardy–Littlewood prediction; deconvolves the "Goldbach comet" |
| `goldbach/lens.py` | **The Farey lens** — the circle method as an instrument: filter the primes' Fourier spectrum to Farey fractions and reconstruct all Goldbach counts from arithmetic structure alone |
| `goldbach/models.py` | Random "fake prime" universes (Cramér-type), for measuring how overwhelmingly a random universe forces Goldbach |
| `goldbach/witnesses.py` | Least-witness records: the evens that resist longest |

The findings, figures, and an honest account of what a proof would require are
in **[REPORT.md](REPORT.md)**.

## Run it

```bash
pip install numpy scipy matplotlib pytest
python -m pytest tests/          # every instrument checked against ground truth
python experiments/exp1_comet.py       # ~2 min: comet + deconvolution at 10^7
python experiments/exp2_lens.py        # ~1 min: spectrum, Farey reconstruction, minor arcs
python experiments/exp3_universes.py   # ~2 min: 200,000 random universes
python experiments/exp4_witnesses.py   # ~1 min: witness records to 10^7
```

Figures land in `figures/`, numeric summaries in `data/`.

## Verification

`tests/test_goldbach.py` checks the FFT counter against OEIS A045917 and
brute force, the singular series against hand-computed values, the lens
against direct weighted counts (major + minor = exact, to 10⁻⁶), and the
witness records against OEIS A025018/A025019.
