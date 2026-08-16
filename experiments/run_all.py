"""Run every experiment, write figures/ and results/summary.json.

Usage: python experiments/run_all.py [N]     (default N = 20_000_000)
"""

import json
import sys
import time
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import numpy as np
import matplotlib.pyplot as plt

from experiments.style import (
    AQUA, BASELINE, BLUE, GRID, INK, INK_2, MUTED, ORANGE, finish, thousands,
)
from goldbach.counts import goldbach_counts, unordered_counts
from goldbach.hardy_littlewood import TWIN_PRIME_C2, li2, singular_series
from goldbach.minimal_prime import minimal_goldbach_prime, records
from goldbach.spectrum import prime_spectrum

FIG = ROOT / "figures"
RES = ROOT / "results"
COMET_LIMIT = 200_000  # plotting window for the comet figures


def comet_figure(g):
    n = np.arange(6, COMET_LIMIT + 1, 2)
    y = g[n]
    div3, div5 = (n % 3 == 0), (n % 5 == 0)
    classes = [
        (~div3, BLUE, "3 ∤ n"),
        (div3 & ~div5, ORANGE, "3 | n,  5 ∤ n"),
        (div3 & div5, AQUA, "15 | n"),
    ]
    fig, ax = plt.subplots()
    for mask, color, label in classes:
        ax.scatter(n[mask], y[mask], s=1.6, c=color, alpha=0.45, lw=0, label=label)
    ax.set_xlim(0, COMET_LIMIT)
    ax.set_ylim(0, y.max() * 1.06)
    ax.xaxis.set_major_formatter(thousands)
    ax.yaxis.set_major_formatter(thousands)
    ax.set_xlabel("even n")
    ax.set_ylabel("g(n) — ordered prime pairs with p + q = n")
    leg = ax.legend(loc="upper left", markerscale=8, title="divisibility of n")
    leg.get_title().set_color(INK_2)
    finish(fig, "The Goldbach comet",
           "Partition counts fan into bands — the bands follow the small prime factors of n, not chance")
    fig.savefig(FIG / "01_comet.png")
    plt.close(fig)


def collapse_figure(g, s):
    n = np.arange(6, COMET_LIMIT + 1, 2)
    y = g[n] / s[n]
    fig, ax = plt.subplots()
    ax.scatter(n, y, s=1.6, c=BLUE, alpha=0.35, lw=0)
    smooth_n = np.linspace(500, COMET_LIMIT, 400)
    smooth = 2.0 * TWIN_PRIME_C2 * li2(smooth_n)
    ax.plot(smooth_n, smooth, c=ORANGE, lw=2)
    ax.annotate("2 C₂ · Li₂(n)   (Hardy–Littlewood, 1923)",
                xy=(smooth_n[290], smooth[290]), xytext=(96_000, 380),
                color=INK_2, fontsize=9,
                arrowprops=dict(arrowstyle="-", color=BASELINE, lw=0.8))
    ax.set_xlim(0, COMET_LIMIT)
    ax.set_ylim(0, y.max() * 1.06)
    ax.xaxis.set_major_formatter(thousands)
    ax.yaxis.set_major_formatter(thousands)
    ax.set_xlabel("even n")
    ax.set_ylabel("g(n) / S(n) — count ÷ singular series")
    finish(fig, "One division collapses the comet",
           "Dividing each count by its singular series S(n) merges every band into a single curve")
    fig.savefig(FIG / "02_comet_collapsed.png")
    plt.close(fig)


def ratio_figure(g, prediction, n_max):
    n = np.arange(1000, n_max + 1, 2)
    rho = g[n] / prediction[n // 2]
    # Percentile band + median in log-spaced bins.
    edges = np.geomspace(1000, n_max, 61)
    idx = np.searchsorted(edges, n, side="right") - 1
    centers, lo_band, hi_band, med = [], [], [], []
    for b in range(60):
        vals = rho[idx == b]
        if vals.size < 50:
            continue
        centers.append(np.sqrt(edges[b] * edges[b + 1]))
        lo_band.append(np.percentile(vals, 1))
        hi_band.append(np.percentile(vals, 99))
        med.append(np.median(vals))
    sample = slice(None, None, 200)
    fig, ax = plt.subplots()
    ax.scatter(n[sample], rho[sample], s=1.5, c=BLUE, alpha=0.22, lw=0)
    ax.fill_between(centers, lo_band, hi_band, color=BLUE, alpha=0.14, lw=0)
    ax.plot(centers, med, c=ORANGE, lw=2)
    ax.axhline(1.0, c=MUTED, lw=1, ls=(0, (4, 3)))
    ax.annotate("median of g(n)/g_HL(n)", xy=(centers[40], med[40]),
                xytext=(3e5, 1.13), color=INK_2, fontsize=9,
                arrowprops=dict(arrowstyle="-", color=BASELINE, lw=0.8))
    ax.annotate("1%–99% envelope", xy=(centers[30], hi_band[30]),
                xytext=(2e4, 1.24), color=INK_2, fontsize=9,
                arrowprops=dict(arrowstyle="-", color=BASELINE, lw=0.8))
    ax.set_xscale("log")
    ax.set_xlim(1000, n_max)
    ax.set_ylim(0.7, 1.35)
    ax.set_xlabel("even n  (log scale)")
    ax.set_ylabel("g(n) / g_HL(n)")
    finish(fig, "Calibrating the prediction",
           "Measured counts over the full Hardy–Littlewood prediction: the envelope tightens onto 1")
    fig.savefig(FIG / "03_hl_ratio.png")
    plt.close(fig)
    return rho, n


def spectrum_figure(n_max):
    alphas, mags, n_primes = prime_spectrum(n_max)
    grid = len(alphas)
    half = grid // 2 + 1
    a, m = alphas[1:half], mags[1:half] / n_primes
    fig, ax = plt.subplots()
    ax.scatter(a, m, s=0.8, c=BLUE, alpha=0.4, lw=0)
    # Predicted major-arc heights |mu(q)|/phi(q) at every reduced a/q, q <= 13,
    # with the measured value re-plotted as a visible dot inside each ring.
    from math import gcd
    mobius = {1: 1, 2: -1, 3: -1, 5: -1, 6: 1, 7: -1, 10: 1, 11: -1, 13: -1}
    phi = {1: 1, 2: 1, 3: 2, 5: 4, 6: 2, 7: 6, 10: 4, 11: 10, 13: 12}
    pa, ph, measured = [], [], []
    for q in mobius:
        for num in range(1, q // 2 + 1):
            if gcd(num, q) == 1:
                pa.append(num / q)
                ph.append(1.0 / phi[q])
                measured.append(mags[round(num / q * grid)] / n_primes)
    ax.scatter(pa, measured, s=12, c=BLUE, lw=0, zorder=3)
    ax.scatter(pa, ph, s=70, facecolors="none", edgecolors=ORANGE, lw=1.4, zorder=4)
    for frac, label, dy in ((Fraction(1, 2), "1/2", 1.4), (Fraction(1, 3), "1/3", 1.4),
                            (Fraction(1, 5), "1/5", 1.5), (Fraction(2, 5), "2/5", 1.5),
                            (Fraction(1, 6), "1/6", 0.6), (Fraction(1, 7), "1/7", 1.55),
                            (Fraction(1, 11), "1/11", 1.55), (Fraction(1, 13), "1/13", 0.52)):
        k = round(float(frac) * grid)
        ax.annotate(label, xy=(float(frac), mags[k] / n_primes * dy),
                    color=INK_2, fontsize=9, ha="center")
    ax.annotate("square-root cancellation noise floor", xy=(0.25, 4e-3),
                xytext=(0.25, 2.8e-2), color=INK_2, fontsize=9, ha="center",
                arrowprops=dict(arrowstyle="-", color=BASELINE, lw=0.8))
    ax.set_yscale("log")
    ax.set_xlim(0, 0.512)
    ax.set_ylim(3e-5, 3)
    ax.set_xlabel("α")
    ax.set_ylabel("|S(α)| / π(N)   (log scale)")
    from matplotlib.lines import Line2D
    handles = [
        Line2D([], [], ls="", marker="o", ms=5, color=BLUE,
               label="measured  |S(α)| / π(N)"),
        Line2D([], [], ls="", marker="o", ms=8, mfc="none", mec=ORANGE, mew=1.4,
               label="predicted  |μ(q)| / φ(q)  at  a/q"),
    ]
    ax.legend(handles=handles, loc="upper left")
    finish(fig, "The power spectrum of the primes",
           "S(α) = Σₚ e(pα) spikes exactly at rationals a/q — and only for squarefree q: no spike at 1/4 or 1/9")
    fig.savefig(FIG / "04_prime_spectrum.png")
    plt.close(fig)
    return {"grid": grid, "pi_N": n_primes,
            "peak_half": float(mags[grid // 2] / n_primes),
            "peak_third": float(mags[round(grid / 3)] / n_primes),
            "peak_quarter": float(mags[round(grid / 4)] / n_primes)}


def minimal_prime_figure(pmin, recs, n_max):
    n = np.arange(6, n_max + 1, 2)
    sample = slice(None, None, 97)
    fig, ax = plt.subplots()
    ax.scatter(n[sample], pmin[sample], s=2, c=BLUE, alpha=0.25, lw=0)
    rn, rp = zip(*recs)
    ax.scatter(rn, rp, s=42, c=ORANGE, zorder=3)
    guide_n = np.geomspace(20, n_max, 200)
    guide = np.log(guide_n) ** 2 * np.log(np.log(guide_n))
    ax.plot(guide_n, guide, c=MUTED, lw=1.2, ls=(0, (4, 3)))
    ax.annotate("(ln n)² ln ln n   heuristic scale", xy=(guide_n[120], guide[120]),
                xytext=(120, max(rp) * 0.72), color=INK_2, fontsize=9,
                arrowprops=dict(arrowstyle="-", color=BASELINE, lw=0.8))
    for x, y in recs[-2:]:
        ax.annotate(str(y), xy=(x, y), xytext=(x, y + max(rp) * 0.07),
                    color=INK_2, fontsize=9, ha="center")
    ax.set_xscale("log")
    ax.set_xlim(6, n_max * 2)
    ax.set_ylim(0, max(rp) * 1.45)
    ax.set_xlabel("even n  (log scale)")
    ax.set_ylabel("smallest prime p with n − p prime")
    from matplotlib.lines import Line2D
    handles = [
        Line2D([], [], ls="", marker="o", ms=5, color=BLUE,
               label="p_min(n)  (1-in-97 sample)"),
        Line2D([], [], ls="", marker="o", ms=8, color=ORANGE, label="record-setters"),
    ]
    ax.legend(handles=handles, loc="upper left")
    finish(fig, "How stubborn can an even number be?",
           "The hardest n still split with a tiny prime: records grow like a power of ln n, not of n")
    fig.savefig(FIG / "05_minimal_prime.png")
    plt.close(fig)


def floor_figure(g, prediction, n_max):
    ks = range(4, int(np.log2(n_max)))
    centers, floors, predicted_floors = [], [], []
    for k in ks:
        lo, hi = 2**k, min(2 ** (k + 1), n_max + 1)
        n = np.arange(lo + (lo % 2), hi, 2)
        centers.append(np.sqrt(lo * (hi - 1)))
        floors.append(int(g[n].min()))
        predicted_floors.append(float(prediction[n // 2].min()))
    fig, ax = plt.subplots()
    ax.plot(centers, predicted_floors, c=ORANGE, lw=2, ls=(0, (4, 3)),
            label="Hardy–Littlewood floor  min g_HL(n)")
    ax.plot(centers, floors, c=BLUE, lw=2, marker="o", ms=5,
            label="measured floor  min g(n)")
    ax.axhline(1, c=MUTED, lw=1)
    ax.annotate("g(n) = 0 would be a counterexample — the floor rides away from it",
                xy=(centers[2], 1), xytext=(centers[2], 2.6), color=INK_2, fontsize=9)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("dyadic window of even n  (log scale)")
    ax.set_ylabel("worst case in window  (log scale)")
    ax.legend(loc="upper left")
    finish(fig, "The floor keeps rising",
           "Worst-case count per dyadic window: the weakest n tracks the predicted minimum, far above zero")
    fig.savefig(FIG / "06_worst_case.png")
    plt.close(fig)
    return list(zip([float(c) for c in centers], floors, predicted_floors))


def main(n_max=20_000_000):
    FIG.mkdir(exist_ok=True)
    RES.mkdir(exist_ok=True)
    t0 = time.time()

    print(f"[1/6] FFT autocorrelation: all Goldbach counts to N = {n_max:,}")
    g, is_prime, residual = goldbach_counts(n_max)
    evens = np.arange(6, n_max + 1, 2)
    assert (g[evens] > 0).all(), "Goldbach counterexample?! check the sieve"
    r = unordered_counts(g, is_prime)
    print(f"      FFT rounding residual {residual:.2e}  (exactness margin vs 0.5)")

    print("[2/6] singular series + Hardy–Littlewood prediction")
    s = singular_series(n_max)
    from goldbach.hardy_littlewood import hl_prediction
    prediction = hl_prediction(evens.astype(np.float64), s[evens])
    pred_by_half = np.zeros(n_max // 2 + 1)
    pred_by_half[evens // 2] = prediction

    print("[3/6] comet + collapse figures")
    comet_figure(g)
    collapse_figure(g, s)

    print("[4/6] calibration ratio over the full range")
    rho, ns = ratio_figure(g, pred_by_half, n_max)
    tail = ns >= 1_000_000
    worst_idx = int(np.argmin(rho))

    print("[5/6] prime spectrum")
    spec = spectrum_figure(n_max)

    print("[6/6] minimal Goldbach primes + worst-case floor")
    pmin = minimal_goldbach_prime(n_max, is_prime)
    recs = records(pmin)
    minimal_prime_figure(pmin, recs, n_max)
    floors = floor_figure(g, pred_by_half, n_max)

    summary = {
        "N": n_max,
        "verified": f"every even 4 <= n <= {n_max:,} has a Goldbach partition "
                    "(4 = 2 + 2; all others with two odd primes)",
        "fft_rounding_residual": residual,
        "min_unordered_count_n_ge_1e6": int(r[evens][evens >= 1_000_000].min()),
        "ratio_worst": {"n": int(ns[worst_idx]), "rho": float(rho[worst_idx])},
        "ratio_tail_n_ge_1e6": {
            "median": float(np.median(rho[tail])),
            "p01": float(np.percentile(rho[tail], 1)),
            "p99": float(np.percentile(rho[tail], 99)),
        },
        "minimal_prime_records": recs,
        "spectrum": spec,
        "dyadic_floors": [
            {"window_center": c, "min_g": f, "min_hl": p} for c, f, p in floors
        ],
        "runtime_seconds": round(time.time() - t0, 1),
    }
    (RES / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2)[:1200])
    print(f"done in {summary['runtime_seconds']}s — figures/ and results/ written")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 20_000_000)
