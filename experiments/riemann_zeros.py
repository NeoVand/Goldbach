"""Experiment: are the Riemann zeros audible in Goldbach partition data?

Hypothesis (from Fujii's 1991 explicit formula): the cumulative weighted
Goldbach count deviates from X^2/2 by a sum of log-periodic oscillations
whose angular frequencies are the imaginary parts of the zeta zeros and
whose amplitudes are 4/|rho(rho+1)|.

Design: compute the deviation signal from the real primes; compute the same
signal from a Cramer random "fake prime" set with the same density (the
control); periodogram both; compare peak locations and heights against the
known zeros.  If the zeros are real structure and not artifact, the prime
signal shows them and the control shows nothing.

Usage: python experiments/riemann_zeros.py [N]     (default 20_000_000)
"""

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import numpy as np
import matplotlib.pyplot as plt

from experiments.style import BASELINE, BLUE, INK_2, MUTED, ORANGE, finish
from goldbach.weighted import (
    ZETA_ZEROS, cumulative_deviation, log_periodogram, mangoldt,
    predicted_amplitudes, weighted_counts,
)

X_MIN = 10_000  # start of the analysis window in X


def cramer_control(n_max: int, seed: int = 42) -> np.ndarray:
    """A Cramer random model of Lambda: weight log n with probability 1/log n.

    Matches the von Mangoldt weighting in mean exactly (E[w(n)] = 1), so the
    control's cumulative sum has the same X^2/2 main term — but its 'primes'
    carry no zeta zeros.  Whatever survives this pipeline is generic noise.
    """
    n = np.arange(n_max + 1, dtype=np.float64)
    with np.errstate(divide="ignore"):
        log_n = np.where(n >= 2, np.log(np.maximum(n, 2)), 1.0)
    rng = np.random.default_rng(seed)
    hits = rng.random(n_max + 1) < 1.0 / log_n
    hits[:2] = False
    return np.where(hits, log_n, 0.0)


def main(n_max: int = 20_000_000):
    fig_dir, res_dir = ROOT / "figures", ROOT / "results"
    fig_dir.mkdir(exist_ok=True)
    res_dir.mkdir(exist_ok=True)

    print(f"[1/4] weighted Goldbach counts for the real primes, N = {n_max:,}")
    x, d_primes = cumulative_deviation(weighted_counts(mangoldt(n_max)))

    print("[2/4] the same signal for a Cramer random control set")
    _, d_control = cumulative_deviation(weighted_counts(cramer_control(n_max)))

    print("[3/4] periodograms in u = log X")
    freqs, amp_primes = log_periodogram(x, d_primes, X_MIN, n_max)
    _, amp_control = log_periodogram(x, d_control, X_MIN, n_max)
    predicted = predicted_amplitudes()

    # Peak match: nearest local maximum of the prime periodogram to each zero.
    resolution = freqs[1] - freqs[0]
    matches = []
    for gamma, a_pred in zip(ZETA_ZEROS, predicted):
        band = (freqs > gamma - 1.5) & (freqs < gamma + 1.5)
        k = np.nonzero(band)[0][np.argmax(amp_primes[band])]
        matches.append({
            "zero": float(gamma),
            "peak_at": float(freqs[k]),
            "offset": float(freqs[k] - gamma),
            "measured_amplitude": float(amp_primes[k]),
            "predicted_amplitude": float(a_pred),
            "control_amplitude_there": float(amp_control[k]),
        })

    # Significance: template statistic T(delta) = sum over first 5 zeros of
    # amplitude at (gamma_k + delta) / predicted_k.  If the peaks are really
    # the zeros, T is maximal at delta = 0 against a fine grid of shifts.
    def template_stat(delta):
        return float(sum(
            amp_primes[np.argmin(np.abs(freqs - (g + delta)))] / p
            for g, p in zip(ZETA_ZEROS[:5], predicted[:5])))

    shifts = np.arange(-6.0, 6.0001, 0.05)
    null_shifts = shifts[np.abs(shifts) > 2 * resolution]
    t0 = template_stat(0.0)
    t_null = np.array([template_stat(d) for d in null_shifts])
    p_value = float((t_null >= t0).mean())

    # Signal-to-floor: the primes' own periodogram noise floor away from zeros.
    in_band = (freqs > 10) & (freqs < 60)
    off_zero = in_band & (np.min(np.abs(freqs[:, None] - ZETA_ZEROS[None, :]),
                                 axis=1) > 1.5)
    noise_floor = float(np.median(amp_primes[off_zero]))

    print("[4/4] figure")
    keep = freqs <= 60.0
    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(9.6, 7.6), height_ratios=[1, 1.6])
    fig.subplots_adjust(hspace=0.42)

    show = (x >= X_MIN)[:: 501]
    ax0.plot(x[::501][show], d_primes[::501][show], c=BLUE, lw=0.9)
    ax0.set_xscale("log")
    ax0.set_xlabel("X  (log scale)")
    ax0.set_ylabel("D(X)")
    ax0.set_title("the deviation signal  D(X) = (Σ G(n) − X²/2) / X^3/2",
                  loc="left", fontsize=10, color=INK_2)

    for gamma in ZETA_ZEROS:
        ax1.axvline(gamma, c=BASELINE, lw=0.9, ls=(0, (4, 3)), zorder=1)
    ax1.plot(freqs[keep], amp_control[keep], c=ORANGE, lw=1.0, alpha=0.55,
             zorder=2, label="Cramér random control (no zeros hidden inside)")
    ax1.plot(freqs[keep], amp_primes[keep], c=BLUE, lw=1.7, zorder=3,
             label="real primes")
    ax1.scatter(ZETA_ZEROS, predicted, s=70, facecolors="none",
                edgecolors=ORANGE, lw=1.4, zorder=4,
                label="Fujii prediction  4/|ρ(ρ+1)|  at each zero")
    top = predicted[0] * 1.75
    label_y = {0: 1.28, 1: 1.5, 2: 1.5, 3: 1.75, 4: 1.32}
    for i, gamma in enumerate(ZETA_ZEROS[:5]):
        ax1.annotate(f"γ{'₁₂₃₄₅'[i]} = {gamma:.2f}", xy=(gamma, predicted[i]),
                     xytext=(gamma, predicted[i] * label_y[i]), color=INK_2,
                     fontsize=8.5, ha="center")
    ax1.annotate("the structureless control is ~10× louder — and tone-deaf:\n"
                 "its broadband noise clips the top of this panel",
                 xy=(9.0, top * 0.97), color=INK_2, fontsize=8.5, ha="center",
                 va="top")
    ax1.set_xlim(0, 60)
    ax1.set_ylim(0, top)
    ax1.set_xlabel("angular frequency γ in u = log X")
    ax1.set_ylabel("tone amplitude")
    ax1.set_title("periodogram of D — dashed lines mark the known zeta zeros",
                  loc="left", fontsize=10, color=INK_2)
    ax1.legend(loc="upper right")

    finish(fig, "The Riemann zeros, heard through Goldbach",
           "The cumulative Goldbach error oscillates in log X at exactly the zeta zeros' frequencies — "
           "the control does not")
    fig.savefig(fig_dir / "07_riemann_zeros.png")
    plt.close(fig)

    out = {
        "N": n_max,
        "x_min": X_MIN,
        "frequency_resolution": float(resolution),
        "matches": matches,
        "mean_abs_offset_first5": float(np.mean([abs(m["offset"]) for m in matches[:5]])),
        "template_stat_at_zero_shift": t0,
        "template_shift_p_value": p_value,
        "primes_noise_floor_10_60": noise_floor,
        "gamma1_snr_vs_own_floor": float(matches[0]["measured_amplitude"] / noise_floor),
        "control_band_amplitude_10_60_median": float(np.median(amp_control[off_zero])),
    }
    (res_dir / "riemann_zeros.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 20_000_000)
