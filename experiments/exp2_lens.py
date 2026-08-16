"""Experiment 2: the Farey lens -- the circle method as an instrument.

The primes up to N = 2^21, weighted by log p, are Fourier-transformed once.
Their spectrum |S(alpha)| is a spiked interference pattern concentrated at
rationals a/q. Filtering the spectrum of S^2 down to neighborhoods of Farey
fractions of level Q and inverse-transforming reconstructs the weighted
Goldbach counts R(n) = sum_{p+q=n} log p log q for every n at once, from
level-Q arcs only. As Q grows, the reconstruction sharpens into the banded
comet; the discarded remainder is the minor-arc term -- the exact quantity
whose size no one can prove. Here we measure it.

Outputs: figures/fig3_spectrum.png, fig4_farey_levels.png, fig5_minor.png,
         data/lens_summary.txt
"""

import pathlib
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))  # repo root
import style

from goldbach.lens import FareyLens
from goldbach.hardy_littlewood import singular_series, TWIN_PRIME_C2

style.apply()
ROOT = pathlib.Path(__file__).resolve().parents[1]
(ROOT / "figures").mkdir(exist_ok=True)
(ROOT / "data").mkdir(exist_ok=True)

N = 1 << 21
BINS = 16  # arc half-width in FFT bins; bin width = 1/2^23 of the circle
print(f"building Farey lens on primes <= {N:,} (FFT size 2^23) ...")
lens = FareyLens(N)

# ---- Figure 3: the interference pattern of the primes ----------------------
alpha, mag = lens.magnitude(stride=1)
# peak-preserving downsample: max over blocks, so spikes survive thinning
block = 64
usable = len(mag) - len(mag) % block
mag_ds = mag[:usable].reshape(-1, block).max(axis=1)
alpha_ds = alpha[:usable:block]

fig, ax = plt.subplots(figsize=(9, 4.6), layout="constrained")
ax.plot(alpha_ds, mag_ds, color=style.BLUE, lw=0.5)
ax.set_yscale("log")
ax.set_xlabel("frequency α (fraction of the circle)")
ax.set_ylabel("|S(α)|   (log scale)")
ax.set_title("The interference pattern of the primes: |Σ log p · e(2πi p α)|, p ≤ 2²¹")
for a, q, ha in [(0, 1, "left"), (1, 2, "right"), (1, 3, "center"),
                 (1, 6, "center"), (1, 5, "center"), (2, 5, "center"),
                 (3, 7, "center")]:
    x = a / q
    ax.annotate(f"{a}/{q}", (x, mag[round(x * lens.m)] * 1.35),
                ha=ha, fontsize=8.5, color=style.INK2)
ax.annotate("no spike at 1/4:  μ(4) = 0", (0.25, 1.2e4),
            ha="center", fontsize=8.5, color=style.INK2,
            arrowprops=dict(arrowstyle="-", color=style.MUTED, lw=0.8),
            xytext=(0.25, 3.2e5), annotation_clip=False)
ax.set_xlim(-0.01, 0.51)
ax.set_ylim(None, mag.max() * 3.5)
fig.savefig(ROOT / "figures/fig3_spectrum.png")
plt.close(fig)

# ---- Figure 4: the comet assembled one Farey level at a time ---------------
full = lens.full_field()
ns = np.arange(len(full))
# stop the window below n = N: R(n) has a corner there (pair truncation
# begins), and the hard spectral mask turns that corner into ringing
N_TOP = int(0.88 * N)
window = (ns >= N // 2) & (ns < N_TOP) & (ns % 2 == 0)
nw = ns[window]
true_norm = full[window] / nw

qs = [1, 3, 5, 15, 64]
fields = {q: lens.major_field(q, BINS)[window] / nw for q in qs}

panels = [(f"Farey level Q = {q}", fields[q]) for q in qs]
panels.append(("exact R(n)/n", true_norm))
fig, axes = plt.subplots(len(panels), 1, figsize=(8.6, 9.2), sharex=True,
                         layout="constrained")
hist_range = (0.4, 2.6)
for axp, (label, vals), shade in zip(
        axes, panels, [1, 2, 3, 4, 5, 6]):
    axp.hist(vals, bins=440, range=hist_range, color=style.SEQ_BLUE[shade],
             edgecolor="none", density=True)
    axp.set_yticks([])
    axp.grid(False)
    axp.spines["left"].set_visible(False)
    axp.text(0.995, 0.82, label, transform=axp.transAxes, ha="right",
             fontsize=9.5, color=style.INK2)
axes[-1].set_xlabel("R(n) / n   (weighted Goldbach count, normalized)")
fig.suptitle("The comet assembled one Farey level at a time\n"
             "distribution of major-arc reconstructions, n ∈ [2²⁰, 1.85·10⁶]",
             fontsize=11, fontweight="semibold")
fig.savefig(ROOT / "figures/fig4_farey_levels.png")
plt.close(fig)

# ---- Figure 5: the minor-arc remainder -------------------------------------
Q_STAR = 64
minor = full - lens.major_field(Q_STAR, BINS)
m_even = minor[window]
rng = np.random.default_rng(2)
idx = rng.choice(len(nw), size=60_000, replace=False)

fig, (ax, axh) = plt.subplots(
    1, 2, figsize=(9, 4.2), layout="constrained",
    gridspec_kw={"width_ratios": [3.2, 1]}, sharey=True)
scaled = m_even / np.sqrt(nw)
ylim = 1.3 * np.percentile(np.abs(scaled), 99.9)
ax.scatter(nw[idx], scaled[idx], s=0.5, c=style.BLUE, alpha=0.25, lw=0)
ax.axhline(0, color=style.INK2, lw=1.0)
ax.set_ylim(-ylim, ylim)
ax.set_xlabel("even number n")
ax.set_ylabel("minor-arc remainder  E(n) / √n")
ax.set_title(f"What the lens discards: E(n) = R(n) − R_major(n),  Q = {Q_STAR}")
axh.hist(scaled, bins=200, range=(-ylim, ylim), orientation="horizontal",
         color=style.BLUE, edgecolor="none")
axh.set_xlabel("count of n")
axh.set_title("distribution", fontsize=10)
fig.savefig(ROOT / "figures/fig5_minor.png")
plt.close(fig)

# ---- Summary ---------------------------------------------------------------
s_full = 2 * TWIN_PRIME_C2 * singular_series(N)  # index: n = 4 + 2i
sw = s_full[(nw - 4) // 2]
main_term = sw * nw
rel = np.abs(m_even) / main_term
rms_ratio = np.sqrt(np.mean((m_even / np.sqrt(nw)) ** 2))
lines = [
    f"primes up to N = 2^21 = {N:,}; FFT size 2^23; arc half-width "
    f"{BINS} bins = {BINS}/2^23 of the circle; window n in [2^20, {N_TOP:,})",
    "",
    f"major arcs Q = {Q_STAR}:",
    f"  RMS of E(n)/sqrt(n)          : {rms_ratio:.2f}",
    f"  max |E(n)| / main term       : {rel.max():.4f}"
    f"   (at n = {nw[rel.argmax()]:,})",
    f"  mean |E(n)| / main term      : {rel.mean():.4f}",
    f"  main term ~ S(n) n           : {main_term.mean():,.0f} (mean)",
    "",
    "reconstruction accuracy R_major vs exact R, correlation by level:",
]
for q in qs:
    c = np.corrcoef(fields[q], true_norm)[0, 1]
    lines.append(f"  Q = {q:>3}: corr = {c:.5f}")
(ROOT / "data/lens_summary.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
