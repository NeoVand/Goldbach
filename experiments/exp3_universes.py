"""Experiment 3: Goldbach in random universes.

200,000 random 'prime' universes are sampled from Cramer-type models and every
Goldbach failure is recorded. The empirical failure probability matches the
exact independence computation, decays like exp(-c n / log^2 n), and its tail
sum -- the expected number of failures ever, beyond any cutoff -- is
astronomically small. In a random universe, Goldbach is not merely true; it is
overwhelmingly forced. Separately, one universe with only the primes' local
(mod 30) structure reproduces the real comet's bands.

Outputs: figures/fig6_universes.png, fig7_model_comet.png,
         data/universe_summary.txt
"""

import pathlib
import sys

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import logsumexp

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))  # repo root
import style

from goldbach.models import (failure_frequencies, representation_counts,
                             sample_universe, theory_curve)
from goldbach.hardy_littlewood import hl_integral

style.apply()
ROOT = pathlib.Path(__file__).resolve().parents[1]
(ROOT / "figures").mkdir(exist_ok=True)
(ROOT / "data").mkdir(exist_ok=True)

LIMIT, TRIALS = 1500, 200_000
rng = np.random.default_rng(2026)
evens = np.arange(6, LIMIT + 1, 2)

print(f"sampling {TRIALS:,} universes up to n = {LIMIT} (naive + parity) ...")
freq_naive = failure_frequencies(LIMIT, TRIALS, rng, modulus=1)
freq_parity = failure_frequencies(LIMIT, TRIALS, rng, modulus=2)
th_naive = theory_curve(LIMIT, modulus=1)
th_parity = theory_curve(LIMIT, modulus=2)

# ---- Figure 6 ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(9, 4.8), layout="constrained")
for freq, th, color, label in [
        (freq_naive, th_naive, style.BLUE, "Cramér model (all k)"),
        (freq_parity, th_parity, style.ORANGE, "parity-aware model (odd k)")]:
    pos = freq > 0
    ax.scatter(evens[pos], freq[pos], s=9, c=color, alpha=0.6, lw=0,
               label=f"{label} — observed")
    ax.plot(evens, th, color=color, lw=1.4, alpha=0.9,
            label=f"{label} — exact model probability")
ax.axhline(1 / TRIALS, color=style.BASELINE, lw=0.8)
ax.annotate("resolution floor: a single failure among 200,000 universes",
            (1240, 1 / TRIALS), ha="right", va="bottom", fontsize=8.5,
            color=style.MUTED)
ax.set_yscale("log")
ax.set_ylim(1 / (TRIALS * 4), 1.2)
ax.set_xlim(0, 1250)
ax.set_xlabel("even number n")
ax.set_ylabel("P(n has no two-prime representation)")
ax.set_title(f"Goldbach failure probability in {TRIALS:,} random universes")
ax.legend(loc="upper right")
fig.savefig(ROOT / "figures/fig6_universes.png")
plt.close(fig)

# ---- Model comet: local structure only -------------------------------------
N_COMET = 1_000_000
print("sampling one mod-30 universe for the model comet ...")
uni = sample_universe(N_COMET, rng, modulus=30)
counts = representation_counts(uni)
ev = np.arange(6, N_COMET + 1, 2)
gm = counts[ev] // 2

classes = [
    ("3 ∤ n", (ev % 3 != 0), style.BLUE),
    ("3 | n, 5 ∤ n", (ev % 3 == 0) & (ev % 5 != 0), style.ORANGE),
    ("15 | n", (ev % 15 == 0), style.AQUA),
]
fig, ax = plt.subplots(figsize=(9, 5.0), layout="constrained")
for label, mask, color in classes:
    idx = np.nonzero(mask)[0]
    idx = rng.choice(idx, size=min(50_000, len(idx)), replace=False)
    ax.scatter(ev[idx], gm[idx], s=0.5, c=color, alpha=0.25, lw=0)
    ax.scatter([], [], s=18, c=color, label=label)
ax.set_xlabel("even number n")
ax.set_ylabel("two-'prime' representations (unordered)")
ax.set_title("A comet from coin flips: random universe with only the primes' "
             "mod-30 structure")
ax.legend(title="divisibility of n", loc="upper left")
ax.set_xlim(0, N_COMET)
ax.set_ylim(0, None)
fig.savefig(ROOT / "figures/fig7_model_comet.png")
plt.close(fig)

# ---- Tail sums: the Borel-Cantelli ledger ----------------------------------
lines = [
    f"trials: {TRIALS:,}; universes up to n = {LIMIT}",
    f"largest n that ever failed (naive model) : "
    f"{evens[np.nonzero(freq_naive)[0][-1]] if freq_naive.any() else '-'}",
    f"largest n that ever failed (parity model): "
    f"{evens[np.nonzero(freq_parity)[0][-1]] if freq_parity.any() else '-'}",
    "",
    "expected number of Goldbach failures at ANY n > x (parity model, "
    "exp(-J(n)/2) tail):",
]
for x in [10**3, 10**4, 10**6, 10**9]:
    # sum_{n > x, even} exp(-J(n)/2) in log space -- the tail underflows
    # float64 spectacularly, which is rather the point
    grid = np.exp(np.linspace(np.log(x), np.log(x) + 14, 4000))
    j = hl_integral(grid)
    dg = np.diff(grid)
    mids = 0.5 * (j[1:] + j[:-1])
    log_terms = np.log(0.5) - mids / 2 + np.log(dg)  # density of evens is 1/2
    log10_tail = logsumexp(log_terms) / np.log(10.0)
    lines.append(f"  x = 10^{int(np.log10(x)):<2}: 10^({log10_tail:,.1f})")
(ROOT / "data/universe_summary.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
