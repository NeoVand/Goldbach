"""Experiment 1: the Goldbach comet and its Hardy-Littlewood deconvolution.

Computes exact partition counts for every even n <= 10^7 (one FFT), then
divides out the predicted count. The famous banded comet collapses to a single
tight ribbon around 1 -- the entire visible structure of the Goldbach counts
is the singular series; what remains is featureless noise.

Outputs: figures/fig1_comet.png, figures/fig2_deconvolved.png,
         data/comet_summary.txt
"""

import pathlib
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))  # repo root
import style

from goldbach.partitions import goldbach_counts
from goldbach.hardy_littlewood import hl_prediction

style.apply()
ROOT = pathlib.Path(__file__).resolve().parents[1]
(ROOT / "figures").mkdir(exist_ok=True)
(ROOT / "data").mkdir(exist_ok=True)

LIMIT = 10_000_000

print(f"counting Goldbach partitions for all even n <= {LIMIT:,} ...")
evens, g = goldbach_counts(LIMIT)
ev, pred = hl_prediction(LIMIT)
r2 = 2 * g[1:].astype(np.float64)  # ordered counts; diagonal term negligible
ratio = r2 / pred

classes = [
    ("3 ∤ n", (ev % 3 != 0), style.BLUE),
    ("3 | n, 5 ∤ n", (ev % 3 == 0) & (ev % 5 != 0), style.ORANGE),
    ("15 | n", (ev % 15 == 0), style.AQUA),
]

# ---- Figure 1: the comet ----------------------------------------------------
rng = np.random.default_rng(1)
fig, ax = plt.subplots(figsize=(9, 5.4), layout="constrained")
for label, mask, color in classes:
    idx = np.nonzero(mask)[0]
    idx = rng.choice(idx, size=min(60_000, len(idx)), replace=False)
    ax.scatter(ev[idx], g[1:][idx], s=0.5, c=color, alpha=0.25, lw=0)
    ax.scatter([], [], s=18, c=color, label=label)  # legend proxies
ax.set_xlabel("even number n")
ax.set_ylabel("Goldbach partitions g(n)")
ax.set_title("The Goldbach comet: g(n) for every even n ≤ 10⁷")
ax.legend(title="divisibility of n", loc="upper left")
ax.set_xlim(0, LIMIT)
ax.set_ylim(0, None)
fig.savefig(ROOT / "figures/fig1_comet.png")
plt.close(fig)

# ---- Figure 2: the deconvolved comet ---------------------------------------
fig, (ax, axh) = plt.subplots(
    1, 2, figsize=(9, 4.2), layout="constrained",
    gridspec_kw={"width_ratios": [3.2, 1]}, sharey=True)
for label, mask, color in classes:
    idx = np.nonzero(mask)[0]
    idx = rng.choice(idx, size=min(60_000, len(idx)), replace=False)
    ax.scatter(ev[idx], ratio[idx], s=0.5, c=color, alpha=0.25, lw=0)
ax.axhline(1.0, color=style.INK2, lw=1.0)
ax.set_xlabel("even number n")
ax.set_ylabel("actual / predicted representations")
ax.set_title("Deconvolved: counts ÷ Hardy–Littlewood prediction")
ax.set_xlim(0, LIMIT)
ax.set_ylim(0.8, 1.2)

sample = ratio[ev > LIMIT // 10]
axh.hist(sample, bins=160, range=(0.8, 1.2), orientation="horizontal",
         color=style.BLUE, edgecolor="none")
axh.axhline(1.0, color=style.INK2, lw=1.0)
axh.set_xlabel("count of n")
axh.set_title("distribution (n > 10⁶)", fontsize=10)
fig.savefig(ROOT / "figures/fig2_deconvolved.png")
plt.close(fig)

# ---- Summary table ----------------------------------------------------------
lines = ["decade            mean(ratio)  sd(ratio)  min g(n)   n at min",
         "-" * 62]
for lo in [10**k for k in range(1, 7)]:
    hi = lo * 10
    m = (ev >= lo) & (ev < hi)
    if not m.any():
        continue
    gm = g[1:][m]
    lines.append(f"[{lo:>9,}, {hi:>10,})  {ratio[m].mean():9.4f}  "
                 f"{ratio[m].std():9.4f}  {gm.min():7d}  {ev[m][gm.argmin()]:>9,}")
mn = g.min()
lines.append(f"\nglobal minimum of g(n): {mn} (n = {evens[g.argmin()]});"
             f" zero partitions found: {(g == 0).sum()}")
(ROOT / "data/comet_summary.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
