"""Experiment 4: first-witness records -- the hardest even numbers.

For every even n <= 10^7, find the least prime p with n - p prime. The
record-setters (numbers whose smallest witness is a new maximum) are the
conjecture's most reluctant cases; their witnesses grow only like log^2 n,
which is why the conjecture "feels" so safe: even the worst case gets easier,
relatively, as n grows.

Outputs: figures/fig8_witnesses.png, data/witness_records.txt
"""

import pathlib
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))  # repo root
import style

from goldbach.witnesses import first_witnesses, witness_records

style.apply()
ROOT = pathlib.Path(__file__).resolve().parents[1]
(ROOT / "figures").mkdir(exist_ok=True)
(ROOT / "data").mkdir(exist_ok=True)

LIMIT = 10_000_000
print(f"finding the least Goldbach witness for every even n <= {LIMIT:,} ...")
evens, w = first_witnesses(LIMIT)
records = witness_records(evens, w)

# ---- Figure 8 ---------------------------------------------------------------
rng = np.random.default_rng(3)
idx = rng.choice(len(evens), size=120_000, replace=False)
rec_n = np.array([n for n, _ in records])
rec_p = np.array([p for _, p in records])

fig, ax = plt.subplots(figsize=(9, 5.0), layout="constrained")
ax.scatter(evens[idx], w[idx], s=0.5, c=style.BLUE, alpha=0.15, lw=0)
ax.scatter([], [], s=18, c=style.BLUE, label="least witness p(n)")
ax.plot(rec_n, rec_p, drawstyle="steps-post", color=style.ORANGE, lw=1.6,
        label="running record")
ax.scatter(rec_n, rec_p, s=22, c=style.ORANGE, zorder=3,
           edgecolor=style.SURFACE, linewidth=1.2)
xs = np.exp(np.linspace(np.log(20), np.log(LIMIT), 300))
ax.plot(xs, 1.2 * np.log(xs) ** 2 * np.log(np.log(xs)), color=style.INK2,
        lw=1.2, label="1.2 · log²n · log log n  (heuristic envelope)")
offsets = {records[-1][0]: (6, 6, "left"), records[-2][0]: (-8, -4, "right"),
           records[-3][0]: (-6, 8, "right")}
for n, p in records[-3:]:
    dx, dy, ha = offsets[n]
    ax.annotate(f"n = {n:,}", (n, p), textcoords="offset points",
                xytext=(dx, dy), ha=ha, fontsize=8.5, color=style.INK2)
ax.set_xscale("log")
ax.set_xlabel("even number n (log scale)")
ax.set_ylabel("least prime witness p(n)")
ax.set_title("The hardest cases: record-setting least Goldbach witnesses, n ≤ 10⁷")
ax.legend(loc="upper left")
fig.savefig(ROOT / "figures/fig8_witnesses.png")
plt.close(fig)

# ---- Records table ----------------------------------------------------------
lines = ["record-setting even n and least witness p (n = p + prime):", ""]
lines += [f"  n = {n:>10,}   p = {p:>4}" for n, p in records]
lines.append(f"\nmean witness over all even n <= {LIMIT:,}: {w.mean():.2f}")
lines.append(f"99.99th percentile witness: {np.percentile(w, 99.99):.0f}")
(ROOT / "data/witness_records.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
