"""Chart style: the skill's validated light palette, applied to matplotlib.

Categorical slots 1–3 (blue, orange, aqua) validate all-pairs for scatter use;
the aqua slot is sub-3:1 on this surface, so every chart using it ships a
visible legend (the relief rule). Ink roles are used for all text — series
color never carries text.
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
BASELINE = "#c3c2b7"

BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"

plt.rcParams.update(
    {
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "figure.figsize": (9.6, 5.6),
        "figure.dpi": 110,
        "savefig.dpi": 200,
        "savefig.bbox": "tight",
        "font.family": "sans-serif",
        "axes.edgecolor": BASELINE,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "axes.axisbelow": True,
        "axes.labelcolor": INK_2,
        "axes.labelsize": 10,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.frameon": False,
        "legend.fontsize": 9,
        "legend.labelcolor": INK_2,
        "text.color": INK,
    }
)

thousands = FuncFormatter(lambda v, _: f"{v:,.0f}".replace(",", " "))


def finish(fig, title: str, subtitle: str):
    """Left-aligned title in primary ink, one-line subtitle in secondary ink."""
    fig.subplots_adjust(top=0.84)
    fig.text(0.07, 0.955, title, fontsize=13, fontweight="bold", color=INK)
    fig.text(0.07, 0.895, subtitle, fontsize=10, color=INK_2)
