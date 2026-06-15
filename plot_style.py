"""
plot_style.py
-------------
Shared aesthetics for the DIF sample-size analysis notebook.
Import at the top of the notebook:  import plot_style as ps

Usage
-----
ps.setup()           # apply rcParams globally
ps.band(ax, x, y, sd, color, label="…")   # ribbon + line
ps.add_thresholds(ax, n_break_p, n_trust_contrast)
ps.logx(ax)          # conditional log-x (reads USE_LOG_X from caller scope indirectly)
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize

# ---------------------------------------------------------------------------
# Palette
# ---------------------------------------------------------------------------
C = {
    "pval":         "#44ff0c",   # p-value method
    "contrast":     "#4a4db5",   # contrast method
    "intersection": "#16a085",   # intersection rule
    "dif":          "#7d3c98",   # true DIF items
    "nodif":        "#7f5539",   # non-DIF items
    "breakdown":    "#c0392b",   # p-value breakdown vline
    "convergence":  "#4a4db5",   # contrast convergence vline
    "neutral":      "#555555",   # reference lines / thresholds
    "star1":        "#e67e22",   # * significance level
    "star2":        "#27ae60",   # ** significance level
    "star3":        "#c0392b",   # *** significance level
}

# ---------------------------------------------------------------------------
# Global setup
# ---------------------------------------------------------------------------
def setup(dpi: int = 150, base_font: int = 11) -> None:
    """Apply rcParams. Call once at notebook top."""
    plt.style.use("seaborn-v0_8-paper")
    mpl.rcParams.update({
        "figure.dpi":           dpi,
        "font.size":            base_font,
        "axes.titlesize":       base_font + 1,
        "axes.labelsize":       base_font,
        "legend.fontsize":      base_font - 1,
        "axes.spines.top":      False,
        "axes.spines.right":    False,
        "axes.grid":            True,
        "grid.alpha":           0.28,
        "grid.linestyle":       "--",
        "legend.frameon":       True,
        "legend.framealpha":    0.88,
        "legend.edgecolor":     "#cccccc",
        "lines.markersize":     4,
    })


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def band(
    ax,
    x,
    y,
    sd,
    color: str,
    alpha: float = 0.15,
    lw: float = 1.8,
    ms: float = 3.5,
    **line_kw,
) -> None:
    """Ribbon (±1 SD) + line plot."""
    ax.fill_between(x, np.maximum(y - sd, 0), y + sd, alpha=alpha, color=color)
    ax.plot(x, y, marker="o", markersize=ms, linewidth=lw, color=color, **line_kw)


def add_thresholds(ax, n_break_p, n_trust_contrast, ymax=None) -> None:
    """Vertical dashed reference lines for the two key n thresholds."""
    if n_break_p is not None:
        ax.axvline(
            n_break_p, color=C["breakdown"], linestyle=":",
            linewidth=1.6, alpha=0.85,
            label=f"p-value breakdown  (n = {n_break_p:,})",
        )
    if n_trust_contrast is not None:
        ax.axvline(
            n_trust_contrast, color=C["convergence"], linestyle=":",
            linewidth=1.6, alpha=0.85,
            label=f"contrast trusted  (n = {n_trust_contrast:,})",
        )


def logx(ax, use_log: bool = True) -> None:
    if use_log:
        ax.set_xscale("log")


def n_colormap(n_values, cmap_name: str = "plasma"):
    """Return (norm, cmap, sm) for colouring points by sample size."""
    norm = Normalize(vmin=min(n_values), vmax=max(n_values))
    cmap = cm.get_cmap(cmap_name)
    sm   = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    return norm, cmap, sm


def save(fig, path, **kw):
    fig.savefig(path, dpi=300, bbox_inches="tight", **kw)
