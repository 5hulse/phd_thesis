# utils.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Fri 19 May 2023 01:06:47 BST

from pathlib import Path
from typing import Iterable
import matplotlib as mpl
from matplotlib import colors
import nmrespy as ne
import numpy as np
from scipy.signal import argrelextrema


DPHIL_DIR = Path("~/Documents/DPhil").expanduser()
DATA_DIR = DPHIL_DIR / "data"
RESULT_DIR = DPHIL_DIR / "results"
THESIS_DIR = DPHIL_DIR / "thesis"
FIGURE_DIR = THESIS_DIR / "figures/output"


def get_pure_shift_labels(estimator: ne.Estimator2DJ, yshift: float = 0., thold=None):
    shifts = estimator.get_shifts(unit="ppm", meshgrid=False)[1]
    cupid_spectrum = estimator.cupid_spectrum().real
    argmaxima = list(argrelextrema(cupid_spectrum, np.greater)[0])
    label_ys = [cupid_spectrum[idx] for idx in argmaxima]
    if thold is not None:
        for i, y in enumerate(label_ys):
            if y < thold:
                label_ys.pop(i)
                argmaxima.pop(i)

    label_xs = [shifts[idx] for idx in argmaxima]
    label_ys = [label_y + yshift for label_y in label_ys]
    label_texts = [f"({chr(65 + i)})" for i, _ in enumerate(argmaxima)]

    return label_xs, label_ys, label_texts


def add_pure_shift_labels(
    axs: np.ndarray[mpl.axes.Axes],
    xs: Iterable[float],
    ys: Iterable[float],
    ss: Iterable[float],
    fs=None,
) -> None:
    active_ax = 0
    for i, (x, y, s) in enumerate(zip(xs, ys, ss)):
        while True:
            l, r = axs[active_ax].get_xlim()
            if l > x > r:
                axs[active_ax].text(
                    x, y, s, ha="center", zorder=1000, fontsize=fs,
                    bbox={"facecolor": "w", "pad": 0.5, "edgecolor": "none"},
                )
                break
            else:
                active_ax += 1


def fix_linewidths(axs, lw):
    for ax_row in axs:
        for ax in ax_row:
            for line in ax.get_lines():
                line.set_linewidth(lw)


def panel_labels(fig, x, ys, start=97):
    for i, y in enumerate(ys):
        fig.text(x, y, f"\\textbf{{{chr(start + i)}.}}")


def raise_axes(axs, scale):
    new_top = scale * axs[0][0].get_ylim()[1]
    for ax in axs[0]:
        ax.set_ylim(top=new_top)


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        f"trunc({cmap.name},{minval:.2f},{maxval:.2f})",
        cmap(np.linspace(minval, maxval, n)),
    )
    return new_cmap
