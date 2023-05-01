# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Mon 01 May 2023 19:53:16 BST

from pathlib import Path
import pickle

import matplotlib as mpl
import nmrespy as ne
import numpy as np
from scipy.signal import argrelextrema

from utils import (
    # RESULT_DIR,
    fix_linewidths,
    panel_labels,
    raise_axes,
)

RESULT_DIR = Path("/tmp/results").resolve()

result_dir = RESULT_DIR / "cupid/sucrose"
estimator_path = result_dir / "estimator"
estimator_dec_path = result_dir / "estimator_dec"
save_path = "figures/sucrose_cupid/sucrose_cupid.pdf"

estimator = ne.Estimator2DJ.from_pickle(estimator_path)
reg3 = estimator._results[3].region
reg4 = estimator._results[4].region
estimator._results[3].region = (None, (reg3[1][0], reg4[1][0]))

colors = mpl.rcParams["axes.prop_cycle"].by_key()["color"]
colors.insert(0, colors.pop())

fig, axs = estimator.plot_result(
    indices=list(range(1, 7)),
    multiplet_vertical_shift=0.,
    axes_left=0.05,
    axes_right=0.99,
    axes_bottom=0.09,
    axes_top=0.98,
    xaxis_label_height=0.015,
    ratio_1d_2d=(3., 1.),
    contour_base=2.,
    contour_nlevels=12,
    contour_factor=1.5,
    contour_color="k",
    contour_lw=0.1,
    multiplet_colors=colors,
    marker_size=2.,
    figsize=(6, 3),
    region_unit="ppm",
    axes_region_separation=0.01,
    multiplet_show_45=False,
)

xlabel = fig.texts[0]
xlabel.set_fontsize(8)

for ax in axs[1]:
    ax.set_ylim(10, -10)
axs[1, 0].set_yticks([-10, -5, 0, 5, 10])

fix_linewidths(axs, mpl.rcParams["lines.linewidth"])
panel_labels(fig, 0.055, (0.945, 0.53, 0.4, 0.285))

estimator_dec = ne.Estimator2DJ.from_pickle(estimator_dec_path)
dec_spec = estimator_dec.spectrum_first_direct.real
shifts = estimator_dec.get_shifts(unit="ppm", meshgrid=False)[-1]

n_axs = len(axs[0])
for i, ax in enumerate(axs[0]):
    idx = -2 if i in (0, n_axs - 1) else -3
    ax.plot(
        shifts + 0.015,
        dec_spec + ax.get_lines()[idx].get_ydata()[0],
        color="#c0c0c0",
        zorder=1,
    )

with open(result_dir / "shifts.pkl", "rb") as fh:
    chem_shifts = pickle.load(fh)
label_texts = [
    f"({chr(i[0] + 65)})"
    for i in reversed(sorted(enumerate(chem_shifts), key=lambda x:x[1]))
]
n = label_texts.pop(11)[1:-1]
label_texts[10] = f"{label_texts[10][:-1]},{n})"
cupid_spectrum = estimator.cupid_spectrum()
argmaxima_cupid = list(argrelextrema(cupid_spectrum, np.greater)[0])
argmaxima_dec = list(argrelextrema(cupid_spectrum, np.greater)[0])
label_ys = [
    55 + max(cupid_spectrum[idx_cupid], dec_spec[idx_dec])
    for idx_cupid, idx_dec in zip(argmaxima_cupid, argmaxima_dec)
]
label_xs = [shifts[idx] for idx in argmaxima_cupid]

mx0, mn0 = axs[0, 0].get_xlim()
mx1, mn1 = axs[0, 1].get_xlim()

for i, (x, y, s) in enumerate(zip(label_xs, label_ys, label_texts)):
    # if i == 2:
    #     x += 0.012
    # elif i == 3:
    #     x -= 0.012
    if i == 7:
        x += 0.014
    elif i == 8:
        x -= 0.014
    if mn0 < x < mx0:
        axs[0, 0].text(
            x, y, s, ha="center",
            bbox={"facecolor": "w", "pad": 0.5, "edgecolor": "none"},
            zorder=1000, fontsize=8,
        )
    elif mn1 < x < mx1:
        axs[0, 1].text(
            x, y, s, ha="center",
            bbox={"facecolor": "w", "pad": 0.5, "edgecolor": "none"},
            zorder=1000, fontsize=8,
        )

raise_axes(axs, 1.07)

for ax_col in axs.T:
    lines = [line for line in ax_col[1].get_lines() if line.get_color() != "k"]
    for line in lines:
        x = line.get_xdata()[0]
        color = line.get_color()
        line.remove()
        for ax in ax_col:
            ax.axvline(x, color=color, alpha=0.4, ls="-", zorder=500, lw=0.5)

fig.savefig(save_path)
