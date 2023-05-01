# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Mon 01 May 2023 20:08:31 BST

from pathlib import Path

import nmrespy as ne
from matplotlib.patches import ConnectionPatch
import matplotlib as mpl
import numpy as np
from utils import (
    THESIS_DIR,
    add_pure_shift_labels,
    get_pure_shift_labels,
    fix_linewidths,
    panel_labels,
    raise_axes,
)

RESULT_DIR = Path("/tmp/results/cupid").resolve()

# === CONFIGURATION ===
estimator_path = RESULT_DIR / "quinine/estimator_postedit"
save_path = "figures/quinine_cupid/quinine_cupid.pdf"

residual_shift = 3e6
multiplet_shift = residual_shift + 1e6
onedim_shift = multiplet_shift + 1e6
ax_top = onedim_shift + 3e6

# =====================


estimator = ne.Estimator2DJ.from_pickle(estimator_path)

colors = mpl.rcParams["axes.prop_cycle"].by_key()["color"]
fig, axs = estimator.plot_result(
    axes_right=0.98,
    axes_bottom=0.09,
    axes_top=0.975,
    axes_left=0.055,
    multiplet_thold=1.2,
    region_unit="ppm",
    contour_base=1.8e4,
    contour_nlevels=10,
    contour_factor=1.6,
    contour_color="k",
    contour_lw=0.1,
    multiplet_colors=colors,
    marker_size=2.,
    multiplet_show_45=False,
    jres_sinebell=True,
    xaxis_label_height=0.015,
    xaxis_ticks=[
        (0, (5.8, 5.7, 5.6)),
        (1, (4.95,)),
        (2, (3.7,)),
        (3, (3.1,)),
        (4, (2.75, 2.65)),
        (5, (1.95, 1.85, 1.75)),
        (6, (1.6, 1.5, 1.4)),
    ],
    ratio_1d_2d=(3., 1.),
    figsize=(6, 3),
)
fig.texts[0].set_fontsize(8)
axs[1, 0].set_yticks([-20, -10, 0, 10, 20])

fix_linewidths(axs, mpl.rcParams["lines.linewidth"])
panel_labels(fig, 0.06, (0.94, 0.5, 0.4, 0.28))

_, shifts = estimator.get_shifts(unit="ppm", meshgrid=False)

# Remove H2O from CUPID spectrum
r1_params = estimator.get_params(indices=[1])
r1_params = np.delete(r1_params, (4), axis=0)
estimator._results[1].params = r1_params
r1_region = estimator.get_results(indices=[1])[0].get_region()
r1_slice = estimator.convert(r1_region, "hz->idx")[1]
cupid_spectrum = estimator.cupid_spectrum().real
shifts = estimator.get_shifts(unit="ppm", meshgrid=False)[1]
shifts_r1 = shifts[r1_slice[0]:r1_slice[1]]
spec_without_h2o = cupid_spectrum[r1_slice[0]:r1_slice[1]]
prev_spec = axs[0][1].get_lines()[-3]
vshift = prev_spec.get_ydata()[0] - spec_without_h2o[0]
prev_spec.remove()
specline = axs[0][1].plot(shifts_r1, spec_without_h2o + vshift, color="k")

raise_axes(axs, 1.07)
n_ax = axs[0].size
for i, (ax0, ax1) in enumerate(zip(axs[0], axs[1])):
    if i in [0, n_ax - 1]:
        lines = ax1.get_lines()[:-1]
    else:
        lines = ax1.get_lines()[:-2]
    for line in lines:
        line.remove()
        x = line.get_xdata()[0]
        con = ConnectionPatch(
            xyA=(x, ax1.get_ylim()[0]),
            xyB=(x, ax0.get_ylim()[0]),
            coordsA="data",
            coordsB="data",
            axesA=ax1,
            axesB=ax0,
            color=line.get_color(),
            lw=0.5,
            zorder=-1,
        )
        ax1.add_patch(con)
        ax0.axvline(x, color=line.get_color(), lw=0.5, zorder=500)

# Label pure shift peaks
xs, ys, ss = get_pure_shift_labels(estimator, yshift=2.1e6)
xs[3] -= 0.008
xs[6] -= 0.018
xs[7] += 0.015
xs[10] -= 0.004
ss[3] = "(D)*"
add_pure_shift_labels(axs[0], xs, ys, ss, fs=8)


fig.text(
    0.32,
    0.41,
    "H\\textsubscript{2}\\hspace{-0.7pt}O",
    color=colors[4],
    fontsize=6,
)

fig.savefig(save_path)
