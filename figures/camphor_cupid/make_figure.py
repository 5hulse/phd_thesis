# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Thu 01 Jun 2023 14:22:58 BST

from pathlib import Path

import nmrespy as ne
from matplotlib.patches import ConnectionPatch
import matplotlib as mpl
import numpy as np
from utils import (
    add_pure_shift_labels,
    get_pure_shift_labels,
    fix_linewidths,
    panel_labels,
    raise_axes,
)

RESULT_DIR = Path("~/Documents/DPhil/results/cupid/camphor").expanduser()

# === CONFIGURATION ===
estimator_path = RESULT_DIR / "estimator_postedit"

residual_shift = 3e6
multiplet_shift = residual_shift + 1e6
onedim_shift = multiplet_shift + 1e6
ax_top = onedim_shift + 3e6

# =====================
estimator = ne.Estimator2DJ.from_pickle(estimator_path)
thold = 2. * (estimator.sw()[1] / estimator.default_pts[1])
colors = mpl.rcParams["axes.prop_cycle"].by_key()["color"]
colors = 5 * [colors[0]] + colors[1:] + colors
fig, axs = estimator.plot_result(
    axes_right=0.99,
    axes_bottom=0.09,
    axes_top=0.975,
    axes_left=0.05,
    multiplet_thold=thold,
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
        (0, (2.5,)),
        (2, (2.05,)),
        (3, (1.9, 1.85, 1.8)),
        (4, (1.675, 1.625)),
    ],
    ratio_1d_2d=(3., 1.),
    figsize=(6, 3),
)
fig.texts[0].set_fontsize(8)
axs[1, 0].set_yticks([-20, -10, 0, 10, 20])

fix_linewidths(axs, 0.8)
panel_labels(fig, 0.055, (0.94, 0.57, 0.42, 0.28))

_, shifts = estimator.get_shifts(unit="ppm", meshgrid=False)

# Remove strong coupling artefacts from CUPID spectrum
to_rm = [2, 5, 8, 9, 12, 13, 16, 18]
r5_params = estimator.get_params(indices=[5])
r5_params = np.delete(r5_params, to_rm, axis=0)
estimator._results[5].params = r5_params
r5_region = estimator.get_results(indices=[5])[0].get_region()
r5_slice = estimator.convert(r5_region, "hz->idx")[1]
cupid_spectrum = estimator.cupid_spectrum().real
shifts = estimator.get_shifts(unit="ppm", meshgrid=False)[1]
shifts_r5 = shifts[r5_slice[0]:r5_slice[1]]
spec_without_sc = cupid_spectrum[r5_slice[0]:r5_slice[1]]
prev_spec = axs[0][5].get_lines()[-2]
vshift = prev_spec.get_ydata()[0] - spec_without_sc[0]
prev_spec.set_color("#b0b0b0")
specline = axs[0][5].plot(shifts_r5, spec_without_sc + vshift, color="k", lw=0.8)

raise_axes(axs, 1.07)
n_ax = axs[0].size
for i, (ax0, ax1) in enumerate(zip(axs[0], axs[1])):
    if i in [0, n_ax - 1]:
        lines = ax1.get_lines()[:-1]
    else:
        lines = ax1.get_lines()[:-2]
    for line in lines:
        line.remove()
        if i == 0:
            continue
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
        ax0.axvline(x, color=line.get_color(), lw=0.5, zorder=-1)

# Label pure shift peaks
xs, ys, ss = get_pure_shift_labels(estimator, yshift=1e6)
xs = [xs[2]] + xs[5:]
ys = [ys[2]] + ys[5:]
ss = ["DMSO"] + ss[:-5]
add_pure_shift_labels(axs[0], xs, ys, ss, fs=7)


fig.savefig("figures/camphor_cupid/camphor_cupid.pdf")
