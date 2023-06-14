# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Tue 13 Jun 2023 10:38:16 BST

from pathlib import Path
import matplotlib as mpl
from matplotlib.patches import ConnectionPatch
import numpy as np

import nmrespy as ne

from utils import (
    add_pure_shift_labels,
    get_pure_shift_labels,
    fix_linewidths,
    panel_labels,
    raise_axes,
)


# === CONFIGURATION ===
data_dir = Path("~/Documents/DPhil/data/estradiol").expanduser()
psyche_dir = data_dir / "6300/pdata/1"
estimator_path = Path(
    "~/Documents/DPhil/results/cupid/estradiol_low_snr/estimator_postedit"
).expanduser()

cupid_shift = 3.5e4
residual_shift = cupid_shift + 1.65e5
multiplet_shift = residual_shift + 1.7e4
onedim_shift = multiplet_shift + 7e3
ax_top = onedim_shift + 6e4

xticks = ([2.25], [2.05], [1.9, 1.8], [1.6], [1.4, 1.3, 1.2, 1.1])
thold_scale = 6.
ax_ylim = (-3.5e4, ax_top)
label_x = 0.017
label_ys = (0.95, 0.83, 0.76, 0.69, 0.4)
# =====================

estimator = ne.Estimator2DJ.from_pickle(estimator_path)
thold = thold_scale * estimator.default_multiplet_thold
fig, axs = estimator.plot_result(
    multiplet_thold=thold,
    region_unit="ppm",
    multiplet_colors=mpl.rcParams["axes.prop_cycle"].by_key()["color"],
    axes_bottom=0.085,
    axes_left=0.01,
    axes_right=0.99,
    axes_top=0.98,
    axes_region_separation=0.035,
    ratio_1d_2d=(100., 1.),
    contour_base=1.e8,
    contour_factor=1.2,
    contour_nlevels=10,
    multiplet_show_center_freq=True,
    figsize=(6, 3),
)
mp_centers = []
n_ax = axs[0].size
for i, ax in enumerate(axs[1]):
    n_lines = len(ax.get_lines())
    center_line_index = (n_lines - (0 if i in [0, n_ax - 1] else 1)) // 2
    mp_centers.append(
        [line.get_xdata()[0] for line in ax.get_lines()[:center_line_index]]
    )
    ax.remove()

fix_linewidths(axs, mpl.rcParams["lines.linewidth"])
axs = axs[0]

cupid_max = 0.
cupid_lines = []
psyche_estimator = ne.Estimator1D.new_bruker(psyche_dir)
psyche_spectra = []
psyche_shifts = []
marker = axs[0].get_lines()[-1].get_marker()
n_axs = axs.size
res = estimator.sw(unit="ppm")[-1] / estimator.default_pts[-1]
for i, (ax, mp_centre, xtick) in enumerate(zip(axs, mp_centers, xticks)):
    ax.spines["bottom"].set_visible(True)
    ax.set_xticks(xtick)
    if i in [0]:
        cupid_idx = -2
        xs = [1]
    elif i == n_axs - 1:
        cupid_idx = -2
        xs = [0]
    else:
        cupid_idx = -3
        xs = [0, 1]

    ys = len(xs) * [0]

    cupid_line = ax.get_lines()[cupid_idx]
    # Number was found by printing relavent shift out in source code for
    # Estimator2DJ.plot_result
    line_max = np.amax(cupid_line.get_ydata()) - 56291.79
    if line_max > cupid_max:
        cupid_max = line_max
    cupid_line.set_ydata(cupid_line.get_ydata() + cupid_shift)
    # hack: data I generate seems to be off by resolution / 2
    cupid_line.set_xdata(cupid_line.get_xdata() + res / 2)

    onedim_idx = cupid_idx - 1
    onedim_line = ax.get_lines()[onedim_idx]
    onedim_ydata = onedim_line.get_ydata()
    residual = onedim_ydata - 23540.499
    onedim_line.set_ydata(onedim_ydata + onedim_shift)
    # hack: data I generate seems to be off by resolution / 2
    onedim_line.set_xdata(onedim_line.get_xdata() + res / 2)

    multiplet_slice = slice(0, onedim_idx)
    shifts = onedim_line.get_xdata()
    slice_ = estimator.convert([None, [shifts[0], shifts[-1]]], "ppm->idx")[1]

    mp_lines = ax.get_lines()[multiplet_slice]
    if i == 0:
        mp_colors = iter([line.get_color() for line in mp_lines])
    for mp_line in mp_lines:
        mp_ydata = mp_line.get_ydata()
        residual -= mp_ydata[slice_[0] : slice_[1] + 1]
        mp_line.set_ydata(mp_ydata + multiplet_shift)
        # hack: data I generate seems to be off by resolution / 2
        mp_line.set_xdata(mp_line.get_xdata() + res / 2)

    for mpc in mp_centre:
        ax.axvline(mpc, color=next(mp_colors), zorder=0, lw=0.5)

    ax.plot(shifts, residual + residual_shift, color="k")

    left, right = ax.get_xlim()
    slc = psyche_estimator.convert([[left, right]], "ppm->idx")[0]
    psyche_shifts.append(psyche_estimator.get_shifts(unit="ppm")[0][slc[0] : slc[1] + 1])
    psyche_spectra.append(psyche_estimator.spectrum.real[slc[0] : slc[1] + 1])

    ax.plot(
        xs,
        ys,
        color="k",
        marker=marker,
        transform=ax.transAxes,
        clip_on=False,
        markersize=6,
        zorder=100,
    )
    ax.set_ylim(ax_ylim)

psyche_max = max([np.amax(psyche_spectrum) for psyche_spectrum in psyche_spectra])
psyche_scale = cupid_max / psyche_max

for ax, psyche_shift, psyche_spectrum in zip(axs, psyche_shifts, psyche_spectra):
    ax.plot(psyche_shift, psyche_scale * psyche_spectrum, color="k")

panel_labels(fig, label_x, label_ys)

xs, ys, ss = get_pure_shift_labels(estimator, yshift=9.9e4, thold=3e4)
xs[7] += 0.001
xs[8] -= 0.001
xs[10] -= 0.006
add_pure_shift_labels(axs, xs, ys, ss, fs=7)

fig.text(
    0.98,
    0.4,
    f"\\times {psyche_scale:.1f}",
    transform=fig.transFigure,
    fontsize=8,
    ha="right",
)

fig.savefig("figures/estradiol_cupid/estradiol_cupid.pdf")
