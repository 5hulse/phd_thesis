# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Fri 21 Jul 2023 13:06:36 BST

from pathlib import Path

from bruker_utils import BrukerDataset
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

RESULT_DIR = Path("~/Documents/DPhil/results/cupid").expanduser()

# === CONFIGURATION ===
estimator_path = RESULT_DIR / "quinine/estimator"
save_path = "figures/quinine_cupid/quinine_cupid.pdf"


residual_shift = 3e6
multiplet_shift = residual_shift + 1e6
onedim_shift = multiplet_shift + 1e6
ax_top = onedim_shift + 3e6

# =====================
colors = mpl.rcParams["axes.prop_cycle"].by_key()["color"]
colors.append("#808080")
mp_colors = [
    colors[i] for i in [
        0, 1, 2, 3, -1, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1,
    ]
]
estimator = ne.Estimator2DJ.from_pickle(estimator_path)
thold = 6. * estimator.default_multiplet_thold
estimator.predict_multiplets(thold=thold, rm_spurious=True, max_iterations=1, check_neg_amps_every=1)
fig, axs = estimator.plot_result(
    axes_right=0.985,
    axes_bottom=0.065,
    axes_top=0.975,
    axes_left=0.05,
    multiplet_thold=thold,
    region_unit="ppm",
    contour_base=1.7e4,
    contour_nlevels=10,
    contour_factor=1.6,
    contour_color="k",
    contour_lw=0.1,
    multiplet_colors=mp_colors,
    marker_size=3.,
    multiplet_show_45=False,
    jres_sinebell=True,
    xaxis_label_height=0.01,
    xaxis_ticks=[
        (0, (5.8, 5.75, 5.7, 5.65, 5.6)),
        (1, (4.95, 4.9,)),
        (2, (3.7, 3.65)),
        (3, (3.15, 3.1,)),
        (4, (2.75, 2.7, 2.65)),
        (5, (1.95, 1.9, 1.85, 1.8, 1.75)),
        (6, (1.6, 1.55, 1.5, 1.45, 1.4)),
    ],
    ratio_1d_2d=(3., 1.),
    figsize=(6, 4),
)
fig.texts[0].set_fontsize(8)
axs[1, 0].set_yticks([-20, -10, 0, 10, 20])

fix_linewidths(axs, 0.8)
panel_labels(fig, 0.053, (0.95, 0.66, 0.405, 0.34, 0.27))

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
specline = axs[0][1].plot(shifts_r1, spec_without_h2o + vshift, color="k", lw=0.8)

raise_axes(axs, 1.9)
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
        ax0.axvline(x, color=line.get_color(), lw=0.5, zorder=-1)

# Label pure shift peaks
xs, ys, ss = get_pure_shift_labels(estimator, yshift=2.1e6, thold=1e5)
xs[6] -= 0.015
ss[3] = "(D)*"
add_pure_shift_labels(axs[0], xs, ys, ss, fs=7)

tilt_dataset = BrukerDataset(Path("~/Documents/DPhil/data/quinine/2/pdata/999").expanduser())
tilt_spectrum = tilt_dataset.data
tilt_shifts, = tilt_dataset.get_samples()

scale = 10
yshift = 5.8e6
for ax in axs[0]:
    ax.plot(tilt_shifts, (scale * tilt_spectrum) + yshift, color="k")

fig.text(
    0.30,
    0.34,
    "H\\textsubscript{2}\\hspace{-0.7pt}O",
    color=colors[-1],
    fontsize=8,
)
fig.text(0.96, 0.7, f"$\\times {scale}$", fontsize=7)
axs[0, 1].text(4.935, 5.9e6, "*")

fig.savefig(save_path)
