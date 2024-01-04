# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Thu 04 Jan 2024 00:17:43 GMT

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

RESULT_DIR = Path("~/Documents/DPhil/results/cupid/camphor").expanduser()

# === CONFIGURATION ===
estimator_path = RESULT_DIR / "estimator"

residual_shift = 3e6
multiplet_shift = residual_shift + 1e6
onedim_shift = multiplet_shift + 1e6
ax_top = onedim_shift + 3e6

# =====================
estimator = ne.Estimator2DJ.from_pickle(estimator_path)
thold = 2. * (estimator.sw()[1] / estimator.default_pts[1])
estimator.predict_multiplets(indices=[1, 2, 3, 4, 5], thold=thold, rm_spurious=True, max_iterations=1, check_neg_amps_every=1)
colors = mpl.rcParams["axes.prop_cycle"].by_key()["color"]
colors.append("#808080")
mp_colors = [
    colors[i] for i in [
        0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 0, -1, -1, -1, 1,
    ]
]
fig, axs = estimator.plot_result(
    axes_right=0.99,
    axes_bottom=0.07,
    axes_top=0.98,
    axes_left=0.05,
    multiplet_thold=thold,
    region_unit="ppm",
    contour_base=1.5e4,
    contour_nlevels=10,
    contour_factor=1.6,
    contour_color="k",
    contour_lw=0.1,
    multiplet_colors=mp_colors,
    marker_size=2.,
    multiplet_show_45=False,
    jres_sinebell=True,
    xaxis_label_height=0.015,
    xaxis_ticks=[
        (0, (2.54, 2.5)),
        (1, (2.34, 2.3, 2.26)),
        (2, (2.08, 2.04)),
        (3, (1.94, 1.9, 1.86, 1.82, 1.78)),
        (4, (1.68, 1.64)),
        (5, (1.36, 1.32, 1.28, 1.24)),
    ],
    ratio_1d_2d=(2.2, 1.),
    figsize=(6, 4.5),
)
fig.texts[0].set_fontsize(8)
axs[1, 0].set_yticks([-20, -10, 0, 10, 20])

fix_linewidths(axs, 0.8)
panel_labels(fig, 0.055, (1., 0.95, 0.72, 0.49, 0.41, 0.32))
del fig.texts[1]
print(fig.texts)

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

raise_axes(axs, 1.76)
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
xs, ys, ss = get_pure_shift_labels(estimator, yshift=1.05e6, thold=4e4)
xs = [xs[2]] + xs[5:]
ys = [ys[2]] + ys[5:]
ss = ["DMSO"] + ss[:-5]
add_pure_shift_labels(axs[0], xs, ys, ss, fs=7)

tilt_dataset = BrukerDataset(Path("~/Documents/DPhil/data/camphor/2/pdata/999").expanduser())
tilt_spectrum = tilt_dataset.data
tilt_shifts, = tilt_dataset.get_samples()

scale = 10
yshift = 2.25e6
for ax in axs[0]:
    ax.plot(tilt_shifts + 0.25 * (estimator.sw(unit="ppm")[1] / estimator.default_pts[1]), (scale * tilt_spectrum) + yshift, color="k")
fig.text(0.96, 0.75, f"$\\times {scale}$", fontsize=7)

axs[0, 1].text(2.253, 2.215e6, "*")
axs[0, 1].text(2.27, 2.215e6, "*")
axs[0, 1].text(2.353, 2.215e6, "*")
axs[0, 3].text(1.915, 2.22e6, "*")
axs[0, 3].text(1.903, 2.26e6, "*")
axs[0, 3].text(1.873, 2.26e6, "*")
axs[0, 3].text(1.855, 2.22e6, "*")
axs[0, 3].text(1.755, 2.22e6, "*")
axs[0, 4].text(1.675, 2.3e6, "*")
axs[0, 4].text(1.69, 2.23e6, "*")
axs[0, 4].text(1.645, 2.22e6, "*")
axs[0, 5].text(1.3175, 2.4e6, "*")
axs[0, 5].text(1.3017, 2.43e6, "*")
axs[0, 5].text(1.286, 2.39e6, "*")

fig.savefig("figures/camphor_cupid/camphor_cupid.pdf")
