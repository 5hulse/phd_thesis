# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Thu 04 Jan 2024 14:22:28 GMT

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
    get_pure_shift_labels,
    add_pure_shift_labels,
)

mpl.rcParams["lines.linewidth"] = 0.6
RESULT_DIR = Path("~/Documents/DPhil/results").expanduser()

result_dir = RESULT_DIR / "cupid/strychnine"
estimator_path = result_dir / "estimator"
estimator_dec_path = result_dir / "estimator_dec"
save_path = "figures/strychnine_cupid/strychnine_cupid.pdf"

estimator = ne.Estimator2DJ.from_pickle(estimator_path)
estimator.predict_multiplets(rm_spurious=True, max_iterations=1, check_neg_amps_every=1)
reg0 = estimator._results[0].region
estimator._results[0].region = (None, (reg0[1][0] + 15., reg0[1][1]))

colors = mpl.rcParams["axes.prop_cycle"].by_key()["color"]
sws = estimator.sw()
pts = estimator.default_pts

fig, axs = estimator.plot_result(
    multiplet_vertical_shift=0.,
    axes_left=0.05,
    axes_right=0.99,
    axes_bottom=0.07,
    axes_top=0.985,
    xaxis_label_height=0.01,
    ratio_1d_2d=(2.5, 1.),
    contour_base=0.6,
    contour_nlevels=12,
    contour_factor=2.2,
    contour_color="k",
    contour_lw=0.1,
    multiplet_colors=colors,
    marker_size=2.,
    figsize=(6., 4.),
    region_unit="ppm",
    axes_region_separation=0.05,
    multiplet_show_45=False,
    xaxis_ticks=[
        (0, [8.1]),
        (1, [7.3, 7.2, 7.1]),
        (2, [5.9]),
        (3, [4.3, 4.2, 4.1, 4.0, 3.9, 3.8, 3.7]),
        (4, [3.2, 3.1]),
        (5, [2.9, 2.8, 2.7, 2.6]),
        (6, [2.4, 2.3]),
        (7, [1.9, 1.8]),
        (8, [1.5, 1.4, 1.3, 1.2]),
    ],
)

for ax in axs[1]:
    ax.set_ylim(20, -20)

xlabel = fig.texts[0]
xlabel.set_fontsize(8)

fix_linewidths(axs, mpl.rcParams["lines.linewidth"])
panel_labels(fig, 0.053, (1., 0.945, 0.52, 0.42, 0.32))
del fig.texts[1]

estimator_dec = ne.Estimator2DJ.from_pickle(estimator_dec_path)
dec_fid = estimator_dec.data[0]
dec_fid = ne.sig.exp_apodisation(dec_fid, 4.)
dec_fid[0] *= 0.5
dec_spec = ne.sig.ft(dec_fid).real
dec_spec = estimator_dec.spectrum_direct
shifts = estimator_dec.get_shifts(unit="ppm", meshgrid=False)[1]
n_axs = len(axs[0])
for i, ax in enumerate(axs[0]):
    idx = -2 if i in (0, n_axs - 1) else -3
    ax.plot(
        shifts + 0.03,
        dec_spec + ax.get_lines()[idx].get_ydata()[0],
        color="#a0a0a0",
        zorder=1,
    )

with open(result_dir.parent.parent / "invrec/strychinine/shifts.pkl", "rb") as fh:
    chem_shifts = pickle.load(fh)

# ss = [
#     f"({chr(i[0] + 65)})"
#     for i in reversed(sorted(enumerate(chem_shifts), key=lambda x:x[1]))
# ]
ss = [f"({chr(i + 65)})" for i in range(22)]
ss[-4] = "(S,T)"
del ss[-3]
xs, ys, _ = get_pure_shift_labels(estimator, yshift=110., thold=10.)
# xs[11] += 0.02
xs[12] += 0.0263
ys[12] += 32.
xs[13] -= 0.03
ys[13] += 5.

most_xs, most_ys, most_ss = xs[:12] + xs[13:], ys[:12] + ys[13:], ss[:12] + ss[13:]

add_pure_shift_labels(axs[0], most_xs, most_ys, most_ss, fs=7)
add_pure_shift_labels(axs[0], [xs[12]], [ys[12]], [ss[12]], fs=6.7)

raise_axes(axs, 1.08)

for ax_col in axs.T:
    lines = [line for line in ax_col[1].get_lines() if line.get_color() != "k"]
    for line in lines:
        x = line.get_xdata()[0]
        color = line.get_color()
        line.remove()
        for ax in ax_col:
            ax.axvline(x, color=color, ls="-", zorder=-1, lw=0.5)

fig.savefig(save_path)
