# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Thu 04 May 2023 20:28:14 BST

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
)

# mpl.rcParams["lines.linewidth"] = 0.5

DATA_DIR = Path("~/Documents/DPhil/data").expanduser()
RESULT_DIR = Path("~/Documents/DPhil/results/cupid").expanduser()
psyche_dir = DATA_DIR / "dexamethasone/1004/pdata/1"
estimator_path = RESULT_DIR / "dexamethasone/estimator_postedit"

estimator = ne.Estimator2DJ.from_pickle(estimator_path)
estimator._results[6].region = (None, (estimator._results[6].region[1][0], 1540.))

thold = (estimator.sw()[1] / estimator.default_pts[1])
print(thold)
colors = mpl.rcParams["axes.prop_cycle"].by_key()["color"]
fig, axs = estimator.plot_result(
    multiplet_thold=thold,
    region_unit="ppm",
    contour_base=3.e3,
    contour_nlevels=10,
    contour_factor=1.8,
    contour_color="k",
    contour_lw=0.1,
    multiplet_colors=colors,
    marker_size=3.,
    multiplet_show_45=False,
    multiplet_show_center_freq=True,
    axes_bottom=0.062,
    axes_left=0.035,
    axes_right=0.995,
    axes_top=0.985,
    xaxis_label_height=0.01,
    jres_sinebell=True,
    ratio_1d_2d=(3., 1.),
    xaxis_ticks=[
        (1, (6.3, 6.1)),
        (2, (5.4, 5.2, 5.0, 4.8, 4.6, 4.4)),
        (3, (4.1,)),
        (4, (2.95,)),
        (5, (2.65,)),
        (6, (2.35,)),
        (7, (2.1,)),
        (8, (1.8, 1.6, 1.4)),
        (9, (1, 0.8)),
    ],
    figsize=(9, 4.5),
)
fig.texts[0].set_fontsize(8)

fix_linewidths(axs, mpl.rcParams["lines.linewidth"])

naxes = len(axs[0])
zero_t1_data = []
cupid_data = []
for i, ax in enumerate(axs[0]):
    if i in (0, naxes - 1):
        idx = -2
    else:
        idx = -3

    cupid_line = ax.lines.pop(idx)
    cupid_data.append((cupid_line.get_xdata(), cupid_line.get_ydata()))
    zero_t1_line = ax.lines.pop(idx)
    zero_t1_data.append((zero_t1_line.get_xdata(), zero_t1_line.get_ydata()))

psyche_factor = 160
psyche = ne.Estimator1D.new_bruker(psyche_dir)
psyche_shifts, = psyche.get_shifts(unit="ppm")
psyche_spectrum = psyche_factor * psyche.spectrum.real

patch_locs = (1.4e5, 3.3e5, 7.4e5, 1.35e6)
zerot1_shift = -6.7e5
cupid_shift = -1.36e6
psyche_shift = 7.9e5

for i, ax in enumerate(axs[0]):
    left, right = ax.get_xlim()
    _, top = ax.get_ylim()

    for j, loc in enumerate(patch_locs):
        patch = mpl.patches.Rectangle(
            (left, loc),
            right - left,
            top - loc,
            facecolor="white",
            edgecolor="none",
            zorder=50 + 3 * j,
        )
        ax.add_patch(patch)

    zerot1_x, zerot1_y = zero_t1_data[i]
    ax.plot(
        zerot1_x,
        zerot1_y + zerot1_shift,
        color="k",
        zorder=52,
    )
    cupid_x, cupid_y = cupid_data[i]
    ax.plot(
        cupid_x,
        cupid_y + cupid_shift,
        color="k",
        zorder=55,
    )
    psyche_slice = slice(*psyche.convert([ax.get_xlim()], "ppm->idx")[0])
    ax.plot(
        psyche_shifts[psyche_slice],
        psyche_spectrum[psyche_slice] + psyche_shift,
        color="k",
        zorder=58,
    )

    ax.set_ylim(top=1.47e6)

    breaklines = [
        line for line in ax.get_lines()
        if np.array_equal(line.get_ydata(), np.array([1]))
    ]
    for line in breaklines:
        line.set_zorder(100)
        line.set_linewidth(mpl.rcParams["axes.linewidth"])

    ax.spines["top"].set_zorder(100)

fig.text(
    0.052,
    0.7,
    f"\\times {psyche_factor}",
    transform=fig.transFigure,
    fontsize=6,
    zorder=100,
)

panel_labels(fig, 0.04, (0.96, 0.635, 0.47, 0.38, 0.29))

axs[1][0].set_yticks([20, 10, 0, -10, -20])

xs, ys, ss = get_pure_shift_labels(estimator, thold=25000)
xs[7] = (xs[7] + xs.pop(7)) / 2
xs.pop(11)
# xs[12] -= 0.01
xs[13] += 0.02
xs[14] -= 0.02
ys[7] = max(ys[7], ys.pop(7))
ys.pop(11)
ys = [8.7e5 + y for y in ys]
y_tweaks = [
    1.92e5,  # A
    1.48e5,  # B
    1.7e5,  # C
    -1e4,  # D
    1.35e5,  # E
    1.25e5,  # F
    2e4,  # G
    1e4,  # H
    0,  # I
    -8e4,  # J
    -4e4,  # K
    -5e4,  # L
    -8e4,  # M
    -5.5e4,  # N
    -6e4,  # O
    -5e4,  # P
    -7e4,  # Q
    -2.9e5,  # R
    -5e4,  # S
    -1e5,  # T
    -5e4,  # U
    -2.25e5,  # V
    7.4e4,  # W
]
ys = [y + y_tweak for y, y_tweak in zip(ys, y_tweaks)]
ss.pop()
add_pure_shift_labels(axs[0], xs, ys, ss, fs=6)

breaks = {
    "ax": [],
    "x": [],
    "color": [],
}
patch_locs = (1.4e5, 3.3e5, 7.4e5, 1.35e6)
idx = 0
for ax_idx, ax_col in enumerate(axs.T):
    lines = [line for line in ax_col[1].get_lines() if line.get_color() != "k"]
    for line in lines:
        x = line.get_xdata()[0]
        color = line.get_color()

        if idx in [4, 22, 26, 27]:
            breaks["ax"].append(ax_idx)
            breaks["x"].append(x)
            breaks["color"].append(color)

        line.remove()
        bottoms = [ax_col[0].get_ylim()[0]] + list(patch_locs)
        tops = list(patch_locs) + [ax_col[0].get_ylim()[1]]
        for i, (bottom, top) in enumerate(zip(bottoms, tops)):
            extra = 0 if i != 3 else 0
            ax_col[0].plot([x, x], [bottom, top], color=color, lw=0.3, zorder=48 + 3 * i)
        ax_col[0].set_ylim(top=tops[-1])
        con = ConnectionPatch(
            xyA=(x, ax_col[1].get_ylim()[0]),
            xyB=(x, ax_col[0].get_ylim()[0]),
            coordsA="data",
            coordsB="data",
            axesA=ax_col[1],
            axesB=ax_col[0],
            color=color,
            lw=0.3,
            zorder=-1,
        )
        ax_col[1].add_patch(con)
        idx += 1


x_offset = 0.02
y_offset = 5e3
shift = 1e4
for i, (ax_idx, x, color) in enumerate(zip(breaks["ax"], breaks["x"], breaks["color"])):
    xs = [x + x_offset, x - x_offset]
    axs[0, ax_idx].plot(
        xs,
        [1.4e5 - y_offset, 1.4e5 + y_offset],
        color=color,
        zorder=1000,
    )
    axs[0, ax_idx].plot(
        xs,
        [1.4e5 - y_offset + shift, 1.4e5 + y_offset + shift],
        color=color,
        zorder=1000,
    )
    axs[0, ax_idx].plot(
        xs,
        [3.3e5 - y_offset, 3.3e5 + y_offset],
        color="k",
        zorder=1000,
    )
    axs[0, ax_idx].plot(
        [x + x_offset, x - x_offset],
        [3.3e5 - y_offset + shift, 3.3e5 + y_offset + shift],
        color="k",
        zorder=1000,
    )
    if i == 0:
        continue
    axs[0, ax_idx].plot(
        xs,
        [7.4e5 - y_offset, 7.4e5 + y_offset],
        color="k",
        zorder=1000,
    )
    axs[0, ax_idx].plot(
        [x + x_offset, x - x_offset],
        [7.4e5 - y_offset + shift, 7.4e5 + y_offset + shift],
        color="k",
        zorder=1000,
    )
    axs[0, ax_idx].plot(
        xs,
        [1.35e6 - y_offset, 1.35e6 + y_offset],
        color="k",
        zorder=1000,
    )
    axs[0, ax_idx].plot(
        [x + x_offset, x - x_offset],
        [1.35e6 - y_offset + shift, 1.35e6 + y_offset + shift],
        color="k",
        zorder=1000,
    )

fig.savefig("figures/dexamethasone_cupid/dexamethasone_cupid.pdf")
