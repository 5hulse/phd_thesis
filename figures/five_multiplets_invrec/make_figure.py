# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Thu 11 May 2023 17:36:46 BST

from pathlib import Path
import re

import nmrespy as ne
import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits
import numpy as np

mpl.rcParams["axes.labelpad"] = 2

xlim = (0.1249, 0.8747)
def transfer(old_ax, new_ax, new_fig):
    dim = "3d" if isinstance(new_ax, mpl_toolkits.mplot3d.axes3d.Axes3D) else "2d"
    for child in old_ax.get_children():
        to_add = False
        if dim == "2d":
            if isinstance(child, mpl.lines.Line2D) and child.get_xdata().shape[0] != 2:
                to_add = True
                func = new_ax.add_line
            elif isinstance(child, mpl.collections.PathCollection):
                to_add = True
                func = new_ax.add_collection

        elif dim == "3d":
            if isinstance(child, mpl_toolkits.mplot3d.art3d.Line3D) and child.get_data_3d()[0].shape[0] != 2:
                to_add = True
                func = new_ax.add_line
                x, y, z = child.get_data_3d()
                length = x.shape[0]
                l, r = [int(length * x) for x in xlim]
                child.set_data_3d([arr[l:r] for arr in (x, y, z)])

        if to_add:
            child.remove()
            child.axes = new_ax
            child.figure = new_fig
            child.set_transform(new_ax.transData)
            func(child)

    coords = ["x", "y"] if dim == "2d" else ["x", "y", "z"]
    for coord in coords:
        for obj in ["ticks", "ticklabels", "label", "lim"]:
            getter = getattr(old_ax, f"get_{coord}{obj}")
            setter = getattr(new_ax, f"set_{coord}{obj}")
            setter(getter())

    if dim == "3d":
        new_ax.view_init(elev=old_ax.elev, azim=old_ax.azim)

estimator_dir = Path("~/Documents/DPhil/results/invrec/five_multiplets").expanduser()

with open(estimator_dir / "shifts_and_couplings.txt", "r") as fh:
    txt = fh.readlines()

shift_t1_regex = re.compile(r"\d: (\d+\.\d+)\n")
couplings_regex = re.compile(r"\d: (\d+\.\d+), (\d+\.\d+), (\d+\.\d+)\n")
shifts = []
widths = []
t1s = []
for i in range(3):
    run_shifts = []
    run_widths = []
    run_t1s = []
    while True:
        line = txt.pop(0)
        if "Shifts" in line:
            break
    txt.pop(0)
    while True:
        line = txt.pop(0)
        match = re.search(shift_t1_regex, line)
        if match is None:
            break
        shift = match.group(1)
        run_shifts.append(float(shift))

    shifts.append(run_shifts)

    txt.pop(0)
    while True:
        line = txt.pop(0)
        match = re.search(couplings_regex, line)
        if match is None:
            break
        width = sum([float(x) for x in match.groups()]) / 2
        run_widths.append(width)

    widths.append(run_widths)

    txt.pop(0)
    while True:
        line = txt.pop(0)
        match = re.search(shift_t1_regex, line)
        if match is None:
            break
        t1 = float(match.group(1))
        run_t1s.append(t1)

    t1s.append(run_t1s)

COLORS = mpl.rcParams["axes.prop_cycle"].by_key()["color"]
COLORS.pop(4)
color_idxs = [
    [0, 0, 1, 2, 1, 0, 0, 1, 0, 0, 1, 1, 2, 1, 2, 1, 2, 1, 0, 0, 2, 2, 2, 2] + 8 * [3] + 8 * [4],
    7 * [0] + [1, 0] + 7 * [1] + 4 * [2] + [3, 2, 3] + 3 * [2] + [4] + 4 * [3] + 4 * [4] + [3, 4, 4, 3, 4],
    5 * [0] + [1, 2, 0, 1, 1, 0, 2, 1, 0, 2, 1, 2, 1, 2,1, 2, 1, 2, 2] + 6 * [3] + [4, 3, 3] + 7 * [4],
]

colors = [
    [COLORS[i] for i in idxs]
    for idxs in color_idxs
]

fig = plt.figure(figsize=(6, 7))
bottom = 0.67
top = 0.995
left = 0.045
right = 0.98
sep = 0.035

axs = []

w = (right - left - 2 * sep) / 3
half_h = (top - bottom) / 2
for i in range(3):
    l = left + i * (w + sep)
    for j in range(1, -1, -1):
        axs.append(fig.add_axes([l, bottom + j * half_h, w, half_h]))

axs.append(
    fig.add_axes(
        [-0.1, -0.08, 1.2, 0.77],
        projection="3d",
    )
)
axs[6].set_box_aspect(aspect=(1, 1.2, 1))

for i, (cols, run_shifts, run_t1s) in enumerate(zip(colors, shifts, t1s)):
    path = estimator_dir / f"estimators/estimator_{i}"
    estimator = ne.EstimatorInvRec.from_pickle(path)
    fig_tmp, axs_tmp = estimator.plot_oscs_vs_fits(
        y_range=(0.7, 5.3),
        oscillator_line_kwargs={"linewidth": 0.4},
        spectrum_line_kwargs={"linewidth": 0.8},
        oscillator_colors=cols,
        xaxis_unit="ppm",
        label_peaks=False,
        contour_base=10.e-5,
        contour_factor=1.8,
        contour_nlevels=10,
        scale=16.,
        contour_kwargs={"linewidths": 0.1},
    )
    for j in range(2):
        transfer(axs_tmp[j], axs[2 * i + j], fig)

    run_shifts = np.array(run_shifts)
    idxs = np.argsort(run_shifts)
    run_shifts = run_shifts[idxs]
    run_t1s = np.array(run_t1s)[idxs]
    for idx, (shift, t1) in enumerate(zip(run_shifts, run_t1s)):
        axs[2 * i + 1].axhline(
            t1,
            color=COLORS[idx],
            ls=":",
        )

    if i == 2:
        fig_tmp, axs_tmp = estimator.plot_result(
            azim=0.,
            elev=80.,
            oscillator_colors=cols,
            xaxis_unit="ppm",
            oscillator_line_kwargs={"linewidth": 0.4},
            spectrum_line_kwargs={"linewidth": 0.8},
        )
        transfer(axs_tmp, axs[6], fig)

for i in range(7):
    axs[i].set_xlim(reversed(xlim))
for i in (3, 5):
    axs[i].set_yticklabels([])
    axs[i].set_ylabel("")

axs[6].set_axis_off()

yax = (-0.25, 4)
axs[6].plot([xlim[1], xlim[1]], yax, [0, 0], color="k")
axs[6].plot([xlim[0], xlim[0]], yax, [0, 0], color="k")
axs[6].plot(xlim, [yax[0], yax[0]], [0, 0], color="k")
for i in range(11):
    t = 0.4 * i
    axs[6].plot([xlim[1] + 0.01, xlim[1]], [t, t], [0, 0], color="k")
    axs[6].text(xlim[1] + 0.02, t, 0, f"{t:.1f}", fontsize=7, ha="right", va="center")
for (x, lab) in zip(axs[6].get_xticks()[:-1], axs[6].get_xticklabels()[:-1]):
    axs[6].plot([x, x], [yax[0] - 0.03, yax[0]], [0, 0], color="k")
    axs[6].plot([x, x], [yax[0], yax[1]], [0, 0], color="k", ls=":", lw=0.6, zorder=-1)
    axs[6].text(x, yax[0] - 0.12, 0, lab.get_text(), fontsize=7, ha="center")
axs[6].text(xlim[1] + 0.11, 2., 0., "$\\tau (\\unit{\\second})$", va="center", fontsize=8)
axs[6].text((xlim[0] + xlim[1]) / 2, yax[0] - 0.2, 0., "¹H (ppm)", ha="center", va="center", fontsize=8)

xs = (0.05, 0.05, 0.2)
ys = (0.98, 0.815, 0.58)
for i, (x, y) in enumerate(zip(xs, ys)):
    fig.text(x, y, f"\\textbf{{{chr(97 + i)}.}}")
    axs[2 * i].text(0.985, 0.98, f"\\textbf{{Run {i + 1}}}", ha="right", va="top", transform=axs[2 * i].transAxes, fontsize=8)

axs[3].text(0.548, 3.15, "*")
axs[3].text(0.545, 1.2, "*")
axs[5].text(0.49, 5., "*", ha="center", va="center")

for i in range(1, 3, 5):
    axs[i].set_xlabel("")

fig.savefig("figures/five_multiplets_invrec/five_multiplets_invrec.pdf")
