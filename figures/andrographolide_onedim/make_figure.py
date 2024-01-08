# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Mon 08 Jan 2024 10:48:34 EST

import copy
from utils import transfer
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch
import matplotlib as mpl
import nmrespy as ne
import numpy as np


colors = mpl.rcParams["axes.prop_cycle"].by_key()["color"]

path = Path("~/Documents/DPhil/results/onedim/andrographolide/estimator_new_new").expanduser()
estimator = ne.Estimator1D.from_pickle(path)

estimator._results[0].region = estimator.convert(((6.665, 6.59),), "ppm->hz")
estimator._results[1].region = estimator.convert(((5.75, 5.665),), "ppm->hz")
estimator._results[2].region = estimator.convert(((4.85, 4.795),), "ppm->hz")
estimator._results[4].region = estimator.convert(((4.67, 4.6),), "ppm->hz")
estimator._results[5].region = estimator.convert(((1.99, 1.83),), "ppm->hz")
estimator._results[6].region = estimator.convert(((2.36, 2.29),), "ppm->hz")
estimator._results[7].region = estimator.convert(((1.415, 1.305),), "ppm->hz")

estimator._results.append(copy.deepcopy(estimator._results[3]))
estimator._results[3].region = estimator.convert(((3.4, 3.2),), "ppm->hz")
estimator._results[-1].region = estimator.convert(((3.473, 3.42),), "ppm->hz")
estimator_params_last = copy.deepcopy(estimator._results[-1].params)
estimator._results[-1].params = np.array([[0.0001, 0., 0., 1.]])

data_ylim = (3e3, 3.7e4)
mpm_ylim = (-3.e3, 4e4)
nlp_ylim = (-2.5e3, 4e4)

mpm_cols = [
    colors[i] for i in [
        1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1,
        1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1,
    ]
]

kwargs = dict(
    xaxis_unit="ppm",
    xaxis_ticks=[
        [0, (6.65, 6.61)],
        [1, (5.74, 5.71, 5.68)],
        [2, (4.84, 4.81)],
        [3, (4.65, 4.62)],
        [4, (3.46, 3.43)],
        [5, [3.39 - i * 0.03 for i in range(7)]],
        [6, (2.34, 2.31)],
        [7, [1.97 - i * 0.03 for i in range(5)]],
        [8, (1.4, 1.37, 1.34, 1.31)],
    ],
    axes_region_separation=0.005,
    oscillator_line_kwargs={"lw": 0.5},
    residual_shift=3.5e3,
    figsize=(9.78, 3.1),
    xaxis_label_height=0.015,
)
_, data_init_axs = estimator.plot_result(**kwargs)
data_init_axs = list(data_init_axs[0])
for ax in data_init_axs:
    for line in ax.lines:
        if line.get_color() != "#000000":
            line.remove()

nlp_params = []
for i in range(8):
    nlp_params.append(estimator._results[i].params)
    estimator._results[i].params = estimator._results[i].trajectory[0]
_, mpm_init_axs = estimator.plot_result(oscillator_colors=mpm_cols, **kwargs)
mpm_init_axs = list(mpm_init_axs[0])

for i in range(8):
    estimator._results[i].params = nlp_params[i]

estimator.edit_result(0, rm_oscs=[2], save_trajectory=True)

fig, nlp_init_axs = estimator.plot_result(oscillator_colors=colors[1], **kwargs)
nlp_init_axs = list(nlp_init_axs[0])
bboxes = [ax.bbox._bbox for ax in nlp_init_axs]


for ax in mpm_init_axs + nlp_init_axs:
    for line in ax.lines:
        if line.get_color() == "#000000":
            line.remove()

fig = plt.figure(figsize=kwargs["figsize"])

top = 0.93
bottom = 0.1
left = 0.005
right = 0.995

data_axs = []
mpm_axs = []
nlp_axs = []

relative_widths = [bbox.x1 - bbox.x0 for bbox in bboxes]
relative_widths_sum = float(sum(relative_widths))
relative_widths = [rwt / relative_widths_sum for rwt in relative_widths]
widths = [
    w * (right - left - (len(relative_widths) - 1) * kwargs["axes_region_separation"])
    for w in relative_widths
]
ylims = (data_ylim, mpm_ylim, nlp_ylim)
h_sum = sum([ylim[1] - ylim[0] for ylim in ylims])
heights = [((ylim[1] - ylim[0]) / h_sum) * (top - bottom) for ylim in ylims]

break_kwargs = {
    "marker": [(-1, -3), (1, 3)],
    "markersize": 8,
    "linestyle": "none",
    "color": "k",
    "mec": "k",
    "mew": 0.8,
    "clip_on": False,
}
for i in range(len(widths)):
    ax_left = left + sum(widths[:i]) + i * kwargs["axes_region_separation"]
    width = widths[i]
    data_axs.append(
        fig.add_axes([ax_left, top - heights[0], width, heights[0]])
    )
    mpm_axs.append(
        fig.add_axes([ax_left, top - sum(heights[:2]), width, heights[1]])
    )
    nlp_axs.append(
        fig.add_axes([ax_left, top - sum(heights[:3]), width, heights[2]])
    )
    transfer(data_init_axs[i], data_axs[i], fig)
    transfer(mpm_init_axs[i], mpm_axs[i], fig)
    transfer(nlp_init_axs[i], nlp_axs[i], fig)

    data_axs[-1].set_ylim(data_ylim)
    mpm_axs[-1].set_ylim(mpm_ylim)
    nlp_axs[-1].set_ylim(nlp_ylim)

    for ax in (data_axs[i], mpm_axs[i], nlp_axs[i]):
        if i != 0:
            ax.spines["left"].set_visible(False)
            data_axs[i].plot([0], [1], transform=data_axs[i].transAxes, **break_kwargs)
            nlp_axs[i].plot([0], [0], transform=nlp_axs[i].transAxes, **break_kwargs)
        if i != len(widths) - 1:
            ax.spines["right"].set_visible(False)
            data_axs[i].plot([1], [1], transform=data_axs[i].transAxes, **break_kwargs)
            nlp_axs[i].plot([1], [0], transform=nlp_axs[i].transAxes, **break_kwargs)

    data_axs[i].spines["bottom"].set_visible(False)
    mpm_axs[i].spines["top"].set_visible(False)
    mpm_axs[i].spines["bottom"].set_visible(False)
    nlp_axs[i].spines["top"].set_visible(False)

    data_axs[i].set_xticks([])
    mpm_axs[i].set_xticks([])

for ax in mpm_axs:
    for line in ax.lines:
        line.set(clip_on=False)

ax4_ylim = (3.6e3, 6e3)
ax4_scale = (ax4_ylim[1] - ax4_ylim[0]) / (data_axs[4].get_ylim()[1] - data_axs[4].get_ylim()[0])
data_axs[4].set_ylim(ax4_ylim)
bottom = -2e2
top = bottom + (ax4_scale * (mpm_axs[4].get_ylim()[1] - mpm_axs[4].get_ylim()[0]))
mpm_axs[4].set_ylim(bottom, top)
nlp_axs[4].set_ylim(bottom, top)
nlp_axs[-1].plot([0], [0], transform=nlp_axs[4].transAxes, **break_kwargs)
nlp_axs[-1].plot([1], [0], transform=nlp_axs[4].transAxes, **break_kwargs)

for ax in (mpm_axs[4], nlp_axs[4]):
    for i, line in enumerate(ax.lines):
        if i <= 4:
            pass
        elif i == 5:
            line.set_ydata(line.get_ydata() - 6.2e3)
        else:
            line.set_ydata(line.get_ydata() - 3.27e3)

data_axs[4].text(
    0.02, 0.5, f"$\\times {round(1 / ax4_scale)}$",
    transform=data_axs[4].transAxes, va="top",
    fontsize=8,
)
fig.text((left + right) / 2, kwargs["xaxis_label_height"], "\\textsuperscript{1}H (ppm)", ha="center")

label_axs = (0, 1, 2, 3, 4, 5, 5, 5, 6, 7, 7, 8)
label_xs = (0.46, 0.48, 0.76, 0.73, 1.05, 0.42, 0.65, 0.82, 0.45, 0.31, 0.735, 0.52)
label_ys = (0.66, 0.47, 0.84, 0.84, 0.3, 0.21, 0.51, 0.41, 0.41, 0.36, 0.44, 0.37)
labels_ss = ("(A)", "(B)", "(E)", "(F)", "EtOH", "H\\textsubscript{2}O", "(K)", "(L)", "(O)", "(P)", "(Q)", "(V)")

for ax_idx, x, y, s in zip(label_axs, label_xs, label_ys, labels_ss):
    data_axs[ax_idx].text(x, y, s, ha="center", transform=data_axs[ax_idx].transAxes)

for i, ax in enumerate((data_axs[0], mpm_axs[0], nlp_axs[0])):
    ax.text(0.03, 0.83, f"\\textbf{{{chr(i + 97)}.}}", va="bottom", transform=ax.transAxes)

arrow_span = (
    left + sum(widths[:4]) + kwargs["axes_region_separation"] * 4,
    left + sum(widths[:6]) + kwargs["axes_region_separation"] * 5,
)
arrow_height = 0.965
arrow = ConnectionPatch(
    xyA=(arrow_span[0], arrow_height),
    xyB=(arrow_span[1], arrow_height),
    coordsA="figure fraction",
    coordsB="figure fraction",
    arrowstyle="|-|",
    facecolor="k",
    lw=0.7,
    mutation_scale=2.,
)
data_axs[4].add_patch(arrow)
fig.text(
    (arrow_span[0] + arrow_span[1]) / 2,
    arrow_height + 0.01,
    "same estimated region",
    fontsize=7,
    ha="center",
)

fig.savefig("figures/andrographolide_onedim/andrographolide_onedim.pdf")
