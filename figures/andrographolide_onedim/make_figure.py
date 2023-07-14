# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Wed 12 Jul 2023 14:58:28 BST

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
print(estimator.sfo)
exit()

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
    figsize=(9, 5),
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

fig = plt.figure(figsize=(9, 5))

top = 0.95
bottom = 0.065
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
fig.text((left + right) / 2, 0.01, "\\textsuperscript{1}H (ppm)", ha="center")

label_axs = (0, 1, 2, 3, 4, 5, 5, 5, 6, 7, 7, 8)
label_xs = (0.46, 0.48, 0.76, 0.73, 1.05, 0.42, 0.65, 0.82, 0.45, 0.31, 0.735, 0.52)
label_ys = (0.65, 0.46, 0.85, 0.85, 0.3, 0.2, 0.5, 0.4, 0.4, 0.35, 0.43, 0.36)
labels_ss = ("(A)", "(B)", "(E)", "(F)", "EtOH", "H\\textsubscript{2}O", "(K)", "(L)", "(O)", "(P)", "(Q)", "(V)")

for ax_idx, x, y, s in zip(label_axs, label_xs, label_ys, labels_ss):
    data_axs[ax_idx].text(x, y, s, ha="center", transform=data_axs[ax_idx].transAxes)

for i, ax in enumerate((data_axs[0], mpm_axs[0], nlp_axs[0])):
    ax.text(0.03, 0.97, f"\\textbf{{{chr(i + 97)}.}}", va="top", transform=ax.transAxes)

# inset_l, inset_r = 0.1, 0.5
# bbox4_x0, bbox4_x1 = bboxes[4].x0, bboxes[4].x1
# l, r = [(bbox4_x1 - bbox4_x0) * x + bbox4_x0 for x in (inset_l, inset_r)]
# inset_xlim = (3.47, 3.41)
# inset_ylims = [(2e3, 5e3), (1e3, 1e4), (1e3, 1e4)]
# for ax, inset_ylim in zip(
#     (data_axs[4], mpm_axs[4], nlp_axs[4]),
#     inset_ylims,
# ):
#     inax = fig.add_axes([l, 0.4, r - l, 0.5], transform=ax.transAxes)


# def cp_lines(from_ax, to_ax):
#     for line in from_ax.lines:
#         to_ax.plot(
#             line.get_xdata(),
#             line.get_ydata(),
#             color=line.get_color(),
#             lw=line.get_lw(),
#         )

# inax.set_xlim(inset_xlim)
# inax.set_ylim(inset_ylim)
# ax.indicate_inset_zoom(inax, edgecolor="k")



# ### OLD CODE FOR 2 LINES OF FIGURES

# top_span = bboxes[3].x1 - bboxes[0].x0
# bottom_span = bboxes[-1].x1 - bboxes[4].x0

# fig = plt.figure(figsize=(6, 8))

# top = 0.99
# bottom = 0.05
# left = 0.01
# right = 0.99
# hspace = 0.02

# if top_span > bottom_span:
#     top_width = right - left
#     bottom_width = top_width * (bottom_span / top_span)
#     top_left = left
#     top_right = right
#     center = (right + left) / 2
#     bottom_left = center - (bottom_width / 2)
#     bottom_right = bottom_left + bottom_width

# else:
#     bottom_width = right - left
#     top_width = bottom_width * (top_span / bottom_span)
#     bottom_left = left
#     bottom_right = right
#     center = (right + left) / 2
#     top_left = center - (top_width / 2)
#     top_right = top_left + top_width

# height = (top - bottom - hspace) / 2
# width = right - left

# data_axs = []
# mpm_axs = []
# nlp_axs = []

# # Top
# relative_widths_top = [bbox.x1 - bbox.x0 for bbox in bboxes[:4]]
# relative_widths_top_sum = float(sum(relative_widths_top))
# relative_widths_top = [rwt / relative_widths_top_sum for rwt in relative_widths_top]
# widths_top = [w * (top_right - top_left - 3 * kwargs["axes_region_separation"]) for w in relative_widths_top]
# for i in range(4):
#     left = top_left + sum(widths_top[:i]) + i * kwargs["axes_region_separation"]
#     width = widths_top[i]
#     data_axs.append(fig.add_axes([left, top - height / 3, width, height / 3]))
#     mpm_axs.append(fig.add_axes([left, top - 2 * height / 3, width, height / 3]))
#     nlp_axs.append(fig.add_axes([left, top - height, width, height / 3]))

# # Bottom
# relative_widths_bottom = [bbox.x1 - bbox.x0 for bbox in bboxes[4:]]
# relative_widths_bottom_sum = float(sum(relative_widths_bottom))
# relative_widths_bottom = [rwb / relative_widths_bottom_sum for rwb in relative_widths_bottom]
# widths_bottom = [w * (bottom_right - bottom_left - 2 * kwargs["axes_region_separation"]) for w in relative_widths_bottom]
# for i in range(3):
#     left = bottom_left + sum(widths_bottom[:i]) + i * kwargs["axes_region_separation"]
#     width = widths_bottom[i]
#     data_axs.append(fig.add_axes([left, bottom + 2 * height / 3, width, height / 3]))
#     mpm_axs.append(fig.add_axes([left, bottom + height / 3, width, height / 3]))
#     nlp_axs.append(fig.add_axes([left, bottom, width, height / 3]))

# for old, new in zip(data_init_axs, data_axs):
#     transfer(old, new, fig)
# for old, new in zip(mpm_init_axs, mpm_axs):
#     transfer(old, new, fig)
# for old, new in zip(nlp_init_axs, nlp_axs):
#     transfer(old, new, fig)

# for ax in mpm_axs + nlp_axs + data_axs:
#     ax.set_ylim(ylim)

# estimator.edit_result(
#     index=-4,
#     split_oscs={
#         0: {"number": 3, "separation": 2.},
#         1: {"number": 3, "separation": 2.},
#         2: {"number": 3, "separation": 2.},
#         3: {"number": 3, "separation": 2.},
#         # 4: {"number": 4, "separation": 1.},
#         # 5: {"number": 4, "separation": 1.},
#         # 6: {"number": 6, "separation": 1.},
#         # 7: {"number": 6, "separation": 1.},
#         # 8: {"number": 4, "separation": 1.},
#         # 9: {"number": 4, "separation": 1.},
#     },
#     initial_trust_radius=0.2,
#     max_iterations=500,
#     check_neg_amps_every=100,
#     save_trajectory=True,
# )
# estimator._results[-4].params = estimator._results[-4].trajectory[0]

arrow_span = (
    left + sum(widths[:4]) + kwargs["axes_region_separation"] * 4,
    left + sum(widths[:6]) + kwargs["axes_region_separation"] * 5,
)
arrow_height = 0.975
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
