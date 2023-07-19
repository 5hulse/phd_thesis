# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Mon 17 Jul 2023 17:18:30 BST

from utils import transfer
import nmrespy as ne
import matplotlib as mpl
import matplotlib.pyplot as plt


colors = mpl.rcParams["axes.prop_cycle"].by_key()["color"]
estimator = ne.Estimator1D.from_pickle("~/Documents/DPhil/results/onedim/cyclosporin/estimator")
# fig, axs = estimator.plot_result(
#     xaxis_unit="ppm",
#     figsize=(6, 4),
# )
# axs = axs[0]
# axs[0].set_yticks(axs[0].get_ylim())
# fig.savefig("figures/cyclosporin/cyclosporin.pdf")
# exit()



data_ylim = (3e3, 2.3e5)
mpm_ylim = (-2e4, 3.4e5)
nlp_ylim = (-1.5e4, 3.1e5)

colors.append("#808080")

mpm_cols = [
    colors[i] for i in [
        -1, 5, -1, 5, 5, 4, 5, 4, 4, 5, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3,
        2, -1, 2, 2, 2, 2, -1, 2, -1, -1, -1,
        -1, -1, 1, 1, 1, 1, -1, 0, 0, 0, 0, -1, -1,
    ]
]
nlp_cols = [
    colors[i] for i in [
        -1, 5, 5, -1, 5, 4, 5, 4, 4, 5, 4, 4, 4, 4, 4, 3, 3, 3, 3,
        2, -1, 2, 2, 2, 2, -1, -1, 2, -1, -1,
        1, 1, 1, 1, 0, 0, 0, 0,
    ]
]

kwargs = dict(
    xaxis_unit="ppm",
    axes_region_separation=0.0075,
    oscillator_line_kwargs={"lw": 0.5},
    xaxis_ticks=[
        [0, [5.52 - 0.02 * i for i in range(5)]],
        [1, [5.28 - 0.02 * i for i in range(5)]],
        [2, [5.02 - 0.02 * i for i in range(8)]],
    ],
    figsize=(6., 4.),
)
_, data_init_axs = estimator.plot_result(**kwargs)
data_init_axs = list(data_init_axs[0])
for ax in data_init_axs:
    for line in ax.lines:
        if line.get_color() != "#000000":
            line.remove()

nlp_params = []
for i, result in enumerate(estimator._results):
    nlp_params.append(result.params)
    estimator._results[i].params = result.trajectory[0]
fig, mpm_init_axs = estimator.plot_result(oscillator_colors=mpm_cols, **kwargs)
_, mpm_init_axs_ins = estimator.plot_result(oscillator_colors=mpm_cols, **kwargs)
mpm_init_axs = list(mpm_init_axs[0])
mpm_init_axs_ins = list(mpm_init_axs_ins[0])

for i, para in enumerate(nlp_params):
    estimator._results[i].params = nlp_params[i]

fig, nlp_init_axs = estimator.plot_result(oscillator_colors=nlp_cols, **kwargs)
fig, nlp_init_axs_ins = estimator.plot_result(oscillator_colors=nlp_cols, **kwargs)
nlp_init_axs = list(nlp_init_axs[0])
nlp_init_axs_ins = list(nlp_init_axs_ins[0])
bboxes = [ax.bbox._bbox for ax in nlp_init_axs]


for ax in mpm_init_axs + nlp_init_axs:
    for line in ax.lines:
        if line.get_color() == "#000000":
            line.remove()

fig = plt.figure(figsize=(6, 4.5))

top = 0.985
bottom = 0.07
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

fig.text((left + right) / 2, 0.01, "\\textsuperscript{1}H (ppm)", ha="center")

l, w, h = 0.25, 0.2, 0.14
insaxs = []
insaxs.append(fig.add_axes([l, 0.61, w, h]))
insaxs.append(fig.add_axes([l, 0.26, w, h]))

transfer(mpm_init_axs_ins[1], insaxs[0], fig)
transfer(nlp_init_axs_ins[1], insaxs[1], fig)
for ax in insaxs:
    for line in ax.lines:
        if line.get_color() == "#000000":
            line.remove()

ins_xlim = (5.27, 5.2)
height = 1.3e4
b0, b1 = 2e4, 8e3
insaxs[0].set_ylim(b0, b0 + height)
insaxs[1].set_ylim(b1, b1 + height)
for insax in insaxs:
    insax.set_xlim(ins_xlim)

for b, ax, insax in zip((b0, b1), (mpm_axs[1], nlp_axs[1]), insaxs):
    rect = mpl.patches.Rectangle(
        xy=(ins_xlim[0], b),
        width=ins_xlim[1] - ins_xlim[0],
        height=height,
        edgecolor="w",
        facecolor="none",
        lw="0.8",
        zorder=500,
    )
    ax.add_patch(rect)
    rect = mpl.patches.Rectangle(
        xy=(ins_xlim[0], b),
        width=ins_xlim[1] - ins_xlim[0],
        height=height,
        edgecolor="k",
        facecolor="none",
        lw="0.4",
        zorder=1500,
    )
    ax.add_patch(rect)

    bbox = insax.bbox._bbox
    xyAs = [
        (ins_xlim[0], b),
        (ins_xlim[1], b + height),
    ]
    xyBs = [
        (bbox.x0, bbox.y0),
        (bbox.x1, bbox.y1),
    ]
    for xyA, xyB in zip(xyAs, xyBs):
        connect = mpl.patches.ConnectionPatch(
            xyA=xyA,
            xyB=xyB,
            coordsA="data",
            coordsB="figure fraction",
            axesA=ax,
            zorder=1000,
            lw=0.8,
            color="w"
        )
        ax.add_patch(connect)
        connect = mpl.patches.ConnectionPatch(
            xyA=xyA,
            xyB=xyB,
            coordsA="data",
            coordsB="figure fraction",
            axesA=ax,
            zorder=2000,
            lw=0.4,
        )
        ax.add_patch(connect)

label_axs = (0, 0, 1, 2, 2, 2)
label_xs = (0.29, 0.69, 0.59, 0.235, 0.54, 0.725)
label_ys = (0.57, 0.55, 0.63, 0.8, 0.75, 0.88)
label_ss = [f"({chr(i + 65)})" for i in range(6)]

for i, (ax_idx, x, y, s) in enumerate(zip(label_axs, label_xs, label_ys, label_ss)):
    data_axs[ax_idx].text(
        x, y, s,
        ha="center", transform=data_axs[ax_idx].transAxes, color=colors[i],
        fontsize=9,
    )

for i, ax in enumerate((data_axs[0], mpm_axs[0], nlp_axs[0])):
    ax.text(0.02, 0.9, f"\\textbf{{{chr(i + 97)}.}}", va="bottom", transform=ax.transAxes)

fig.savefig("figures/cyclosporin/cyclosporin.pdf")
