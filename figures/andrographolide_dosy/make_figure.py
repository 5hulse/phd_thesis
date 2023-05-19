# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Fri 19 May 2023 18:36:05 BST

import nmrespy as ne
from utils import RESULT_DIR


estimator = ne.EstimatorDiffOneshot.from_pickle(
    RESULT_DIR / "diffusion/andrographolide/estimator_postedit",
)
# estimator.sigma = 1.
ud_regions = {
    1: ((5.800050185695225, 5.6),),
    # 6: ((3.700084219682121, 3.),),
}

for i, res in enumerate(estimator._estimators[0]._results):
    print(i, res.get_region(unit="ppm"))
# exit()

for i, region in ud_regions.items():
    for est in estimator._estimators:
        est._results[i].region = estimator.convert(region, "ppm->hz")
        est._results.pop(8)

fig, axs = estimator.plot_dosy(
    y_range=(2.e-10, 5.e-10),
    y_pts=128,
    # oscillator_colors="viridis",
    xaxis_unit="ppm",
    contour_base=8.e8,
    contour_factor=1.35,
    contour_nlevels=14,
    distribution_width=2.5,
    region_separation=0.005,
    label_peaks=False,
    xaxis_ticks=[(0, [6.65]), (1, [5.7]), (2, list([5.1 - i * 0.2 for i in range(18)])), (4, [0.65])],
    figsize=(9, 3),
    gridspec_kwargs={
        "width_ratios": [1, 6],
        "left": 0.04,
        "right": 0.995,
        "bottom": 0.089,
    },
    spectrum_line_kwargs={"linewidth": .8},
    oscillator_line_kwargs={"linewidth": .0},
    contour_kwargs={"linewidths": 0.05},
)
axs = list(axs)
ax1_box = axs[1].get_position()
distribution = axs[1].lines[0]
xdata = distribution.get_xdata()
ydata = distribution.get_ydata()
xlim = axs[1].get_xlim()
ylim = axs[1].get_ylim()
h = xlim[1] - xlim[0]
xlims = [(xlim[0] + 0.025 * h, xlim[0] + 0.06 * h), (xlim[0] + 0.94 * h, xlim[0] + 0.975 * h)]
w_ratios = [x[1] - x[0] for x in xlims]
w_ratios = [wr / sum(w_ratios) for wr in w_ratios]
ws = [(ax1_box.x1 - ax1_box.x0) * wr for wr in w_ratios]
axs[1].remove()
axs.pop(1)
space = 0.002
lefts = [ax1_box.x0 + i * space + sum(ws[:i]) for i in range(2)]


for i, (left, width) in enumerate(zip(lefts, ws)):
    axs.insert(
        1 + i,
        fig.add_axes([left, ax1_box.y0, width - space, ax1_box.y1 - ax1_box.y0]),
    )

axs[1].spines["right"].set_visible(False)
axs[2].spines["left"].set_visible(False)
axs[1].plot(xdata, ydata, color="k")
axs[1].set_xlim(xlims[0])
axs[1].set_xticks([])
axs[1].set_ylim(ylim)
axs[2].plot(xdata, ydata, color="k")
axs[2].set_xlim(xlims[1])
axs[2].set_xticks([])
axs[2].set_ylim(ylim)
axs[2].set_yticks([])


ybot, ytop = axs[0].get_ylim()
axs[0].set_ylim(0.3 * ybot, 0.22 * ytop)
axs[1].set_ylabel("$D (10^{10} \\unit{\\meter\\squared\\per\\second})$", labelpad=2)
yticks = [2e-10 + i * 5e-11 for i in range(7)]
axs[1].set_yticks(yticks)
axs[1].set_yticklabels([f"{yt * 1e10:.1f}" for yt in yticks])
axs[3].set_xlabel(axs[3].get_xlabel(), labelpad=-1)

break_kwargs = {
    "marker": [(-1, -3), (1, 3)],
    "markersize": 6,
    "linestyle": "none",
    "color": "k",
    "mec": "k",
    "mew": 1,
    "clip_on": False,
}
axs[1].plot([1], [1], transform=axs[1].transAxes, **break_kwargs)
axs[1].plot([1], [0], transform=axs[1].transAxes, **break_kwargs)
axs[2].plot([0], [1], transform=axs[2].transAxes, **break_kwargs)
axs[2].plot([0], [0], transform=axs[2].transAxes, **break_kwargs)

fig.savefig("figures/andrographolide_dosy/andrographolide_dosy.pdf")
