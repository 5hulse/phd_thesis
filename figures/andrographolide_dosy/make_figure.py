# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Sun 21 May 2023 21:21:48 BST

import matplotlib as mpl
import nmrespy as ne
from utils import RESULT_DIR, fix_linewidths


region_separation = 0.005

estimator = ne.EstimatorDiffOneshot.from_pickle(
    RESULT_DIR / "diffusion/andrographolide/estimator_postedit",
)
ud_regions = {
    1: ((5.800050185695225, 5.6),),
    # 6: ((3.700084219682121, 3.1),),
    # 7: ((2.7996756813633845, 2.4),),
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
    xaxis_unit="ppm",
    contour_base=8.e8,
    contour_factor=1.35,
    contour_nlevels=14,
    distribution_width=2.5,
    region_separation=region_separation,
    label_peaks=False,
    xaxis_ticks=[(0, [6.65]), (1, [5.7]), (2, list([5.1 - i * 0.2 for i in range(18)])), (4, [0.65])],
    figsize=(6, 3),
    gridspec_kwargs={
        "width_ratios": [1, 6],
        "left": 0.057,
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
xlims = [xlim, (xlim[0] + 0.93 * h, xlim[0] + 0.975 * h)]
w_ratios = [x[1] - x[0] for x in xlims]
w_ratios[1] *= 10.
w_ratios = [wr / sum(w_ratios) for wr in w_ratios]
ws = [(ax1_box.x1 - ax1_box.x0) * wr for wr in w_ratios]
axs[1].remove()
axs.pop(1)
space = 0.00
lefts = [ax1_box.x0 + i * space + sum(ws[:i]) for i in range(2)]


for i, (left, width) in enumerate(zip(lefts, ws)):
    axs.insert(
        1 + i,
        fig.add_axes([left, ax1_box.y0, width - space, ax1_box.y1 - ax1_box.y0]),
    )

axs[1].spines["right"].set_visible(False)
axs[2].spines["left"].set_ls(":")
axs[2].spines["left"].set_lw(0.5)

axs[1].plot(xdata, ydata, color="k")
axs[1].set_xlim(xlims[0])
axs[1].set_xticks([])
axs[1].set_ylim(ylim)
axs[2].plot(xdata, ydata, color="k")
axs[2].set_xlim(xlims[1])
axs[2].set_xticks([])
axs[2].set_ylim(ylim)
axs[2].set_yticks([])
axs[2].text(0.05, 0.94, "$\\times 10$", fontsize=6, transform=axs[2].transAxes)
axs[1].text(0.9, 0.84, "H\\textsubscript{2}O", fontsize=6, transform=axs[1].transAxes, ha="right")
axs[1].text(0.9, 0.67, "(C, H)", fontsize=6, transform=axs[1].transAxes, ha="right")
axs[1].text(0.9, 0.41, "(B)", fontsize=6, transform=axs[1].transAxes, ha="right")


ybot, ytop = axs[0].get_ylim()
axs[0].set_ylim(0.3 * ybot, 0.223 * ytop)
axs[1].set_ylabel("$D (10^{10} \\unit{\\meter\\squared\\per\\second})$", labelpad=2)
yticks = [2e-10 + i * 5e-11 for i in range(7)]
axs[1].set_yticks(yticks)
axs[1].set_yticklabels([f"{yt * 1e10:.1f}" for yt in yticks])
axs[3].set_xlabel(axs[3].get_xlabel(), labelpad=-1)

# break_kwargs = {
#     "marker": [(-1, -3), (1, 3)],
#     "markersize": 6,
#     "linestyle": "none",
#     "color": "k",
#     "mec": "k",
#     "mew": 1,
#     "clip_on": False,
# }
# axs[1].plot([1], [1], transform=axs[1].transAxes, **break_kwargs)
# axs[1].plot([1], [0], transform=axs[1].transAxes, **break_kwargs)
# axs[2].plot([0], [1], transform=axs[2].transAxes, **break_kwargs)
# axs[2].plot([0], [0], transform=axs[2].transAxes, **break_kwargs)

label_xs = [
    6.63,  # A
    5.71,  # B
    5.048,  # C
    4.92,  # D
    4.823,  # E
    4.632,  # F
    4.4,  # G
    4.13,  # H
    4.045,  # I
    3.848,  # J
    3.25,  # K, L
    2.49,  # M, N
    2.327,  # O
    1.94,  # P
    1.87,  # Q
    1.78,  # R
    1.71,  # S
    1.65,  # T, U
    1.36,  # V
    1.22,  # W, X
    1.16,  # Y
    0.73,  # Z
    3.38,  # H2O
]
label_ys = [
    3.3e5,  # A
    2.2e5,  # B
    1.05e5,  # C
    2.9e5,  # D
    6.7e5,  # E
    6.6e5,  # F
    5.e5,  # G
    0.96e5,  # H
    5.1e5,  # I
    4.4e5,  # J
    2.9e5,  # K, L
    2.3e5,  # M, N
    2.3e5,  # O
    1.9e5,  # P
    2.4e5,  # Q
    1.9e5,  # R
    2.9e5,  # S
    4.1e5,  # T, U
    2.e5,  # V
    4.3e5,  # W, X
    6.7e5,  # Y
    6.5e5,  # Z
    0.94e5,  # H2O
]
labels = [
    "(A)",
    "(B)",
    "(C)",
    "(D)",
    "(E)",
    "(F)",
    "(G)",
    "(H)",
    "(I)",
    "(J)",
    "(K, L)",
    "(M, N)",
    "(O)",
    "(P)",
    "(Q)",
    "(R)",
    "(S)",
    "(T, U)",
    "(V)",
    "(W, X)",
    "(Y)",
    "(Z)",
    "H\\textsubscript{2}O",
]

regions, merge_regions = estimator._plot_regions(list(range(14)), "ppm")
merge_region_spans = estimator._get_3d_xaxis_spans(merge_regions, region_separation)
convert = lambda x: estimator._transform_freq_to_xaxis(
    [x], merge_regions, merge_region_spans,
)

for i, (x, y, s) in enumerate(zip(label_xs, label_ys, labels)):
    x = convert(x)[0]
    axs[0].text(
        x,
        y,
        s,
        fontsize=6,
        va="center",
        ha="center",
    )

for char, idx in zip(["a", "b", "c"], [0, 1, 3]):
    axs[idx].text(
        0.006 if idx != 1 else 0.03,
        0.91,
        f"\\textbf{{{char}.}}",
        transform=axs[idx].transAxes,
    )

fig.savefig("figures/andrographolide_dosy/andrographolide_dosy.pdf")
