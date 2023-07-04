# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Mon 03 Jul 2023 16:08:24 BST


import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import nmrespy as ne

COLORS = mpl.rcParams["axes.prop_cycle"].by_key()["color"]

estimator = ne.EstimatorDiffBi.from_pickle(
    "/home/simon/Documents/DPhil/results/diffusion/gluc_val_thre/estimator_2",
)
estimator.sigma = 1.

# ==============
region_separation = 0.01
y_ranges = [(4.e-10, 7.5e-10), (18.e-10, 19.5e-10)]
ratio_1d_2d = [1, 3]
width_ratios = [1, 5]
bottom, top, left, right = 0.08, 0.99, 0.06, 0.99
# ==============

fig, axs = estimator.plot_dosy(
    y_range=(4.5e-10, 7.e-10),
    y_pts=250,
    xaxis_unit="ppm",
    contour_base=6.e8,
    contour_factor=1.7,
    contour_nlevels=12,
    distribution_width=1.5,
    region_separation=region_separation,
    label_peaks=False,
    xaxis_ticks=[
        (0, [5.15]),
        (1, [4.7]),
        (2, [4.56]),
        (4, [3.8, 3.75, 3.7, 3.65, 3.6]),
        (5, [3.5, 3.45, 3.4, 3.35, 3.3]),
        (6, [3.16]),
        (7, [2.25, 2.2, 2.15]),
        (8, [1.24]),
    ],
    figsize=(6, 3.5),
    gridspec_kwargs={
        "width_ratios": [1, 6],
        "height_ratios": [1, 2],
        "left": 0.057,
        "right": 0.995,
        "top": 0.982,
        "bottom": 0.082,
    },
    spectrum_line_kwargs={"linewidth": .6},
    oscillator_line_kwargs={"linewidth": .0},
    contour_kwargs={"linewidths": 0.1},
)

axs[0].set_ylim(-1e5, 2.e6)

patch_ys = [
    (6.3e-10, 6.5e-10),
    (6.1e-10, 6.3e-10),
    (5.e-10, 5.75e-10),
]
for i, y in enumerate(patch_ys):
    axs[1].axhspan(y[0], y[1], facecolor=COLORS[i], zorder=-1, alpha=0.1)
    axs[2].axhspan(y[0], y[1], facecolor=COLORS[i], zorder=-1, alpha=0.1)

label_xs = [0.05, 0.355, 0.485]
label_ys = [0.77, 0.63, 0.25]
label_ss = ["threonine", "valine", "glucose"]
for i, (x, y, s) in enumerate(zip(label_xs, label_ys, label_ss)):
    axs[1].text(x, y, s, color=COLORS[i], transform=axs[1].transAxes, fontsize=8)

axs[0].text(0.072, 0.86, "H\\textsubscript{2}O", transform=axs[0].transAxes, fontsize=7)
axs[1].set_xlim(-7.1e13, 4.e12)

xs = [0.195, 0.06, 0.195]
ys = [0.95, 0.645, 0.645]
for i, (x, y) in enumerate(zip(xs, ys)):
    fig.text(x, y, f"\\textbf{{{chr(97 + i)}.}}")


ticks = np.linspace(4.5e-10, 7e-10, 6)
axs[1].set_yticks(ticks)
axs[1].set_yticklabels([str(t * 1e10) for t in ticks])
axs[1].set_ylabel("$D (10^{10} \\unit{\\meter\\squared\\per\\second})$", labelpad=2)

label_xs = [0.84, 0.575, 0.52, 0.235, 0.164, 0.102]
label_ys = [1.7e6, 1.7e6, 1.7e6, 5e5, 1.7e6, 1.7e6]
label_ss = ["(A)", "(A)", "(B)", "(B)", "(C)", "(C)"]
label_colors = [0, 1, 0, 1, 0, 1]

for x, y, s, col in zip(label_xs, label_ys, label_ss, label_colors):
    axs[0].text(x, y, s, color=COLORS[col], fontsize=7)

fig.savefig("figures/glu_thre_val_diffusion/glu_thre_val_diffusion.pdf")
