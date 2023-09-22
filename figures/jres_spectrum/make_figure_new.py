# make_figure_new.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Tue 19 Sep 2023 20:33:38 BST

import matplotlib as mpl
from matplotlib import pyplot as plt, patches, rcParams
import numpy as np
import nmrespy as ne



fig, axs = plt.subplots(
    ncols=3,
    nrows=2,
    gridspec_kw=dict(
        height_ratios=(1, 2),
        left=0.05,
        bottom=0.11,
        right=0.995,
        top=0.99,
        hspace=0.,
        wspace=0.02,
    ),
    figsize=(6., 2.5)
)

for ax in axs[0]:
    ax.set_xticks([])
    ax.set_yticks([])
for ax in axs[1, 1:]:
    ax.set_xticks([])
    ax.set_yticks([])

for ax in axs[0]:
    ax.spines["bottom"].set_ls(":")
for ax in axs[1]:
    ax.spines["top"].set_visible(False)

estimator = ne.Estimator2DJ.from_pickle("~/Documents/DPhil/results/cupid/strychnine/estimator")
f1, f2 = estimator.get_shifts(unit="ppm")
spec = estimator.spectrum_sinebell
spec_twist = spec.real
spec_abs = np.abs(spec).real

xlim = (7.31, 7.05)
ylim = (25, -25)
for ax in axs[0]:
    ax.set_xlim(xlim)
for ax in axs[1]:
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

sw1, sw2 = estimator.sw(unit="hz")
n1, n2 = estimator.default_pts
tilt_factor = (sw1 * n2) / (sw2 * n1)
spec_tilt = np.zeros_like(spec_abs)
for i, row in enumerate(spec_abs):
    spec_tilt[i] = np.roll(row, shift=int(tilt_factor * (n1 // 2 - i)))

axs[0,0].plot(f2[0], np.sum(spec_twist, axis=0), color="k")
axs[0,1].plot(f2[0], np.sum(spec_abs, axis=0), color="k")
axs[0,2].plot(f2[0], np.sum(spec_tilt, axis=0), color="k")

base, factor, nlev = 1., 2., 10
levels = [base * factor ** i for i in range(nlev)]
levels = [-x for x in reversed(levels)] + levels
colors = nlev * ["#808080"] + nlev * ["k"]
axs[1,0].contour(f2, f1, spec_twist, colors=colors, linewidths=0.25, levels=levels)
axs[1,1].contour(f2, f1, spec_abs, colors=colors, linewidths=0.25, levels=levels)
axs[1,2].contour(f2, f1, spec_tilt, colors=colors, linewidths=0.25, levels=levels)

ylim = (-100, 1200)
for ax in axs[0]:
    ax.set_ylim(ylim)

xlim = axs[1, 0].get_xlim()
ylim = axs[1, 0].get_ylim()
sfo = estimator.sfo[1]
for f in estimator.predict_multiplets(indices=[1]):
    xs = [x / sfo for x in [(f + ylim[0]), (f + ylim[1])]]
    ys = [ylim[0], ylim[1]]
    axs[1, 0].plot(xs, ys)

axs[1, 0].set_xlabel("\\textsuperscript{1}H (ppm)", labelpad=-1.2)
axs[1, 0].set_ylabel("Hz")

xs = [7.277, 7.21, 7.134, 7.08]
ys = [50, 10, 210, 50]
for x, y in zip(xs, ys):
    axs[0, 2].text(x, y, "*")

for i, ax in enumerate(axs[0]):
    ax.text(0.01, 0.83, f"\\textbf{{{chr(97 + i)}.}}", transform=ax.transAxes)

xys = [
    (0.205, 0.705),
    (0.54, 0.12),
    (0.82, 0.86),
    (0.825, 0.295),
]
heights = [0.42, 0.4, 0.4, 0.42]
widths = [0.1, 0.1, 0.1, 0.1]
angles = [-20, -20, -20, -20]
for xy, width, height, angle in zip(xys, widths, heights, angles):
    ellipse = patches.Ellipse(
        xy=xy, width=width, height=height, angle=angle,
        facecolor="none", edgecolor="k", ls=":", transform=axs[1, 0].transAxes,
        lw=0.6,
    )
    axs[1, 0].add_artist(ellipse)

fig.savefig("figures/jres_spectrum/jres_spectrum_new.pdf")
