# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Wed 19 Jul 2023 12:31:32 BST

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

for ax in axs.flatten():
    ax.set_xticks([])
    ax.set_yticks([])

for ax in axs[0]:
    ax.spines["bottom"].set_ls(":")
for ax in axs[1]:
    ax.spines["top"].set_visible(False)

red = mpl.rcParams["axes.prop_cycle"].by_key()["color"][0]

# 2DJ spectrum
expinfo = ne.ExpInfo(dim=2, sw=(1., 16.), default_pts=(64, 512))
params = np.array([
    [1, 0, 0.32, 6.32, 0.3, 0.7],
    [1, 0, 0.08, 6.08, 0.3, 0.7],
    [1, 0, -0.08, 5.92, 0.3, 0.7],
    [1, 0, -0.32, 5.68, 0.3, 0.7],
    [1, 0, 0.27, 0.27, 0.3, 0.7],
    [1, 0, 0.03, 0.03, 0.3, 0.7],
    [1, 0, -0.03, -0.03, 0.3, 0.7],
    [1, 0, -0.27, -0.27, 0.3, 0.7],
    [1, 0, 0.35, -5.65, 0.3, 0.7],
    [1, 0, 0.05, -5.95, 0.3, 0.7],
    [1, 0, -0.05, -6.05, 0.3, 0.7],
    [1, 0, -0.35, -6.35, 0.3, 0.7],
])
fid = expinfo.make_fid(params)
spec = ne.sig.ft(ne.sig.sinebell_apodisation(fid))
spec_twist = spec.real
spec_abs = np.abs(spec).real

sw1, sw2 = expinfo.sw(unit="hz")
n1, n2 = expinfo.default_pts
tilt_factor = (sw1 * n2) / (sw2 * n1)
spec_tilt = np.zeros_like(spec_abs)
for i, row in enumerate(spec_abs):
    spec_tilt[i] = np.roll(row, shift=int(tilt_factor * (n1 // 2 - i)))

f1_shifts, f2_shifts = expinfo.get_shifts()
base, factor, nlev = 0.1, 1.8, 5
levels = [base * factor ** i for i in range(nlev)]
levels = [-x for x in reversed(levels)] + levels

axs[0,0].plot(f2_shifts[0], np.sum(spec_twist, axis=0), color="k")
axs[0,1].plot(f2_shifts[0], np.sum(spec_abs, axis=0), color="k")
axs[0,2].plot(f2_shifts[0], np.sum(spec_tilt, axis=0), color="k")

axs[1,0].contour(f2_shifts, f1_shifts, spec_twist, colors="k", linewidths=0.2, levels=levels)
axs[1,1].contour(f2_shifts, f1_shifts, spec_abs, colors="k", linewidths=0.2, levels=levels)
axs[1,2].contour(f2_shifts, f1_shifts, spec_tilt, colors="k", linewidths=0.2, levels=levels)

axs[1,1].set_xlim(f2_shifts[0][0], f2_shifts[0][-1])
axs[1,1].set_ylim(f1_shifts[0][0], f1_shifts[-1][0])
axs[0,1].set_xlim(f2_shifts[0][0], f2_shifts[0][-1])
ylim = axs[1,1].get_ylim()
axs[1,1].plot([6 + ylim[0], 6 + ylim[1]], [ylim[0], ylim[1]], color=red)
arc = patches.Arc(
    xy=(6 + ylim[0], f1_shifts[0][0]),
    width=2,
    height=0.16,
    angle=180,
    theta1=0,
    theta2=52.5,
    edgecolor=red,
    lw=1,
)
axs[1,1].add_patch(arc)

axs[1,0].annotate(
    text="", xy=(0, -0.04), xytext=(1, -0.04), xycoords="axes fraction", clip_on=False,
    arrowprops={"arrowstyle": "<|-|>", "shrinkA": 0, "shrinkB": 0, "mutation_scale": 6, "lw": 0.6, "facecolor": "k"},
)
for ax in axs[1]:
    ax.scatter(0.5, 0, marker=2, color="k", transform=ax.transAxes, clip_on=False)
    ax.scatter(0, 0.5, marker=1, color="k", transform=ax.transAxes, clip_on=False)

axs[1, 0].text(0.5, -0.15, "$f_{\\mathrm{sw}}^{(2)}$", transform=axs[1,0].transAxes, ha="center")
axs[1,0].text(0.52, 0.05, "$f_{\\mathrm{off}}$", transform=axs[1,0].transAxes)
axs[1,0].annotate(
    text="", xy=(-0.03, 0), xytext=(-0.03, 1), xycoords="axes fraction", clip_on=False,
    arrowprops={"arrowstyle": "<|-|>", "shrinkA": 0, "shrinkB": 0, "mutation_scale": 6, "lw": 0.6, "facecolor": "k"},
)
axs[1,0].text(-0.15, 0.5, "$f_{\\mathrm{sw}}^{(1)}$", transform=axs[1,0].transAxes, va="center")
axs[1,0].text(0.02, 0.52, "$0$", transform=axs[1,0].transAxes)



axs[1,2].set_xlim(f2_shifts[0][0], f2_shifts[0][-1])
axs[1,2].set_ylim(f1_shifts[0][0], f1_shifts[-1][0])
axs[0,2].set_xlim(f2_shifts[0][0], f2_shifts[0][-1])


axs[1,1].text(0.15, 0.05, "\\ang{45}", color=red, transform=axs[1,1].transAxes)

ylim = (-2, 22)
scale = 2.1
axs[0,0].set_xlim(f2_shifts[0][0], f2_shifts[0][-1])
axs[0, 0].set_ylim(ylim)
axs[0, 1].set_ylim(ylim)
axs[0, 2].set_ylim([scale * y for y in ylim])
axs[0, 1].text(0.99, 0.95, f"\\times {scale}", ha="right", va="top", transform=axs[0, 1].transAxes, fontsize=7)

for i, ax in enumerate(axs[0]):
    ax.text(0.01, 0.83, f"\\textbf{{{chr(97 + i)}.}}", transform=ax.transAxes)

fig.savefig("figures/jres_spectrum/jres_spectrum.pdf")
