# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Wed 26 Apr 2023 22:20:27 BST

import colorsys
import matplotlib as mpl
from matplotlib import colors as mc, patches, pyplot as plt
import numpy as np
import nmrespy as ne


def lighten_color(color, amount=0.5):
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = np.array(colorsys.rgb_to_hls(*mc.to_rgb(c)))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])


estimator_path = "~/Documents/DPhil/results/cupid/sucrose/estimator"
COLORS = mpl.rcParams["axes.prop_cycle"].by_key()["color"]
red = COLORS[0]
red_pale = lighten_color(red, 0.2)
blue = COLORS[1]
blue_pale = lighten_color(blue, 0.2)

t1_time = 1.5
t1_pts = 5
t2_time = 2.
t2_pts = 20
marker_size = 20
xpad = 0.15
ypad = 0.15

ax_h_ratios = 2.

t1 = np.linspace(0, t1_time, t1_pts)
t2 = np.linspace(0, t2_time, t2_pts)
pos_t1_x, pos_t1_y = np.meshgrid(t2, t1, indexing="ij")
pos_t1_x, pos_t1_y = pos_t1_x.flatten(), pos_t1_y.flatten()
neg_t1_x, neg_t1_y = np.meshgrid(t2, -t2, indexing="ij")
neg_t1_x, neg_t1_y = neg_t1_x.flatten(), neg_t1_y.flatten()

fig, axs = plt.subplots(
    nrows=1,
    ncols=2,
    gridspec_kw=dict(
        left=0.02,
        right=0.98,
        bottom=0.02,
        top=0.98,
        wspace=0.02,
    ),
    figsize=(4.5, 4.5),
)

for ax in axs:
    for orientation in ("top", "bottom", "left", "right"):
        ax.spines[orientation].set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])

axs[0].set_xlim(-0.1, t2_time + xpad + 0.1)
axs[0].set_ylim(-t2_time - ypad - 0.04, t2_time + ypad + 0.04)
neg_scatter = axs[0].scatter(neg_t1_x, neg_t1_y, color=blue_pale, edgecolor="none", s=marker_size)  # noqa: E501
axs[0].scatter(pos_t1_x, pos_t1_y, color=red_pale, edgecolor="none", s=marker_size)
axs[0].scatter(t2, -t2, color=blue, edgecolor="none", s=marker_size)
axs[0].scatter(pos_t1_x[0::t1_pts], pos_t1_y[0::t1_pts], color=red, s=marker_size, edgecolor="none")  # noqa: E501
axs[0].plot([0, 0], [-t2_time - ypad, t2_time + ypad], color="k", zorder=0)
axs[0].plot([-xpad, t2_time + xpad], [0, 0], color="k", zorder=0)
axs[0].scatter([0], [t2_time + ypad], marker="^", color="k")
axs[0].scatter([0], [-t2_time - ypad], marker="v", color="k")
axs[0].scatter([t2_time + xpad], [0], marker=">", color="k")

axs[0].text(
    x=0.07,
    y=t2_time + ypad + 0.01,
    s="$t^{(1)}$",
    color="k",
    verticalalignment="center",
)
axs[0].text(
    x=0.07,
    y=-t2_time - ypad + 0.01,
    s="$-t^{(1)}$",
    color="k",
    verticalalignment="center",
)
axs[0].text(
    x=t2_time + xpad,
    y=0.07,
    s="$t^{(2)}$",
    color="k",
    horizontalalignment="center",
)

axs[0].text(
    x=0.1,
    y=0.15,
    s="$\\symbf{y}_{\\text{1D}}(t) = \\symbf{Y}_{\\text{2DJ}}\\left(0, t\\right)$",
    color=red,
    bbox={"facecolor": "w", "edgecolor": red, "boxstyle": "round"},
)
axs[0].text(x=0.7, y=0.67, s="\\textbf{FT}", color=red, transform=axs[0].transAxes)  # noqa: E501
axs[0].text(x=0.7, y=0.345, s="\\textbf{FT}", color=blue, transform=axs[0].transAxes)  # noqa: E501

axs[0].text(
    x=0.35,
    y=-1.9,
    s="$\\symbf{y}_{\\ang{-45}}(t) = \\symbf{Y}_{\\text{2DJ}}\\left(-t, t\\right)$",
    color=blue,
    bbox={"facecolor": "w", "edgecolor": blue, "boxstyle": "round"},
)
arrow = mpl.patches.FancyArrowPatch(
    (1.365, -1.4),
    (2.45, -0.5),
    arrowstyle="-|>",
    connectionstyle=mpl.patches.ConnectionStyle.Arc3(rad=-0.4),
    color=blue,
    mutation_scale=10,
    clip_on=False,
)
axs[0].add_patch(arrow)

arrow = mpl.patches.FancyArrowPatch(
    (1.365, 0),
    (2.45, 0.9),
    arrowstyle="-|>",
    connectionstyle=mpl.patches.ConnectionStyle.Arc3(rad=-0.4),
    color=red,
    mutation_scale=10,
    clip_on=False,
)
axs[0].add_patch(arrow)

diameter = 5 * np.sqrt(2) * t2[1]
arc = patches.Arc(
    xy=(0., 0.),
    width=diameter,
    height=diameter,
    angle=315.,
    theta1=0.,
    theta2=45.,
)
axs[0].add_patch(arc)
axs[0].text(
    x=0.37,
    y=-0.2,
    s="\\textbf{$-$45°}",
)

pts = (64, 2 ** 17)
estimator = ne.Estimator2DJ.from_pickle(estimator_path)
estimator._pts = pts
estimator._default_pts = pts
slice_ = slice(*estimator.convert((None, (4.7, 3.97)), "ppm->idx")[1])
shifts = estimator.get_shifts(unit="ppm", meshgrid=False)[1][slice_]
cupid_spectrum = estimator.cupid_spectrum().real[slice_]
cupid_max_height = np.max(cupid_spectrum)
coupled_fid = estimator.make_fid_from_result()[0]
coupled_fid[0] *= 0.5
coupled_spectrum = ne.sig.ft(coupled_fid).real[slice_]

axs[1].plot(shifts, coupled_spectrum + 1.1 * cupid_max_height, color=red, lw=1.)
axs[1].plot(shifts, cupid_spectrum, color=blue, lw=1.)
axs[1].set_xlim(reversed(axs[1].get_xlim()))

fig.savefig("figures/neg_45_signal/neg_45_signal.pdf")
