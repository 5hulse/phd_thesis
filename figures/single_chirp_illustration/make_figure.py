# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Mon 05 Jun 2023 14:11:25 BST

import nmrespy as ne
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


def chirp(ax, left, width, base, height):
    x = np.linspace(left, left + width, 100000)
    x_saltire = np.linspace(0, 1, 100000)
    envelope = height * (1 - np.abs(np.cos(np.pi * x_saltire)) ** 50)
    delta_F = 30
    delta_f = -delta_F / 2
    phase = np.cos(
        np.pi * delta_F * (x_saltire - 0.5) ** 2 -
        2 * np.pi * delta_f * (x_saltire - 0.5)
    )
    y = np.abs(envelope * phase) + base
    ax.plot(x, y, color="k", lw=0.7)
    ax.plot(x, envelope + base, color="k")
    mid = (2 * left + width) / 2


fig, ax = plt.subplots(
    nrows=1,
    ncols=1,
    gridspec_kw={
        "left": 0.,
        "right": 1.,
        "bottom": 0.2,
        "top": 0.99,
        "hspace": 0.,
    },
    figsize=(5, 2.5),
)

for x in ("left", "right", "top", "bottom"):
    ax.spines[x].set_visible(False)

ax.set_xticks([])
ax.set_yticks([])

top = 0.9
tau_p = 0.37
tau_del = 0.1
t1 = 1. - tau_p - tau_del

freq_locs = [0.2, 0.6, 0.8]
sw = 20.
pts = 32000
freqs = [fl * sw for fl in freq_locs]
starts = [tau_p * fl for fl in freq_locs]
base_params = np.array([[1., 0., 0., 0.01]])

expinfo = ne.ExpInfo(dim=1, sw=sw, offset=sw / 2., default_pts=pts)

chirp(ax, 0., tau_p, 0.5, top - 0.5)

for i, (freq, start) in enumerate(zip(freqs, starts)):
    base_params[0, 2] = freq + 0.9 * ((sw / 2) - freq)
    fid = 0.23 * expinfo.make_fid(base_params) + (0.25 + 0.25j)
    tp = np.linspace(start, start + 80., pts)
    tp_acqu_idx = np.argmin(np.abs(tp - (tau_p + tau_del)))
    line = ax.plot(tp[:tp_acqu_idx + 1], fid[:tp_acqu_idx + 1].real, alpha=0.4)[0]
    color = line.get_color()
    ax.plot(tp[:tp_acqu_idx + 1], fid[:tp_acqu_idx + 1].imag, alpha=0.4, ls=":", color=color)[0]
    ax.plot(tp[tp_acqu_idx:], fid[tp_acqu_idx:].real, color=color)
    ax.plot(tp[tp_acqu_idx:], fid[tp_acqu_idx:].imag, color=color, ls=":")
    ax.plot([start, start], [0., top], color=color, lw=1.)

    arrow_y = 0.04 + i * 0.07
    ax.plot(
        [start, start],
        [arrow_y, 0.2],
        color=color,
        lw=1.,
        transform=fig.transFigure,
        clip_on=False,
        solid_capstyle="round",
    )


    ax.add_patch(
        mpl.patches.ConnectionPatch(
            xyA=[start, arrow_y],
            xyB=[tau_p + tau_del, arrow_y],
            coordsA="figure fraction",
            coordsB="figure fraction",
            lw=1.,
            arrowstyle="<->",
            color=color,
        )
    )

    midpoint = (start + tau_p + tau_del) / 2
    bbox = {"facecolor": "w", "edgecolor": "none", "pad": 1.5}
    if i == 0:
        color_txt = "red"
    elif i == 1:
        color_txt = "blue"
    elif i == 2:
        color_txt = "yellow"

    ax.text(
        midpoint,
        arrow_y,
        f"$t_0\\left(f^{{(1)}}_{{\\text{{{color_txt}}}}}\\right)$",
        bbox=bbox,
        ha="center",
        va="center",
        color=color,
        transform=fig.transFigure,
        fontsize=7,
    )

ax.set_xlim(0, 1)
ax.set_ylim(0, 1)

ax.plot([0., 1.], [0.5, 0.5], color="k", transform=ax.transAxes, clip_on=False)

ax.plot(2 * [tau_p + tau_del], [0., top], color="k")
ax.plot(2 * [tau_p + tau_del], [0.04, 0.2], color="k", transform=fig.transFigure, clip_on=False)

arrow = mpl.patches.ConnectionPatch(
    xyA=[0., (1. + top) / 2.],
    xyB=[tau_p, (1. + top) / 2.],
    coordsA="data",
    coordsB="data",
    lw=1.,
    arrowstyle="<->",
)
ax.add_patch(arrow)

arrow = mpl.patches.ConnectionPatch(
    xyA=[tau_p, (0.5 + top) / 2.],
    xyB=[tau_p + tau_del, (0.5 + top) / 2.],
    coordsA="data",
    coordsB="data",
    lw=1.,
    arrowstyle="<->",
)
ax.add_patch(arrow)

bbox = {"facecolor": "w", "edgecolor": "none", "pad": 1.5}
ax.text(tau_p / 2, (1 + top) / 2, "$\\tau_{\\text{p}}$", bbox=bbox, ha="center", va="center")
ax.text(tau_p + tau_del / 2, (0.5 + top) / 2, "$\\tau_{\\text{del}}$", bbox=bbox, ha="center", va="center")
ax.text(tau_p + tau_del + t1 / 2, (0.5 + top) / 2, "$t^{(1)}$", bbox=bbox, ha="center", va="center")

fig.savefig("figures/single_chirp_illustration/single_chirp_illustration.pdf")
