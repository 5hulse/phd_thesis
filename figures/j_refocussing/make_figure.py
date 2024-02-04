# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Sun 04 Feb 2024 15:16:10 EST

from dataclasses import dataclass
import matplotlib as mpl
from matplotlib.patches import FancyArrowPatch, Rectangle
import matplotlib.pyplot as plt
import numpy as np


fig, ax = plt.subplots(
    nrows=1,
    ncols=1,
    gridspec_kw={
        "left": 0.05,
        "right": 0.98,
        "bottom": 0.,
        "top": 0.94,
        "hspace": 0.,
    },
    figsize=(6, 1.5),
)

for x in ("top", "bottom", "left", "right"):
    ax.spines[x].set_visible(False)
ax.set_xticks([])
ax.set_yticks([])
ax.set_ylim(0, 1)
ax.axhline(0.5, color="k")
ax.text(
    -0.01,
    0.52,
    "\\textsuperscript{1}H",
    transform=ax.transAxes,
    clip_on=False,
    va="center",
    ha="right",
)

def rectangular(ax, left, width, height):
    right = left + width
    xs = [left, left, right, right]
    ys = [0.5, 0.5 + height, 0.5 + height, 0.5]
    ax.fill(xs, ys, "k", transform=ax.transAxes)


def acquisition(ax, left, one_over_sw, extra_acqu, height):
    right = left + one_over_sw + extra_acqu
    n_points = 512
    xs = np.linspace(left, right, n_points)
    xs_fid = np.linspace(0, 30, n_points)
    ys = np.cos(xs_fid) * np.exp(0.03 * -xs_fid)
    scale = (height) / ys[0]
    ys = ys * scale + 0.5
    frac = one_over_sw / (one_over_sw + extra_acqu)
    split = int(frac * n_points)
    ax.plot(xs[:split], ys[:split], color="k", transform=ax.transAxes)
    ax.plot(xs[split:], ys[split:], color="k", ls=":", transform=ax.transAxes)

    arrow_width = one_over_sw
    arrow_height = 0.05
    arrow = FancyArrowPatch(
        posA=(left, arrow_height),
        posB=(left + arrow_width, arrow_height),
        arrowstyle="<->",
        shrinkA=0,
        shrinkB=0,
        mutation_scale=10,
        transform=ax.transAxes,
        clip_on=False,
        lw=0.7,
    )
    ax.add_patch(arrow)

    ax.text(
        left + 0.5 * arrow_width,
        arrow_height + 0.07,
        "$\\nicefrac{1}{f_{\\text{sw}}^{(1)}}$",
        ha="center",
        transform=ax.transAxes,
        clip_on=False,
    )
    ax.text(
        left + 0.5 * arrow_width,
        0.85,
        "$t^{(2)}$",
        ha="center",
        transform=ax.transAxes,
        clip_on=False,
    )


class Widths:
    def __init__(
        self,
        left_pad,
        ninty,
        t_one_over_two,
        one_over_sw,
        j_refocus,
        extra_acqu,
        right_pad,
    ):
        total_width = left_pad + ninty + 2 * t_one_over_two + 1.5 * one_over_sw + j_refocus + extra_acqu + right_pad
        self.left_pad = left_pad / total_width
        self.ninty = ninty / total_width
        self.t_one_over_two = t_one_over_two / total_width
        self.one_over_sw = one_over_sw / total_width
        self.j_refocus = j_refocus / total_width
        self.extra_acqu = extra_acqu / total_width
        self.right_pad = right_pad / total_width
        self.pre_post_one_over_twoeighty_delay = 0.5 * (0.5 * self.one_over_sw - 2 * self.ninty)
        assert self.pre_post_one_over_twoeighty_delay > 0.0


widths = Widths(2.0, 0.5, 12.0, 30.0, 22.0, 10.0, 0.0)
pulse_height = 0.43

left = widths.left_pad

# 90
rectangular(
    ax,
    left,
    widths.ninty,
    pulse_height,
)
ax.text(
    left + widths.ninty / 2,
    0.99,
    "$\\nicefrac{\\pi}{2}$",
    ha="center",
    transform=ax.transAxes,
)

left += widths.ninty + widths.t_one_over_two / 2

ax.text(
    left,
    0.75,
    "$\\nicefrac{t^{(1)}}{2}$",
    ha="center",
    va="center",
    transform=ax.transAxes,
)

left += widths.t_one_over_two / 2

ax.plot(
    [left, left],
    [0.5 - pulse_height, 0.5 + pulse_height],
    ls=":",
    color="k",
    lw=0.6,
    transform=ax.transAxes,
    clip_on=False,
)

arrow_width = 0.5 * widths.one_over_sw
arrow_height = 0.3
arrow = FancyArrowPatch(
    posA=(left, arrow_height),
    posB=(left + arrow_width, arrow_height),
    arrowstyle="<->",
    shrinkA=0,
    shrinkB=0,
    mutation_scale=10,
    transform=ax.transAxes,
    clip_on=False,
    lw=0.7,
)
ax.add_patch(arrow)

ax.text(
    left + 0.5 * arrow_width,
    arrow_height + 0.07,
    "$\\nicefrac{1}{2 f_{\\text{sw}}^{(1)}}$",
    ha="center",
    transform=ax.transAxes,
    clip_on=False,
)

left += widths.pre_post_one_over_twoeighty_delay

rectangular(
    ax,
    left,
    2 * widths.ninty,
    pulse_height,
)
ax.text(
    left + widths.ninty,
    0.99,
    "$\\pi$",
    ha="center",
    transform=ax.transAxes,
)

left += 2 * widths.ninty + widths.pre_post_one_over_twoeighty_delay

j_block = Rectangle(
    (left, 0.5),
    widths.j_refocus,
    pulse_height,
    facecolor="#e0e0e0",
    edgecolor="none",
    transform=ax.transAxes,
    zorder=-1
)
ax.add_patch(j_block)
left += widths.j_refocus / 2.

ax.text(
    left,
    0.7,
    "J-refocussing element",
    va="center",
    ha="center",
    zorder=2000,
    clip_on=False,
    transform=ax.transAxes,
    fontsize=9,
)
left += widths.j_refocus / 2.
left += widths.t_one_over_two / 2

ax.text(
    left,
    0.75,
    "$\\nicefrac{t^{(1)}}{2}$",
    ha="center",
    va="center",
    transform=ax.transAxes,
)

left += widths.t_one_over_two / 2

acquisition(
    ax,
    left,
    widths.one_over_sw,
    widths.extra_acqu,
    0.37,
)

fig.savefig("figures/j_refocussing/j_refocussing.pdf")
