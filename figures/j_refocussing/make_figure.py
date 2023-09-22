# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Tue 19 Sep 2023 14:37:41 BST

from dataclasses import dataclass
import matplotlib as mpl
from matplotlib.patches import FancyArrowPatch, Rectangle
import matplotlib.pyplot as plt
import numpy as np


fig, axs = plt.subplots(
    nrows=2,
    ncols=1,
    gridspec_kw={
        "left": 0.05,
        "right": 0.98,
        "bottom": 0.,
        "top": 0.96,
        "hspace": 0.,
    },
    figsize=(6, 2),
)

for i, ax in enumerate(axs):
    for x in ("top", "bottom", "left", "right"):
        ax.spines[x].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_ylim(0, 1)
    ax.axhline(0.5, color="k")
    if i == 0:
        ax.text(
            -0.01,
            0.52,
            "\\textsuperscript{1}H",
            transform=ax.transAxes,
            clip_on=False,
            va="center",
            ha="right",
        )
    else:
        ax.text(
            -0.01,
            0.5,
            "$g_z$",
            transform=ax.transAxes,
            clip_on=False,
            va="center",
            ha="right",
        )


def normalise_widths(*widths):
    total_w = sum(widths)
    return [w / total_w for w in widths]


def rectangular(ax, left, width, height):
    right = left + width
    xs = [left, left, right, right]
    ys = [0.5, 0.5 + height, 0.5 + height, 0.5]
    ax.fill(xs, ys, "k", transform=ax.transAxes)


def sine(ax, left, width, height, neg=False):
    right = left + width
    # mean = (left + right) / 2
    xs = np.linspace(left, right, 100)
    ys = np.sin(np.linspace(0, np.pi, 100))
    scale = height / np.amax(ys)
    if neg:
        ys *= -1
    ys = ys * scale + 0.5
    ax.fill(xs, ys, "k", transform=ax.transAxes)


def acquisition(ax, left, width, height):
    right = left + width
    xs = np.linspace(left, right, 1000)
    xs_fid = np.linspace(0, 30, 1000)
    ys = np.cos(xs_fid) * np.exp(0.03 * -xs_fid)
    scale = (height) / ys[0]
    ys = ys * scale + 0.5
    split = 500
    ax.plot(xs[:split], ys[:split], color="k", transform=ax.transAxes)
    ax.plot(xs[split:], ys[split:], color="k", ls=":", transform=ax.transAxes)

    arrow_width = width * 0.5
    arrow_height = 0.9
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

    axs[0].text(
        left + 0.5 * arrow_width,
        0.97,
        "$\\nicefrac{1}{f_{\\text{sw}}^{(1)}}$",
        ha="center",
        transform=ax.transAxes,
        clip_on=False,
    )


class Widths:
    ninty = 0.01
    acqu = 0.2
    left_pad = 0.02
    right_pad = 0.

    def __init__(self, n_ninty, **kwargs):
        names = []
        widths = []
        for n, w in kwargs.items():
            names.append(n)
            widths.append(w)
        widths = normalise_widths(*widths)
        extra_scale = 1 - (n_ninty * self.ninty + self.acqu + self.left_pad + self.right_pad)
        widths = [w * extra_scale for w in widths]

        for n, w in zip(names, widths):
            setattr(self, n, w)


class JWidths(Widths):

    def __init__(self):
        super().__init__(n_ninty=3, t_one=0.1, grad=0.03, j_elem=0.1, gap=0.05)
        self.t_one /= 2.
        self.grad /= 3.
        self.gap /= 2.


widths = JWidths()

left = widths.left_pad

# 90
rectangular(
    axs[0],
    left,
    widths.ninty,
    0.43,
)
axs[0].text(
    left + widths.ninty / 2,
    0.99,
    "$\\nicefrac{\\pi}{2}$",
    ha="center",
    transform=axs[0].transAxes,
)

left += widths.ninty + widths.t_one / 2

axs[0].text(
    left,
    0.75,
    "$\\nicefrac{t^{(1)}}{2}$",
    ha="center",
    va="center",
    transform=axs[0].transAxes,
)

left += widths.t_one / 2

sine(
    axs[1],
    left,
    widths.grad,
    0.2,
)

arrow_width = 2 * widths.grad + 2 * widths.gap + 2 * widths.ninty
arrow_height = 1.
arrow = FancyArrowPatch(
    posA=(left, arrow_height),
    posB=(left + arrow_width, arrow_height),
    arrowstyle="<->",
    shrinkA=0,
    shrinkB=0,
    mutation_scale=10,
    transform=axs[1].transAxes,
    clip_on=False,
    lw=0.7,
)
axs[1].add_patch(arrow)

axs[1].text(
    left + 0.5 * arrow_width,
    1.07,
    "$\\nicefrac{1}{2 f_{\\text{sw}}^{(1)}}$",
    ha="center",
    transform=axs[1].transAxes,
    clip_on=False,
)

left += widths.grad + widths.gap


rectangular(
    axs[0],
    left,
    2 * widths.ninty,
    0.43,
)
axs[0].text(
    left + widths.ninty,
    0.99,
    "$\\pi$",
    ha="center",
    transform=axs[0].transAxes,
)

left += 2 * widths.ninty + widths.gap

sine(
    axs[1],
    left,
    widths.grad,
    -0.2,
)

left += widths.grad

j_block = Rectangle(
    (left, 0.),
    widths.j_elem,
    0.9,
    facecolor="#e0e0e0",
    edgecolor="none",
    transform=axs[0].transAxes,
    zorder=1000
)
axs[0].add_patch(j_block)
j_block = Rectangle(
    (left, 0.1),
    widths.j_elem,
    0.9,
    facecolor="#e0e0e0",
    edgecolor="none",
    transform=axs[1].transAxes,
    zorder=1000
)
axs[1].add_patch(j_block)
left += widths.j_elem / 2.

axs[1].text(
    left,
    1,
    "J-refocussing element",
    va="center",
    ha="center",
    zorder=2000,
    clip_on=False,
    transform=axs[1].transAxes,
    fontsize=9,
)
left += widths.j_elem / 2.
sine(
    axs[1],
    left,
    widths.grad,
    -0.4,
)

left += widths.grad + widths.t_one / 2

axs[0].text(
    left,
    0.75,
    "$\\nicefrac{t^{(1)}}{2}$",
    ha="center",
    va="center",
    transform=axs[0].transAxes,
)

left += widths.t_one / 2

acquisition(
    axs[0],
    left,
    widths.acqu,
    0.37,
)

# # 180
# rectangular(
#     axs[0],
#     left,
#     2 * widths.ninty,
#     0.4,
# )
# axs[0].text(
#     left + widths.ninty,
#     0.99,
#     "$\\pi$",
#     ha="center",
#     transform=axs[0].transAxes,
# )

# left += (2 * widths.ninty) + (widths.t_one / 2)

# axs[0].text(
#     left,
#     0.75,
#     "$\\nicefrac{t^{(1)}}{2}$",
#     ha="center",
#     va="center",
#     transform=axs[0].transAxes,
# )

# left += widths.t_one / 2

# # acquistion
# acquisition(
#     axs[0],
#     left,
#     widths.acqu,
#     0.4,
# )

# axs[0].text(
#     left + 0.66 * widths.acqu,
#     0.75,
#     "$t^{(2)}$",
#     ha="center",
#     va="center",
#     transform=axs[0].transAxes,
# )


# # Zangger-Sterk
# # =============
# class ZS_Widths(Widths):

#     def __init__(
#         self,
#         selective: float,
#         t_one: float,
#         grad: float,
#         one_over_sw: float,
#         middle_gap: float,
#     ):
#         super().__init__(
#             n_ninty=3,
#             selective=selective,
#             t_one=t_one,
#             grad=grad,
#             one_over_sw=one_over_sw,
#             middle_gap=middle_gap,
#         )
#         self.grad /= 3
#         self.t_one /= 2
#         self.one_over_sw /= 2
#         self.middle_gap /= 2


# zs_widths = ZS_Widths(2., 20., 1., 4., 1.)

# left = zs_widths.left_pad

# # 90
# rectangular(
#     axs[1],
#     left,
#     zs_widths.ninty,
#     0.4,
# )
# axs[1].text(
#     left + zs_widths.ninty / 2,
#     0.99,
#     "$\\nicefrac{\\pi}{2}$",
#     ha="center",
#     transform=axs[1].transAxes,
# )
# left += zs_widths.ninty + (zs_widths.t_one / 2)

# axs[1].text(
#     left,
#     0.75,
#     "$\\nicefrac{t^{(1)}}{2}$",
#     ha="center",
#     va="center",
#     transform=axs[1].transAxes,
# )

# left += zs_widths.t_one / 2

# axs[1].plot([left, left], [0.5, 1], ls=":", color="k", transform=axs[1].transAxes)

# # g1
# sine(
#     axs[2],
#     left,
#     zs_widths.grad,
#     0.2,
# )
# axs[2].text(
#     left + zs_widths.grad / 2,
#     0.8,
#     "$g_1$",
#     ha="center",
#     transform=axs[2].transAxes,
# )

# left += zs_widths.grad

# arrow_width = 2 * zs_widths.one_over_sw + 2 * zs_widths.ninty
# arrow_height = 1.1
# arrow = FancyArrowPatch(
#     posA=(left, arrow_height),
#     posB=(left + arrow_width, arrow_height),
#     arrowstyle="<->",
#     shrinkA=0,
#     shrinkB=0,
#     mutation_scale=10,
#     transform=axs[2].transAxes,
#     clip_on=False,
#     lw=0.7,
# )
# axs[2].add_patch(arrow)

# axs[2].text(
#     left + arrow_width / 2,
#     arrow_height + 0.12,
#     "$\\nicefrac{1}{2f_{\\text{sw}}^{(2)}}$",
#     ha="center",
#     va="center",
#     transform=axs[2].transAxes,
# )

# left += zs_widths.one_over_sw
# # 180
# rectangular(
#     axs[1],
#     left,
#     2 * zs_widths.ninty,
#     0.4,
# )
# axs[1].text(
#     left + zs_widths.ninty,
#     0.99,
#     "$\\pi$",
#     ha="center",
#     transform=axs[1].transAxes,
# )

# left += 2 * zs_widths.ninty + zs_widths.one_over_sw


# # g2
# sine(
#     axs[2],
#     left,
#     zs_widths.grad,
#     0.2,
#     neg=True,
# )
# axs[2].text(
#     left + zs_widths.grad / 2,
#     0.15,
#     "$g_2$",
#     ha="center",
#     va="bottom",
#     transform=axs[2].transAxes,
# )

# left += zs_widths.grad + zs_widths.middle_gap

# # 180 selective
# gaussian(
#     axs[1],
#     left,
#     zs_widths.selective,
#     0.4,
# )
# axs[1].text(
#     left + zs_widths.selective / 2,
#     0.99,
#     "$π_{\\text{sel}}$",
#     ha="center",
#     transform=axs[1].transAxes,
# )

# # g_sel
# rectangular(
#     axs[2],
#     left,
#     zs_widths.selective,
#     0.1,
# )
# axs[2].text(
#     left + zs_widths.selective / 2,
#     0.7,
#     "$g_{\\text{sel}}$",
#     ha="center",
#     transform=axs[2].transAxes,
# )

# left += zs_widths.selective + zs_widths.middle_gap

# # g3
# sine(
#     axs[2],
#     left,
#     zs_widths.grad,
#     0.35,
#     neg=True,
# )
# axs[2].text(
#     left + zs_widths.grad / 2,
#     0.,
#     "$g_3$",
#     ha="center",
#     va="bottom",
#     transform=axs[2].transAxes,
# )

# left += zs_widths.grad

# axs[1].plot([left, left], [0.5, 1], ls=":", color="k", transform=axs[1].transAxes)

# left += zs_widths.t_one / 2

# axs[1].text(
#     left,
#     0.75,
#     "$\\nicefrac{t^{(1)}}{2}$",
#     ha="center",
#     va="center",
#     transform=axs[1].transAxes,
# )

# left += zs_widths.t_one / 2

# # acquistion
# acquisition(
#     axs[1],
#     left,
#     zs_widths.acqu,
#     0.4,
# )

# axs[1].text(
#     left + 0.66 * zs_widths.acqu,
#     0.75,
#     "$t^{(2)}$",
#     ha="center",
#     va="center",
#     transform=axs[1].transAxes,
# )

# arrow_width = zs_widths.acqu
# arrow_height = 1.1
# arrow = FancyArrowPatch(
#     posA=(left, arrow_height),
#     posB=(left + arrow_width, arrow_height),
#     arrowstyle="<->",
#     shrinkA=0,
#     shrinkB=0,
#     mutation_scale=10,
#     transform=axs[2].transAxes,
#     clip_on=False,
#     lw=0.7,
# )
# axs[2].add_patch(arrow)

# axs[2].text(
#     left + arrow_width / 2,
#     arrow_height + 0.12,
#     "$\\nicefrac{1}{f_{\\text{sw}}^{(2)}}$",
#     ha="center",
#     va="center",
#     transform=axs[2].transAxes,
# )

# # BIRD
# # =============
# class BIRD_Widths(Widths):

#     def __init__(
#         self,
#         t_one: float,
#         grad: float,
#         one_over_2J: float,
#         one_over_sw: float,
#     ):
#         super().__init__(
#             n_ninty=9,
#             t_one=t_one,
#             grad=grad,
#             one_over_2J=one_over_2J,
#             one_over_sw=one_over_sw,
#         )
#         self.grad /= 3
#         self.t_one /= 2
#         self.one_over_2J /= 4
#         self.one_over_sw /= 2


# bird_widths = BIRD_Widths(12., 1., 8., 8.)

# left = bird_widths.left_pad

# # 90
# rectangular(
#     axs[3],
#     left,
#     bird_widths.ninty,
#     0.4,
# )
# axs[3].text(
#     left + bird_widths.ninty / 2,
#     0.99,
#     "$\\nicefrac{\\pi}{2}$",
#     ha="center",
#     transform=axs[3].transAxes,
# )

# left += bird_widths.ninty

# axs[3].text(
#     left + bird_widths.one_over_2J / 2,
#     0.75,
#     "$\\nicefrac{1}{2J_{\\text{CH}}}$",
#     ha="center",
#     va="center",
#     transform=axs[3].transAxes,
#     fontsize=7,
# )
# left += bird_widths.one_over_2J

# rectangular(
#     axs[3],
#     left,
#     2 * bird_widths.ninty,
#     0.4,
# )
# axs[3].text(
#     left + bird_widths.ninty,
#     0.99,
#     "$\\pi$",
#     ha="center",
#     transform=axs[3].transAxes,
# )

# rectangular(
#     axs[4],
#     left,
#     2 * bird_widths.ninty,
#     0.4,
# )
# axs[4].text(
#     left + bird_widths.ninty,
#     0.99,
#     "$\\pi$",
#     ha="center",
#     transform=axs[4].transAxes,
# )

# left += 2 * bird_widths.ninty
# axs[3].text(
#     left + bird_widths.one_over_sw / 2,
#     0.75,
#     "$\\nicefrac{1}{4f_{\\text{sw}}^{(2)}}$",
#     ha="center",
#     va="center",
#     transform=axs[3].transAxes,
# )

# left += bird_widths.one_over_sw
# axs[3].plot([left, left], [0.5, 1], color="k", ls=":", transform=axs[3].transAxes)

# axs[3].text(
#     left + bird_widths.t_one / 2,
#     0.75,
#     "$\\nicefrac{t^{(1)}}{2}$",
#     ha="center",
#     va="center",
#     transform=axs[3].transAxes,
# )

# left += bird_widths.t_one

# sine(
#     axs[5],
#     left,
#     bird_widths.grad,
#     0.1,
# )
# axs[5].text(
#     left + bird_widths.grad / 2,
#     0.8,
#     "$g$",
#     ha="center",
#     transform=axs[5].transAxes,
# )

# rectangular(
#     axs[3],
#     left,
#     2 * bird_widths.ninty,
#     0.4,
# )
# axs[3].text(
#     left + bird_widths.ninty,
#     0.99,
#     "$\\pi$",
#     ha="center",
#     transform=axs[3].transAxes,
# )
# left += 2 * bird_widths.ninty

# axs[3].text(
#     left + bird_widths.one_over_sw / 2,
#     0.75,
#     "$\\nicefrac{1}{4f_{\\text{sw}}^{(2)}}$",
#     ha="center",
#     va="center",
#     transform=axs[3].transAxes,
# )
# left += bird_widths.one_over_sw

# sine(
#     axs[5],
#     left,
#     bird_widths.grad,
#     0.4,
# )
# axs[5].text(
#     left + bird_widths.grad / 2,
#     0.8,
#     "$4g$",
#     ha="center",
#     transform=axs[5].transAxes,
# )

# left += bird_widths.grad

# rectangular(
#     axs[3],
#     left,
#     bird_widths.ninty,
#     0.4,
# )
# axs[3].text(
#     left + bird_widths.ninty / 2,
#     0.99,
#     "$\\nicefrac{\\pi}{2}$",
#     ha="center",
#     transform=axs[3].transAxes,
# )
# left += bird_widths.ninty

fig.savefig("figures/j_refocussing/j_refocussing.pdf")
