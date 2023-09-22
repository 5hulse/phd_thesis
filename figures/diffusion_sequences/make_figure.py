# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Wed 20 Sep 2023 19:16:28 BST

from dataclasses import dataclass
import matplotlib as mpl
from matplotlib.patches import FancyArrowPatch
import matplotlib.pyplot as plt
import numpy as np

mpl.rcParams["font.size"] = 8

fig, axs = plt.subplots(
    nrows=4,
    ncols=1,
    gridspec_kw={
        "left": 0.05,
        "right": 0.98,
        "bottom": 0.,
        "top": 0.96,
        "hspace": 0.08,
    },
    figsize=(
        mpl.rcParams["figure.figsize"][0],
        4,
    ),
)

for i, ax in enumerate(axs):
    for x in ("top", "bottom", "left", "right"):
        ax.spines[x].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axhline(0.5, color="k")
    ax.text(
        -0.045,
        0.98,
        f"\\textbf{{{chr(97 + i)}}}.",
        fontsize=12,
        transform=ax.transAxes,
    )


def normalise_widths(widths):
    total_w = sum(widths)
    return [w / total_w for w in widths]


def rectangular(ax, left, width, height):
    right = left + width
    xs = [left, left, right, right]
    ys = [0.5, 0.5 + height, 0.5 + height, 0.5]
    ax.fill(xs, ys, "k", transform=ax.transAxes)


def sinebell(ax, left, width, height, neg=False, color="k"):
    right = left + width
    # mean = (left + right) / 2
    xs = np.linspace(left, right, 100)
    ys = np.sin(np.linspace(0, np.pi, 100))
    scale = height / np.amax(ys)
    if neg:
        ys *= -1
    ys = ys * scale + 0.5
    ax.fill(xs, ys, "k", transform=ax.transAxes, color=color)


def acquisition(ax, left, width, height):
    right = left + width
    xs = np.linspace(left, right, 100)
    xs_fid = np.linspace(0, 30, 100)
    ys = np.cos(xs_fid) * np.exp(0.1 * -xs_fid)
    scale = (height) / ys[0]
    ys = ys * scale + 0.5
    ax.plot(xs, ys, color="k", transform=ax.transAxes)


class PGSE:
    def __init__(
        self,
        ninty: float,
        oneeighty: float,
        pre_grad_delay: float,
        grad: float,
        post_grad_delay: float,
        acqu: float,
        left_pad: float,
        right_pad: float,
    ):
        (
            self.ninty,
            self.oneeighty,
            self.pre_grad_delay,
            self.grad,
            self.post_grad_delay,
            self.acqu,
            self.left_pad,
            self.right_pad,
        ) = normalise_widths([
            ninty,
            oneeighty,
            pre_grad_delay,
            grad,
            post_grad_delay,
            acqu,
            left_pad,
            right_pad,
        ])
        self.pre_grad_delay /= 2
        self.grad /= 2
        self.post_grad_delay /= 2


pgse_widths = PGSE(*normalise_widths([1., 2., 1., 2., 29., 4., 0.5, 0.5]))

left = pgse_widths.left_pad

# 90
rectangular(
    axs[0],
    left,
    pgse_widths.ninty,
    0.3,
)
axs[0].text(
    left + pgse_widths.ninty / 2,
    0.84,
    "$\\nicefrac{\\pi}{2}$",
    ha="center",
    transform=axs[0].transAxes,
)

left += pgse_widths.ninty + pgse_widths.pre_grad_delay

Delta_w = (
    pgse_widths.grad +
    pgse_widths.post_grad_delay +
    pgse_widths.oneeighty +
    pgse_widths.pre_grad_delay
)
arrow = FancyArrowPatch(
    (left, 0.96),
    (left + Delta_w, 0.96),
    arrowstyle="|-|",
    lw=0.8,
    shrinkA=0,
    shrinkB=0,
    transform=axs[0].transAxes,
)
axs[0].add_patch(arrow)
axs[0].text(
    left + (Delta_w / 2),
    0.99,
    "$\\Updelta$",
    ha="center",
    transform=axs[0].transAxes,
)

for i in range(3, 0, -1):
    h = 0.1 * i
    color = 3 * [(i - 1) * 0.33]
    sinebell(
        axs[0], left,
        pgse_widths.grad,
        h,
        color=color,
    )
    if i == 3:
        arrowstyle = "->"
    elif i == 2:
        arrowstyle = "-"
    elif i == 1:
        arrowstyle = "<-"
    arrow = FancyArrowPatch(
        (left + pgse_widths.grad + 0.01, 0.5),
        (left + pgse_widths.grad + 0.01, 0.5 + h),
        arrowstyle=arrowstyle,
        mutation_scale=5,
        lw=0.8,
        shrinkA=0,
        shrinkB=0,
        transform=axs[0].transAxes,
        color=color,
    )
    axs[0].add_patch(arrow)
axs[0].text(
    left + pgse_widths.grad + 0.015,
    0.65,
    "$g$",
    va="center",
    transform=axs[0].transAxes,
)

axs[0].text(
    left + pgse_widths.grad / 2,
    0.3,
    "$\\delta$",
    ha="center",
    transform=axs[0].transAxes,
)

arrow = FancyArrowPatch(
    (left, 0.45),
    (left + pgse_widths.grad, 0.45),
    arrowstyle="|-|",
    lw=0.8,
    shrinkA=0,
    shrinkB=0,
    transform=axs[0].transAxes,
)
axs[0].add_patch(arrow)


left += pgse_widths.grad + pgse_widths.post_grad_delay

rectangular(
    axs[0],
    left,
    pgse_widths.oneeighty,
    0.3,
)
axs[0].text(
    left + pgse_widths.oneeighty / 2,
    0.84,
    "$\\pi$",
    ha="center",
    transform=axs[0].transAxes,
)

left += pgse_widths.oneeighty + pgse_widths.pre_grad_delay

for i in range(3, 0, -1):
    h = 0.1 * i
    color = 3 * [(i - 1) * 0.33]
    sinebell(
        axs[0],
        left,
        pgse_widths.grad,
        h,
        color=color,
    )
    if i == 3:
        arrowstyle = "->"
    elif i == 2:
        arrowstyle = "-"
    elif i == 1:
        arrowstyle = "<-"
    arrow = FancyArrowPatch(
        (left + pgse_widths.grad + 0.01, 0.5),
        (left + pgse_widths.grad + 0.01, 0.5 + h),
        arrowstyle=arrowstyle,
        mutation_scale=5,
        lw=0.8,
        shrinkA=0,
        shrinkB=0,
        transform=axs[0].transAxes,
        color=color,
    )
    axs[0].add_patch(arrow)
axs[0].text(
    left + pgse_widths.grad + 0.015,
    0.65,
    "$g$",
    va="center",
    transform=axs[0].transAxes,
)

axs[0].text(
    left + pgse_widths.grad / 2,
    0.3,
    "$\\delta$",
    ha="center",
    transform=axs[0].transAxes,
)

arrow = FancyArrowPatch(
    (left, 0.45),
    (left + pgse_widths.grad, 0.45),
    arrowstyle="|-|",
    lw=0.8,
    shrinkA=0,
    shrinkB=0,
    transform=axs[0].transAxes,
)
axs[0].add_artist(arrow)

left += pgse_widths.grad + pgse_widths.post_grad_delay

acquisition(
    axs[0],
    left,
    pgse_widths.acqu,
    0.4,
)


class PGSTE:

    def __init__(
        self,
        ninty: float,
        pre_grad_delay: float,
        grad: float,
        post_grad_outer_delay: float,
        post_grad_inner_delay: float,
        acqu: float,
        left_pad: float,
        right_pad: float,
    ):
        (
            self.ninty,
            self.pre_grad_delay,
            self.grad,
            self.post_grad_outer_delay,
            self.post_grad_inner_delay,
            self.acqu,
            self.left_pad,
            self.right_pad,
        ) = normalise_widths([
            ninty,
            pre_grad_delay,
            grad,
            post_grad_outer_delay,
            post_grad_inner_delay,
            acqu,
            left_pad,
            right_pad,
        ])
        self.ninty /= 3
        self.pre_grad_delay /= 3
        self.grad /= 3
        self.post_grad_outer_delay /= 2


pgste_widths = PGSTE(2., 2., 2., 10., 19., 4., 0.5, 0.5)

left = pgste_widths.left_pad
# 90
rectangular(
    axs[1],
    left,
    pgste_widths.ninty,
    0.3,
)
axs[1].text(
    left + pgste_widths.ninty / 2,
    0.84,
    "$\\nicefrac{\\pi}{2}$",
    ha="center",
    transform=axs[1].transAxes,
)

left += pgste_widths.ninty + pgste_widths.pre_grad_delay

for i in range(3, 0, -1):
    h = 0.1 * i
    color = 3 * [(i - 1) * 0.33]
    sinebell(
        axs[1],
        left,
        pgse_widths.grad,
        h,
        color=color,
    )
    if i == 3:
        arrowstyle = "->"
    elif i == 2:
        arrowstyle = "-"
    elif i == 1:
        arrowstyle = "<-"
    arrow = FancyArrowPatch(
        (left + pgse_widths.grad + 0.01, 0.5),
        (left + pgse_widths.grad + 0.01, 0.5 + h),
        arrowstyle=arrowstyle,
        mutation_scale=5,
        lw=0.8,
        shrinkA=0,
        shrinkB=0,
        transform=axs[1].transAxes,
        color=color,
    )
    axs[1].add_patch(arrow)
axs[1].text(
    left + pgse_widths.grad + 0.015,
    0.65,
    "$g$",
    va="center",
    transform=axs[1].transAxes,
)

axs[1].text(
    left + pgse_widths.grad / 2,
    0.3,
    "$\\delta$",
    ha="center",
    transform=axs[1].transAxes,
)

arrow = FancyArrowPatch(
    (left, 0.45),
    (left + pgse_widths.grad, 0.45),
    arrowstyle="|-|",
    lw=0.8,
    shrinkA=0,
    shrinkB=0,
    transform=axs[1].transAxes,
)
axs[1].add_artist(arrow)

Delta_w = (
    pgste_widths.post_grad_outer_delay +
    2 * pgste_widths.ninty +
    2 * pgste_widths.pre_grad_delay +
    2 * pgste_widths.grad +
    pgste_widths.post_grad_inner_delay
)
arrow = FancyArrowPatch(
    (left, 0.96),
    (left + Delta_w, 0.96),
    arrowstyle="|-|",
    lw=0.8,
    shrinkA=0,
    shrinkB=0,
    transform=axs[1].transAxes,
)
axs[1].add_patch(arrow)
axs[1].text(
    left + (Delta_w / 2),
    0.99,
    "$\\Updelta$",
    ha="center",
    transform=axs[1].transAxes,
)

left += pgste_widths.grad + pgste_widths.post_grad_outer_delay
# 90
rectangular(
    axs[1],
    left,
    pgste_widths.ninty,
    0.3,
)
axs[1].text(
    left + pgste_widths.ninty / 2,
    0.84,
    "$\\nicefrac{\\pi}{2}$",
    ha="center",
    transform=axs[1].transAxes,
)

left += pgste_widths.ninty + pgste_widths.pre_grad_delay
sinebell(
    axs[1],
    left,
    pgste_widths.grad,
    0.066,
)

axs[1].text(
    left + pgste_widths.grad / 2,
    0.64,
    "$g_{\\text{spoil}}$",
    ha="center",
    transform=axs[1].transAxes,
)


left += pgste_widths.grad + pgste_widths.post_grad_inner_delay
# 90
rectangular(
    axs[1],
    left,
    pgste_widths.ninty,
    0.3,
)
axs[1].text(
    left + pgste_widths.ninty / 2,
    0.84,
    "$\\nicefrac{\\pi}{2}$",
    ha="center",
    transform=axs[1].transAxes,
)

left += pgste_widths.ninty + pgste_widths.pre_grad_delay

for i in range(3, 0, -1):
    h = 0.1 * i
    color = 3 * [(i - 1) * 0.33]
    sinebell(
        axs[1],
        left,
        pgse_widths.grad,
        h,
        color=color,
    )
    if i == 3:
        arrowstyle = "->"
    elif i == 2:
        arrowstyle = "-"
    elif i == 1:
        arrowstyle = "<-"
    arrow = FancyArrowPatch(
        (left + pgse_widths.grad + 0.01, 0.5),
        (left + pgse_widths.grad + 0.01, 0.5 + h),
        arrowstyle=arrowstyle,
        mutation_scale=5,
        lw=0.8,
        shrinkA=0,
        shrinkB=0,
        transform=axs[1].transAxes,
        color=color,
    )
    axs[1].add_patch(arrow)

axs[1].text(
    left + pgse_widths.grad + 0.015,
    0.65,
    "$g$",
    va="center",
    transform=axs[1].transAxes,
)

axs[1].text(
    left + pgse_widths.grad / 2,
    0.3,
    "$\\delta$",
    ha="center",
    transform=axs[1].transAxes,
)

arrow = FancyArrowPatch(
    (left, 0.45),
    (left + pgse_widths.grad, 0.45),
    arrowstyle="|-|",
    lw=0.8,
    shrinkA=0,
    shrinkB=0,
    transform=axs[1].transAxes,
)
axs[1].add_artist(arrow)

left += pgste_widths.grad + pgste_widths.post_grad_outer_delay

acquisition(
    axs[1],
    left,
    pgste_widths.acqu,
    0.4,
)


class PGSTEBP:

    def __init__(
        self,
        ninty: float,
        oneeighty: float,
        pre_grad_delay: float,
        grad: float,
        post_grad_outer_delay: float,
        post_grad_inner_delay: float,
        acqu: float,
        left_pad: float,
        right_pad: float,
    ):
        (
            self.ninty,
            self.oneeighty,
            self.pre_grad_delay,
            self.grad,
            self.post_grad_outer_delay,
            self.post_grad_inner_delay,
            self.acqu,
            self.left_pad,
            self.right_pad,
        ) = normalise_widths([
            ninty,
            oneeighty,
            pre_grad_delay,
            grad,
            post_grad_outer_delay,
            post_grad_inner_delay,
            acqu,
            left_pad,
            right_pad,
        ])
        self.ninty /= 3
        self.oneeighty /= 2
        self.grad /= 5
        self.pre_grad_delay /= 5
        self.post_grad_outer_delay /= 4


# PGSETBP and Oneshot
pgstebp_widths = PGSTEBP(2., 3., 3., 3., 15., 9., 4., 0.5, 0.5)
for i, ax in enumerate(axs[2:]):
    if i == 0:
        scale_up, scale_down = 0, 0
        label_up, label_down = "$g$", "$g$"
    elif i == 1:
        scale_up, scale_down = 0.02, -0.02
        label_up, label_down = "$g(1+\\alpha)$", "$g(1-\\alpha)$"

    left = pgstebp_widths.left_pad

    # 90
    rectangular(
        ax,
        left,
        pgstebp_widths.ninty,
        0.3,
    )
    ax.text(
        left + (pgstebp_widths.ninty) / 2,
        0.84,
        "$\\nicefrac{\\pi}{2}$",
        ha="center",
        transform=ax.transAxes,
    )

    left += pgstebp_widths.ninty + pgstebp_widths.pre_grad_delay

    Delta_w = (
        3 * pgstebp_widths.grad +
        2 * pgstebp_widths.post_grad_outer_delay +
        2 * pgstebp_widths.ninty +
        pgstebp_widths.oneeighty +
        3 * pgstebp_widths.pre_grad_delay +
        pgstebp_widths.post_grad_inner_delay
    )
    arrow = FancyArrowPatch(
        (left, 0.96),
        (left + Delta_w, 0.96),
        arrowstyle="|-|",
        lw=0.8,
        shrinkA=0,
        shrinkB=0,
        transform=ax.transAxes,
    )
    ax.add_patch(arrow)
    ax.text(
        left + (Delta_w / 2),
        0.99,
        "$\\Updelta$",
        ha="center",
        transform=ax.transAxes,
    )

    for i in range(3, 0, -1):
        h = (0.1 + scale_up) * i
        color = 3 * [(i - 1) * 0.33]
        sinebell(
            ax,
            left,
            pgstebp_widths.grad,
            h,
            color=color,
        )
        if i == 3:
            arrowstyle = "->"
        elif i == 2:
            arrowstyle = "-"
        elif i == 1:
            arrowstyle = "<-"
        arrow = FancyArrowPatch(
            (left + pgse_widths.grad + 0.0025, 0.5),
            (left + pgse_widths.grad + 0.0025, 0.5 + h),
            arrowstyle=arrowstyle,
            mutation_scale=5,
            lw=0.8,
            shrinkA=0,
            shrinkB=0,
            transform=ax.transAxes,
            color=color,
        )
        ax.add_patch(arrow)
    ax.text(
        left + pgse_widths.grad + 0.0075,
        (0.5 + 0.5 + (0.1 + scale_up) * 3) / 2,
        label_up,
        va="center",
        transform=ax.transAxes,
    )

    arrow = FancyArrowPatch(
        (left, 0.45),
        (left + pgstebp_widths.grad, 0.45),
        arrowstyle="|-|",
        lw=0.8,
        shrinkA=0,
        shrinkB=0,
        transform=ax.transAxes,
    )
    ax.add_patch(arrow)

    ax.text(
        left + pgstebp_widths.grad / 2,
        0.3,
        "$\\nicefrac{\\delta}{2}$",
        ha="center",
        transform=ax.transAxes,
    )

    left += pgstebp_widths.grad
    tau_w = pgstebp_widths.post_grad_outer_delay

    arrow = FancyArrowPatch(
        (left, 0.2),
        (left + tau_w, 0.2),
        arrowstyle="|-|",
        lw=0.8,
        shrinkA=0,
        shrinkB=0,
        transform=ax.transAxes,
    )
    ax.add_patch(arrow)
    ax.text(
        left + (tau_w / 2),
        0.1,
        "$\\tau$",
        ha="center",
        transform=ax.transAxes,
    )

    left += pgstebp_widths.post_grad_outer_delay

    rectangular(
        ax,
        left,
        pgstebp_widths.oneeighty,
        0.3,
    )
    ax.text(
        left + (pgstebp_widths.oneeighty) / 2,
        0.84,
        "$\\pi$",
        ha="center",
        transform=ax.transAxes,
    )

    left += pgstebp_widths.oneeighty + pgstebp_widths.pre_grad_delay

    for i in range(3, 0, -1):
        h = (0.1 + scale_down) * i
        color = 3 * [(i - 1) * 0.33]
        sinebell(
            ax,
            left,
            pgstebp_widths.grad,
            h,
            color=color,
            neg=True,
        )
        if i == 3:
            arrowstyle = "->"
        elif i == 2:
            arrowstyle = "-"
        elif i == 1:
            arrowstyle = "<-"
        arrow = FancyArrowPatch(
            (left + pgse_widths.grad + 0.0025, 0.5),
            (left + pgse_widths.grad + 0.0025, 0.5 - h),
            arrowstyle=arrowstyle,
            mutation_scale=5,
            lw=0.8,
            shrinkA=0,
            shrinkB=0,
            transform=ax.transAxes,
            color=color,
        )
        ax.add_patch(arrow)
    ax.text(
        left + pgse_widths.grad + 0.0075,
        (0.5 + 0.5 - (0.1 + scale_down) * 3) / 2,
        label_down,
        va="center",
        transform=ax.transAxes,
    )



    arrow = FancyArrowPatch(
        (left, 0.55),
        (left + pgstebp_widths.grad, 0.55),
        arrowstyle="|-|",
        lw=0.8,
        shrinkA=0,
        shrinkB=0,
        transform=ax.transAxes,
        )
    ax.add_patch(arrow)

    ax.text(
        left + pgstebp_widths.grad / 2,
        0.7,
        "$\\nicefrac{\\delta}{2}$",
        ha="center",
        va="top",
        transform=ax.transAxes,
    )
    left += pgstebp_widths.grad + pgstebp_widths.post_grad_outer_delay

    rectangular(
        ax,
        left,
        pgstebp_widths.ninty,
        0.3,
    )
    ax.text(
        left + pgstebp_widths.ninty / 2,
        0.84,
        "$\\nicefrac{\\pi}{2}$",
        ha="center",
        transform=ax.transAxes,
    )

    left += pgstebp_widths.ninty + pgstebp_widths.pre_grad_delay

    sinebell(
        ax,
        left,
        pgstebp_widths.grad,
        0.066,
    )
    ax.text(
        left + pgste_widths.grad / 2,
        0.64,
        "$g_{\\text{spoil}}$",
        ha="center",
        transform=ax.transAxes,
    )

    left += pgstebp_widths.grad + pgstebp_widths.post_grad_inner_delay

    rectangular(
        ax,
        left,
        pgstebp_widths.ninty,
        0.3,
    )
    ax.text(
        left + (pgstebp_widths.ninty) / 2,
        0.84,
        "$\\nicefrac{\\pi}{2}$",
        ha="center",
        transform=ax.transAxes,
    )

    left += pgstebp_widths.ninty + pgstebp_widths.pre_grad_delay

    for i in range(3, 0, -1):
        h = (0.1 + scale_up) * i
        color = 3 * [(i - 1) * 0.33]
        sinebell(
            ax,
            left,
            pgstebp_widths.grad,
            h,
            color=color,
        )
        if i == 3:
            arrowstyle = "->"
        elif i == 2:
            arrowstyle = "-"
        elif i == 1:
            arrowstyle = "<-"
        arrow = FancyArrowPatch(
            (left + pgse_widths.grad + 0.0025, 0.5),
            (left + pgse_widths.grad + 0.0025, 0.5 + h),
            arrowstyle=arrowstyle,
            mutation_scale=5,
            lw=0.8,
            shrinkA=0,
            shrinkB=0,
            transform=ax.transAxes,
            color=color,
        )
        ax.add_patch(arrow)
    ax.text(
        left + pgse_widths.grad + 0.0075,
        (0.5 + 0.5 + (0.1 + scale_up) * 3) / 2,
        label_up,
        va="center",
        transform=ax.transAxes,
    )



    arrow = FancyArrowPatch(
        (left, 0.45),
        (left + pgstebp_widths.grad, 0.45),
        arrowstyle="|-|",
        lw=0.8,
        shrinkA=0,
        shrinkB=0,
        transform=ax.transAxes,
        )
    ax.add_patch(arrow)

    ax.text(
        left + pgstebp_widths.grad / 2,
        0.3,
        "$\\nicefrac{\\delta}{2}$",
        ha="center",
        transform=ax.transAxes,
    )

    left += pgstebp_widths.grad
    tau_w = pgstebp_widths.post_grad_outer_delay

    arrow = FancyArrowPatch(
        (left, 0.2),
        (left + tau_w, 0.2),
        arrowstyle="|-|",
        lw=0.8,
        shrinkA=0,
        shrinkB=0,
        transform=ax.transAxes,
    )
    ax.add_patch(arrow)
    ax.text(
        left + (tau_w / 2),
        0.1,
        "$\\tau$",
        ha="center",
        transform=ax.transAxes,
    )

    left += pgstebp_widths.post_grad_outer_delay

    rectangular(
        ax,
        left,
        pgstebp_widths.oneeighty,
        0.3,
    )
    ax.text(
        left + (pgstebp_widths.oneeighty) / 2,
        0.84,
        "$\\pi$",
        ha="center",
        transform=ax.transAxes,
    )

    left += pgstebp_widths.oneeighty + pgstebp_widths.pre_grad_delay

    for i in range(3, 0, -1):
        h = (0.1 + scale_down) * i
        color = 3 * [(i - 1) * 0.33]
        sinebell(
            ax,
            left,
            pgstebp_widths.grad,
            h,
            color=color,
            neg=True,
        )
        if i == 3:
            arrowstyle = "->"
        elif i == 2:
            arrowstyle = "-"
        elif i == 1:
            arrowstyle = "<-"
        arrow = FancyArrowPatch(
            (left + pgse_widths.grad + 0.0025, 0.5),
            (left + pgse_widths.grad + 0.0025, 0.5 - h),
            arrowstyle=arrowstyle,
            mutation_scale=5,
            lw=0.8,
            shrinkA=0,
            shrinkB=0,
            transform=ax.transAxes,
            color=color,
        )
        ax.add_patch(arrow)
    ax.text(
        left + pgse_widths.grad + 0.0075,
        (0.5 + 0.5 - (0.1 + scale_down) * 3) / 2,
        label_down,
        va="center",
        transform=ax.transAxes,
    )


    arrow = FancyArrowPatch(
        (left, 0.55),
        (left + pgstebp_widths.grad, 0.55),
        arrowstyle="|-|",
        lw=0.8,
        shrinkA=0,
        shrinkB=0,
        transform=ax.transAxes,
    )
    ax.add_patch(arrow)

    ax.text(
        left + pgstebp_widths.grad / 2,
        0.7,
        "$\\nicefrac{\\delta}{2}$",
        ha="center",
        va="top",
        transform=ax.transAxes,
    )
    left += pgstebp_widths.grad + pgstebp_widths.post_grad_outer_delay

    acquisition(
        ax,
        left,
        pgste_widths.acqu,
        0.4,
    )

fig.savefig("figures/diffusion_sequences/diffusion_sequences.pdf")
