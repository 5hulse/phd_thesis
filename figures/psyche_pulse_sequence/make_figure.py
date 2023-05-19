# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Fri 19 May 2023 00:32:24 BST

import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.patches import ConnectionPatch, FancyArrowPatch

mpl.rcParams["axes.xmargin"] = 0.
mpl.rcParams["font.size"] = 7
mpl.rcParams["lines.linewidth"] = 0.5

AX0_YLIM = (-0.6, 1)
AX1_YLIM = (0., 0.47)
HARD_PULSE_HEIGHT = 0.84


def rect(ax, left, width, height, label):
    right = left + width
    ax.fill(
        [left, right, right, left],
        [0, 0, height, height],
        color="k",
        edgecolor="none",
    )
    mid = (2 * left + width) / 2
    ax.text(mid, height + 0.06, label, ha="center")


def delay(axs, text, left, width, left_line=True, right_line=True):
    mid = (2 * left + width) / 2
    axs[0].text(mid, HARD_PULSE_HEIGHT / 2, text, ha="center", va="center")
    right = left + width
    if left_line:
        line = ConnectionPatch(
            xyA=[left, HARD_PULSE_HEIGHT], xyB=[left, AX1_YLIM[0]], coordsA="data",
            coordsB="data", axesA=axs[0], axesB=axs[1], color="k", ls=":",
            lw=mpl.rcParams["lines.linewidth"]
        )
        axs[0].add_artist(line)
    if right_line:
        line = ConnectionPatch(
            xyA=[right, HARD_PULSE_HEIGHT], xyB=[right, AX1_YLIM[0]], coordsA="data",
            coordsB="data", axesA=axs[0], axesB=axs[1], color="k", ls=":",
            lw=mpl.rcParams["lines.linewidth"]
        )
        axs[0].add_artist(line)


def saltire(ax, left, width, height, label):
    x = np.linspace(left, left + width, 10000)
    x_saltire = np.linspace(0, 1, 10000)
    y = np.abs(
        height * (
            (1 - np.abs(np.cos(np.pi * x_saltire)) ** 50) *
            np.cos(
                60 * np.pi * (x_saltire - 0.5) ** 2 -
                2 * np.pi * 0.5 * (x_saltire - 0.5)
            )
        )
    )
    ax.plot(x, y, color="k")
    mid = (2 * left + width) / 2
    ax.text(mid, height + 0.06, label, ha="center")


def sin(ax, left, width, height, label):
    x = np.linspace(left, left + width, 100)
    y = height * (np.sin(np.linspace(0, np.pi, 100)))
    ax.fill(x, y, color="k", clip_on=False, edgecolor="none")
    mid = (2 * left + width) / 2
    ax.text(mid, height + 0.06, label, ha="center")


def acqu(ax, left, width, height):
    x = np.linspace(left, left + width, 100)
    y = height * (np.cos(np.linspace(0, 10 * np.pi, 100)) * np.exp(np.linspace(0, -4, 100)))  # noqa: E501
    ax.plot(x, y, color="k", lw=0.8, solid_capstyle="round")
    mid = (2 * left + width) / 2
    ax.text(mid, HARD_PULSE_HEIGHT / 2, "$t^{(2)}$ ($\\Phi_{\\text{rx}}$)", ha="center", va="center")


# === Set up axes ===
fig, axs = plt.subplots(
    nrows=2,
    ncols=1,
    gridspec_kw={
        "left": 0.03,
        "right": 1.,
        "top": 0.99,
        "bottom": 0.01,
        "height_ratios": [1.6, 0.47],
        "hspace": 0.,
    },
    figsize=(6, 2),
)

for i, ax in enumerate(axs):
    for pos in ("top", "bottom", "left", "right"):
        ax.spines[pos].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(0, 1)
    ax.axhline(0, color="k", clip_on=False)

axs[0].set_ylim(AX0_YLIM)
axs[1].set_ylim(AX1_YLIM)


# === Configure element widths ===
class Widths:
    def __init__(self, d0, p1, p16, p50, t1, taua, taud, d16, d18, d19, acqu):
        total = (
            d0 +
            3 * p1 +
            2 * t1 +
            2 * taua +
            4 * p16 +
            6 * d16 +
            taud +
            2 * d18 +
            2 * d19 +
            2 * p50 +
            acqu
        )
        for k, v in locals().items():
            if k == "self":
                continue
            setattr(self, k, v / total)

    @property
    def p2(self):
        return 2 * self.p1

    @property
    def p32(self):
        return 2 * self.d19 + 2 * self.p50


widths = Widths(2, 0.5, 1.2, 4, 3, 2, 2, 2, 2, 2, 7)

delay(axs, "$d_1$", 0., widths.d0, left_line=False)
left = widths.d0
rect(axs[0], left, widths.p1, HARD_PULSE_HEIGHT, "$\\nicefrac{\\pi}{2}$ ($\\Phi_1$)")
left += widths.p1
delay(axs, "$\\nicefrac{t^{(1)}}{2}$", left, widths.t1)
left += widths.t1
delay(axs, "$\\tau_a$", left, widths.taua)
left += widths.taua
sin(axs[1], left, widths.p16, 0.31, "$g_1$")
left += widths.p16
delay(axs, "$d_{16}$", left, widths.d16)
left += widths.d16
rect(axs[0], left, widths.p2, HARD_PULSE_HEIGHT, "$\\pi$ ($\\Phi_2$)")
left += widths.p2
delay(axs, "", left, 0., right_line=False)
sin(axs[1], left, widths.p16, 0.31, "$g_1$")
left += widths.p16
delay(axs, "$d_{16}$", left, widths.d16)
left += widths.d16
delay(axs, "$\\tau_a$", left, widths.taua)
left += widths.taua
delay(axs, "$\\tau_d$", left, widths.taua)
left += widths.taud
delay(axs, "$d_{16}$", left, widths.d16)
left += widths.d16
sin(axs[1], left, widths.p16, 0.43, "$g_2$")
left += widths.p16
delay(axs, "$d_{16}$", left, widths.d16)
left += widths.d16
delay(axs, "$d_{18}$", left, widths.d18)
left += widths.d18
rect(axs[1], left, widths.p32, 0.05, "$g_3$")
delay(axs, "$d_{19}$", left, widths.d19)
left += widths.d19
saltire(axs[0], left, widths.p50, 0.3, "")
left += widths.p50
axs[0].text(left, 0.36, "$\\beta$ ($\\Phi_3$)", ha="center")
saltire(axs[0], left, widths.p50, 0.3, "")
left += widths.p50
delay(axs, "$d_{19}$", left, widths.d19)
left += widths.d19
delay(axs, "$d_{18}$", left, widths.d18)
left += widths.d18
delay(axs, "$d_{16}$", left, widths.d16)
left += widths.d16
sin(axs[1], left, widths.p16, 0.43, "$g_2$")
left += widths.p16
delay(axs, "$d_{16}$", left, widths.d16)
left += widths.d16
delay(axs, "$\\nicefrac{t^{(1)}}{2}$", left, widths.t1)
left += widths.t1
acqu(axs[0], left, widths.acqu, 0.41)

# === Misc additions ===
axs[0].text(
    left + widths.acqu / 2, -0.45, "$\\nicefrac{1}{f_{\\text{sw}}^{(1)}}$", ha="center",
    va="center", bbox={"facecolor": "white", "edgecolor": "none", "pad": 2},
)
arrow = FancyArrowPatch(
    (left, -0.45),
    (left + widths.acqu, -0.45),
    arrowstyle="<->",
    mutation_scale=8,
    lw=0.8,
)
axs[0].add_artist(arrow)

left = widths.d0 + widths.p1 + widths.t1
width = 2 * (widths.taua + widths.d16 + widths.p16 + widths.p1)
axs[0].text(
    left + width / 2, -0.45, "$\\nicefrac{1}{2f_{\\text{sw}}^{(1)}}$", ha="center", va="center",
    bbox={"facecolor": "white", "edgecolor": "none", "pad": 2},
)
arrow = FancyArrowPatch(
    (left, -0.45),
    (left + width, -0.45),
    arrowstyle="<->",
    mutation_scale=8,
    lw=0.8,
)
axs[0].add_artist(arrow)

axs[0].text(-0.017, 0.43, "\\textsuperscript{1}H", va="center", ha="center", fontsize=8)
axs[1].text(-0.017, 0.27, "$g_z$", va="center", ha="center", fontsize=8)
fig.savefig("figures/psyche_pulse_sequence/psyche_pulse_sequence.pdf")
