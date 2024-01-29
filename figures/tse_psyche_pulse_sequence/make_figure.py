# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Mon 29 Jan 2024 17:07:07 EST

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
        edgecolor="none",
        lw=0.,
        color="k",
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


def chirp(ax, left, width, height, label, type_):
    x = np.linspace(left, left + width, 10000)
    x_saltire = np.linspace(0, 1, 10000)
    envelope = height * (1 - np.abs(np.cos(np.pi * x_saltire)) ** 50)

    if type_ == "saltire":
        delta_F = 30
        delta_f = 0.
        phase = np.cos(
            np.pi * delta_F * (x_saltire - 0.5) ** 2 -
            2 * np.pi * delta_f * (x_saltire - 0.5)
        )
        y = np.abs(envelope * phase)
        ax.plot(x, y, color="k", lw=0.6)

    else:
        y = np.abs(envelope)
        ax.plot(x, y, color="k", lw=0.6)
        if type_ == "lohi":
            lft, right, bottom, top = (
                x[600],
                x[8500],
                0.02,
                height - 0.02,
            )
            xyA, xyB = (lft, bottom), (right, top)
        if type_ == "hilo":
            lft, right, bottom, top = (
                x[1500],
                x[9400],
                0.02,
                height - 0.02,
            )
            xyA, xyB = (lft, top), (right, bottom)

        arrow = ConnectionPatch(
            xyA,
            xyB,
            coordsA="data",
            coordsB="data",
            arrowstyle="-|>",
            facecolor="k",
            edgecolor="k",
            lw=0.6,
            mutation_scale=6,
        )
        ax.add_patch(arrow)

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
    def __init__(self, p1, p16, p20, p21, p22, t1, taua, taub, d1, d16, acqu):
        total = (
            p1 +
            6 * p16 +
            2 * p20 +
            p21 +
            p22 +
            2 * t1 +
            2 * taua +
            taub +
            d1 +
            6 * d16 +
            acqu
        )
        for k, v in locals().items():
            if k == "self":
                continue
            setattr(self, k, v / total)


widths = Widths(0.5, 1.2, 2.5, 5, 5, 5, 2, 2, 1.7, 1.7, 5)

delay(axs, "$d_1$", 0., widths.d1, left_line=False)
left = widths.d1
rect(axs[0], left, widths.p1, HARD_PULSE_HEIGHT, "$\\nicefrac{\\pi}{2}$ ($\\Phi_1$)")
left += widths.p1
delay(axs, "$\\tau_b$", left, widths.taub)
left += widths.taub
delay(axs, "$\\tau_a$", left, widths.taua)
left += widths.taua
sin(axs[1], left, widths.p16, 0.49, "$g_1$")
left += widths.p16
delay(axs, "$d_{16}$", left, widths.d16)
left += widths.d16
chirp(axs[0], left, widths.p21, 0.5, "$\\pi$ ($\\Phi_2$)", "lohi")
rect(axs[1], left, widths.p21, 0.02, "$g_2$")
left += widths.p21
sin(axs[1], left, widths.p16, 0.49, "$g_1$")
left += widths.p16
delay(axs, "$d_{16}$", left, widths.d16)
left += widths.d16
delay(axs, "$\\tau_a$", left, widths.taua)
left += widths.taua
delay(axs, "$\\nicefrac{t^{(1)}}{2}$", left, widths.t1)
left += widths.t1
sin(axs[1], left, widths.p16, 0.35, "$g_3$")
left += widths.p16
delay(axs, "$d_{16}$", left, widths.d16)
left += widths.d16
chirp(axs[0], left, widths.p20, 0.3, "", "saltire")
rect(axs[1], left, 2 * widths.p20, 0.03, "$g_4$")
left += widths.p20
axs[0].text(left, 0.36, "$\\beta$ ($\\Phi_3$)", ha="center")
chirp(axs[0], left, widths.p20, 0.3, "", "saltire")
left += widths.p20
sin(axs[1], left, widths.p16, 0.35, "$g_3$")
left += widths.p16
delay(axs, "$d_{16}$", left, widths.d16)
left += widths.d16
sin(axs[1], left, widths.p16, 0.77, "$g_5$")
left += widths.p16
delay(axs, "$d_{16}$", left, widths.d16)
left += widths.d16
chirp(axs[0], left, widths.p22, 0.5, "$\\pi$ ($\\Phi_4$)", "hilo")
rect(axs[1], left, widths.p22, 0.02, "$g_6$")
left += widths.p22
sin(axs[1], left, widths.p16, 0.77, "$g_5$")
left += widths.p16
delay(axs, "$d_{16}$", left, widths.d16)
left += widths.d16
delay(axs, "$\\nicefrac{t^{(1)}}{2}$", left, widths.t1)
left += widths.t1
acqu(axs[0], left, widths.acqu, 0.41)

# === Misc additions ===
axs[0].text(
    left + widths.acqu / 2, -0.35, "$\\nicefrac{1}{f_{\\text{sw}}^{(1)}}$", ha="center",
    va="center", bbox={"facecolor": "white", "edgecolor": "none", "pad": 2},
)
arrow = FancyArrowPatch(
    (left, -0.35),
    (left + widths.acqu, -0.35),
    arrowstyle="<->",
    mutation_scale=8,
    lw=mpl.rcParams["lines.linewidth"],
)
axs[0].add_artist(arrow)

axs[0].text(-0.015, 0.44, "\\textsuperscript{1}H", ha="center", va="center", fontsize=8)
axs[1].text(-0.015, 0.385, "$g_z$", ha="center", va="center", fontsize=8)
fig.savefig("figures/tse_psyche_pulse_sequence/tse_psyche_pulse_sequence.pdf")
