# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Fri 30 Jun 2023 16:44:09 BST

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch, FancyArrowPatch

mpl.rcParams["axes.xmargin"] = 0.
mpl.rcParams["font.size"] = 7
mpl.rcParams["lines.linewidth"] = 0.5

AX0_YLIM = (-0.6, 1.15)
AX1_YLIM = (-0.7, 0.7)
HARD_PULSE_HEIGHT = 0.84
GP6 = 0.5
GP7 = -0.1713 * GP6
GP8 = -0.1317 * GP6

with open("/home/simon/Documents/DPhil/data/andrographolide/2/gpnam1", "r") as fh:
    smsq = fh.readlines()[23:-1]
SMSQ = np.array([float(s.strip()) for s in smsq])

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


def sin(ax, left, width, height, label):
    x = np.linspace(left, left + width, 100)
    y = height * SMSQ
    ax.fill(x, y, color="k", clip_on=False, edgecolor="none")
    mid = (2 * left + width) / 2
    va, shift = ("bottom", 0.04) if height > 0 else ("top", -0.04)
    ax.text(
        mid, height + shift, label, ha="center", va=va,
        bbox={"facecolor": "w", "edgecolor": "none", "pad": 0.3},
    )

def diffgrad(ax, left, width, height, label):
    x = np.linspace(left, left + width, 100)
    y = height * SMSQ
    for i in (3, 2, 1):
        y *= (i / 3)
        col = str((i - 1) / 3)
        ax.fill(x, y, color=col, clip_on=False, edgecolor="none")
    mid = (2 * left + width) / 2
    va, shift = ("bottom", 0.02) if height > 0 else ("top", -0.02)
    ax.text(
        mid, height + shift, label, ha="center", va=va,
        bbox={"facecolor": "w", "edgecolor": "none", "pad": 0},
    )


def acqu(ax, left, width, height):
    x = np.linspace(left, left + width, 100)
    y = height * (np.cos(np.linspace(0, 10 * np.pi, 100)) * np.exp(np.linspace(0, -4, 100)))  # noqa: E501
    ax.plot(x, y, color="k", lw=0.8, solid_capstyle="round")
    mid = (2 * left + width) / 2
    ax.text(mid, HARD_PULSE_HEIGHT / 2, "$t^{(1)}$ ($\\Phi_{\\text{rx}}$)", ha="center", va="center")


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

# === Set up axes ===
fig, axs = plt.subplots(
    nrows=2,
    ncols=1,
    gridspec_kw={
        "left": 0.03,
        "right": 1.,
        "top": 0.99,
        "bottom": 0.01,
        "height_ratios": [AX0_YLIM[1] - AX0_YLIM[0], AX1_YLIM[1] - AX1_YLIM[0]],
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


class Widths:
    def __init__(self, p1, p19, p30, delta, led, d1, d16, acqu):
        total = (
            9 * p1 +
            2 * p19 +
            4 * p30 +
            delta +
            led +
            d1 +
            6 * d16 +
            acqu
        )
        for k, v in locals().items():
            if k == "self":
                continue
            setattr(self, k, v / total)


widths = Widths(0.1, 0.3, 0.3, 7., 3., 1., 1., 3.)

delay(axs, "$d_1$", 0., widths.d1, left_line=False)
left = widths.d1
rect(axs[0], left, widths.p1, HARD_PULSE_HEIGHT, "$\\nicefrac{\\pi}{2}$ ($\\Phi_1$)")
left += widths.p1
sin(axs[1], left, widths.p30, GP6, "$g_6$")
left += widths.p30
delay(axs, "$d_{16}$", left, widths.d16)
left += widths.d16
rect(axs[0], left, 2 * widths.p1, HARD_PULSE_HEIGHT, "$\\pi$ ($\\Phi_1$)")
left += 2 * widths.p1
sin(axs[1], left, widths.p30, -GP6, "$-g_6$")
left += widths.p30
delay(axs, "$d_{16}$", left, widths.d16)
left += widths.d16
rect(axs[0], left, widths.p1, HARD_PULSE_HEIGHT, "$\\nicefrac{\\pi}{2}$ ($\\Phi_2$)")
left += widths.p1
sin(axs[1], left, widths.p19, GP7, "$g_7$")
left += widths.p19
delay(axs, "$d_{16}$", left, widths.d16)
left += widths.d16
delay(axs, "", left, widths.delta)
left += widths.delta
rect(axs[0], left, widths.p1, HARD_PULSE_HEIGHT, "$\\nicefrac{\\pi}{2}$ ($\\Phi_3$)")
left += widths.p1
sin(axs[1], left, widths.p30, GP6, "$g_6$")
left += widths.p30
delay(axs, "$d_{16}$", left, widths.d16)
left += widths.d16
rect(axs[0], left, 2 * widths.p1, HARD_PULSE_HEIGHT, "$\\pi$ ($\\Phi_1$)")
left += 2 * widths.p1
sin(axs[1], left, widths.p30, -GP6, "$-g_6$")
left += widths.p30
delay(axs, "$d_{16}$", left, widths.d16)
left += widths.d16
rect(axs[0], left, widths.p1, HARD_PULSE_HEIGHT, "$\\nicefrac{\\pi}{2}$ ($\\Phi_4$)")
left += widths.p1
sin(axs[1], left, widths.p19, GP8, "$g_8$")
left += widths.p19
delay(axs, "$d_{16}$", left, widths.d16)
left += widths.d16
delay(axs, "", left, widths.led)
left += widths.led
rect(axs[0], left, widths.p1, HARD_PULSE_HEIGHT, "$\\nicefrac{\\pi}{2}$ ($\\Phi_5$)")
left += widths.p1
acqu(axs[0], left, widths.acqu, 0.41)

# === Misc additions ===
# lefts = [
#     widths.d1 + widths.p19 + widths.d16 + widths.p1 + 0.5 * widths.p30,
#     widths.d1 + 2 * widths.p19 + 6 * widths.d16 + 5 * widths.p1 + 4 * widths.p30 + widths.delta,
# ]
# tau = 2 * widths.p1 + widths.d16 + widths.p30
# for left in lefts:
#     axs[0].text(
#         left + 0.5 * tau, -0.35, "$\\tau$", ha="center",
#         va="center", bbox={"facecolor": "white", "edgecolor": "none", "pad": 2},
#     )
#     arrow = FancyArrowPatch(
#         (left, -0.35),
#         (left + tau, -0.35),
#         arrowstyle="<->",
#         mutation_scale=8,
#         lw=mpl.rcParams["lines.linewidth"],
#     )
#     axs[0].add_artist(arrow)

left = widths.d1
difftime = widths.delta + 4 * widths.p1 + 2 * widths.p30 + 3 * widths.d16 + widths.p19
axs[0].text(
    left + 0.5 * difftime, 1.1, "$d_{20}$", ha="center",
    va="center", bbox={"facecolor": "white", "edgecolor": "none", "pad": 2},
)
arrow = FancyArrowPatch(
    (left, 1.1),
    (left + difftime, 1.1),
    arrowstyle="<->",
    mutation_scale=8,
    lw=mpl.rcParams["lines.linewidth"],
    shrinkA=0.,
    shrinkB=0.,
)
axs[0].add_artist(arrow)

left = 1 - widths.acqu -  widths.p1 - widths.led - widths.d16 - widths.p19
ledtime = widths.led + widths.d16 + widths.p19
axs[0].text(
    left + 0.5 * ledtime, 1.1, "$d_{21}$", ha="center",
    va="center", bbox={"facecolor": "white", "edgecolor": "none", "pad": 2},
)
arrow = FancyArrowPatch(
    (left, 1.1),
    (left + ledtime, 1.1),
    arrowstyle="<->",
    mutation_scale=8,
    lw=mpl.rcParams["lines.linewidth"],
    shrinkA=0.,
    shrinkB=0.,
)
axs[0].add_artist(arrow)

axs[0].text(-0.015, 0.46, "\\textsuperscript{1}H", ha="center", va="center", fontsize=8)
axs[1].text(-0.015, 0., "$g_z$", ha="center", va="center", fontsize=8)




fig.savefig("figures/ledbpgp_pulse_sequence/ledbpgp_pulse_sequence.pdf")

