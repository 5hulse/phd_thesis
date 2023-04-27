# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Mon 20 Mar 2023 19:41:50 GMT

import numpy as np
import matplotlib as mpl
from matplotlib import cm, lines, pyplot as plt, transforms


hbar = 6.62607E-34
NA = 6.02214E23
# 1H, 7Li, 17O
gammas = [2.6752E8, 1.0398E8, -3.6281E7]
multiplicities = [2, 4, 6]
B0 = np.array([0., 23.5])


def energy(mag_mom, gamma, field):
    """Energy caused by Zeeman splitting, Jmol-1"""
    return -mag_mom * hbar * gamma * field * NA


fig, ax = plt.subplots(figsize=(4.8, 4))
ax.set_xlim(B0)
ax.set_ylim(-1.5, 1.5)
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

props = dict(facecolor="w", edgecolor="none")
label_B0 = 21.3
label_shifts = iter([
    -0.05,
    0.05,
    -0.04,
    -0.05,
    0.05,
    0.04,
    -0.13,
    -0.1,
    0.07,
    -0.07,
    0.1,
    0.13,
])
for gamma, mult, color in zip(gammas, multiplicities, colors):
    for i in range(-mult + 1, mult + 1, 2):
        m = 0.5 * i
        energies = energy(m, gamma, B0)
        ax.plot(B0, energies, color=color)
        label = f"$m = <SIGN>\\nicefrac{{{abs(i)}}}{2}$"
        label = label.replace("<SIGN>", "+" if i > 0 else "-")
        label_y = energy(m, gamma, label_B0)
        try:
            label_y += next(label_shifts)
        except:
            pass
        ax.text(
            label_B0,
            label_y,
            label,
            color=color,
            fontsize="x-small",
            va="center",
        )

# Plot chosen magnetic fields (400, 600, 800MHz)
ymin, ymax = ax.get_ylim()

fields = [9.398, 14.095, 18.793]
labels = [f"\\qty{{{sfo}}}{{\\mega\\hertz}}" for sfo in (400, 600, 800)]
label_height = 0.02
label_color = "#c0c0c0"

trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
for field, label in zip(fields, labels):
    ax.axvline(field, color=label_color, ls="--", zorder=0)
    ax.text(
        field,
        label_height,
        label,
        transform=trans,
        color=label_color,
        bbox=props,
        ha="center",
        fontsize="small",
        clip_on=True,
        zorder=1,
    )

ax.set_xlabel("$B_0$ (\\unit{\\tesla})")
ax.set_ylabel("$E_m$ (\\unit{\\joule\\per\\mole})")

legend_labels = [
    "\\textsuperscript{1}H ($I = \\nicefrac{1}{2}$)",
    "\\textsuperscript{7}Li ($I = \\nicefrac{3}{2}$)",
    "\\textsuperscript{17}O ($I = \\nicefrac{5}{2}$)",
]
legend_lines = [
    lines.Line2D([0], [0], color=colors[0]),
    lines.Line2D([0], [0], color=colors[1]),
    lines.Line2D([0], [0], color=colors[2]),
]
fig.legend(legend_lines, legend_labels, loc='upper left', bbox_to_anchor=(0.12, 0.98))

fig.savefig("figures/energy_levels/energy_levels.pdf")
