# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Fri 07 Jul 2023 17:55:18 BST

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import colors
import nmrespy as ne
# mpl.use("tkAgg")
# mpl.rcParams["text.usetex"] = False

left = -0.1
right = 1.1
bottom = -0.1
top = 1.1
width = 0.5
spacing = 0.2


fig = plt.figure(figsize=(6, 2))
axs = []
fig, axs = plt.subplots(
    nrows=2,
    ncols=4,
    gridspec_kw={
        "left": 0.09,
        "right": 0.995,
        "bottom": 0.07,
        "top": 0.99,
        "wspace": 0.06,
        "hspace": 0.06,
    },
    figsize=(6, 3),
)
axs = list(axs.flatten())

params = np.array([[1., 0., 5., 0., 5., 4.]])
expinfo = ne.ExpInfo(
    dim=2,
    sw=(22., 15.),
    default_pts=(512, 256),
)
shifts_f1, shifts_f2 = expinfo.get_shifts()
fid_cos, fid_sin = expinfo.make_fid(params, indirect_modulation="amp")
phase_twist_double = ne.sig.ft(fid_cos).real
spec_f2_cos = ne.sig.ft(fid_cos, axes=[1])
spec_f2_cos.imag = 0.
double_absorption_cos = ne.sig.ft(spec_f2_cos, axes=[0]).real
spec_f2_sin = ne.sig.ft(fid_sin, axes=[1])
spec_f2_sin.imag = 0.
double_absorption_sin = ne.sig.ft(spec_f2_sin, axes=[0]).imag
freq_discrim = double_absorption_cos - double_absorption_sin

fid_pos, fid_neg = expinfo.make_fid(params, indirect_modulation="phase")
phase_twist_pos = ne.sig.ft(fid_pos)
phase_twist_neg = ne.sig.ft(fid_neg)
phase_twist_neg_flip = phase_twist_neg[::-1]
pos_plus_neg = phase_twist_pos + phase_twist_neg_flip

COLORS = mpl.rcParams["axes.prop_cycle"].by_key()["color"]
cmap = colors.LinearSegmentedColormap.from_list(
    "custom",
    [COLORS[1], "#f8f8f8", COLORS[0]],
)

base = 0.1
nlevs = 18
factor = 1.3
levels = [base * factor ** i for i in range(nlevs)]
levels = [-x for x in reversed(levels)] + levels
double_levels = [2 * x for x in levels]
quad_levels = [4 * x for x in levels]
kwargs = {
    "levels": levels,
    "cmap": cmap,
    "linewidths": 0.6,
}

for i, (ax, z) in enumerate(zip(
    axs,
    (
        phase_twist_double,
        double_absorption_cos,
        double_absorption_sin,
        freq_discrim,
        phase_twist_pos,
        phase_twist_neg,
        phase_twist_neg_flip,
        pos_plus_neg,
    ),
)):
    if i == 3:
        kwargs["levels"] = double_levels
    elif i == 7:
        kwargs["levels"] = quad_levels

    ax.contour(shifts_f2, shifts_f1, z, **kwargs)

    ax.set_xticks([0])
    ax.set_yticks([-5, 0, 5])

    if i == 0:
        ax.set_yticklabels([
            "$2f_{\\text{off}}^{(1)} - f_{\\vphantom{\\text{off}}}^{(1)}$",
            "$f_{\\text{off}}^{(1)}$",
            "$f_{\\vphantom{\\text{off}}}^{(1)}$",
        ])
        ax.set_xticklabels([])
    elif i == 4:
        ax.set_xticklabels(["$f_{\\vphantom{\\text{off}}}^{(2)}$"])
        ax.set_yticklabels([])
    else:
        ax.set_xticklabels([])
        ax.set_yticklabels([])

    ax.text(0.02, 0.91, f"\\textbf{{{chr(97 + i)}.}}", transform=ax.transAxes)
    if i == 3:
        ax.text(0.98, 0.93, "$\\times 0.5$", transform=ax.transAxes, fontsize=8, ha="right")
    if i == 7:
        ax.text(0.98, 0.93, "$\\times 0.25$", transform=ax.transAxes, fontsize=8, ha="right")

fig.savefig("figures/2d_lineshapes/2d_lineshapes.pdf")
