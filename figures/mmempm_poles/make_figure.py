# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Thu 03 Aug 2023 18:21:35 BST

# I made a temporary branch in NMR-EsPy to produce this figure.
# After creating the new branch, I made the `Z` and `groupings`
# variables in `MatrixPencil._mpm_2d` attributes, so they could be accessed.
# The result was also return unsorted

import re
import numpy as np
import nmrespy as ne
import matplotlib as mpl
from matplotlib import colors
import matplotlib.pyplot as plt


def ud_groupings(groupings, order):
    dic = {i: j for i, j in enumerate(order)}
    new_groupings = []
    for group in groupings:
        new_group = []
        for idx in group:
            new_group.append(dic[idx])
        new_groupings.append(new_group)
    return new_groupings


cols = mpl.rcParams["axes.prop_cycle"].by_key()["color"]
red = cols[0]
cols.pop(4)

fig, axs = plt.subplots(
    nrows=2,
    ncols=2,
    gridspec_kw={
        "left": 0.08,
        "right": 0.995,
        "bottom": 0.09,
        "width_ratios": [2, 1],
        "hspace": 0.25,
    },
    figsize=(5, 3.2),
)

expinfo = ne.ExpInfo(dim=2, sw=(200., 200.), default_pts=(128, 128))
params_match = np.array(
    [
        [1., 0., 50., 50., 5., 5.],
        [1., 0., 50., -50., 5., 5.],
        [1., 0., -50., 20., 5., 5.],
        [1., 0., -50., -70., 5., 5.],
        [1., 0., 20., -20., 5., 5.],
    ]
)
params_unmatch = np.array(
    [
        [1., 0., 70., 50., 5., 5.],
        [1., 0., 50., -50., 5., 5.],
        [1., 0., -50., 20., 5., 5.],
        [1., 0., -70., -70., 5., 5.],
        [1., 0., 20., -20., 5., 5.],
    ]
)

fid_match = fid_unmatch = expinfo.make_fid(params=params_unmatch)
f1, f2 = expinfo.get_shifts()
cmap = colors.LinearSegmentedColormap.from_list("wr_custom", ["#ffffff", "#000000"], N=128)
base = 60
nlevs = 10
factor = 2.
levels = [base * factor ** i for i in range(nlevs)]
norm = colors.Normalize(vmin=0, vmax=1)

# HACK: not sure why the groupings are mismatched (see commented line below)
gs = ([[0], [1], [2], [3], [4]], [[2], [0, 3], [1, 4]])
for i, (params, groupings) in enumerate(zip((params_unmatch, params_match), gs)):
    fid = expinfo.make_fid(params=params)
    mpm = ne.mpm.MatrixPencil(expinfo, fid, oscillators=5)
    order = mpm.order
    Z = (np.abs(mpm.Z) > 1.e-10).astype(int)
    Z = Z[order]
    Z = Z[:, order]
    # groupings = ud_groupings(mpm.groupings, order)
    spec = np.abs(ne.sig.ft(fid)).real
    axs[i, 0].contour(f2, f1, spec, levels=levels, colors="k", linewidths=0.5)
    axs[i, 1].matshow(Z, cmap=cmap, norm=norm)
    params = mpm.params

    # F2 freq positions
    for j, para in enumerate(params):
        freq2 = para[3]
        axs[i, 0].axvline(freq2, color=cols[j], ls=(0, (5, 2)), lw=0.8)

    # F1 freq positions
    for group in groupings:
        if len(group) == 1:
            idx = group[0]
            freq1 = params[idx, 2]
            c = cols[idx]
            axs[i, 0].axhline(freq1, color=c, ls=(0, (5, 2)), lw=0.8)
        else:
            for incr, idx in enumerate(group):
                freq1 = params[idx, 2]
                c = cols[idx]
                axs[i, 0].axhline(freq1, color=c, ls=(7 * incr, (5, 9)), lw=0.8)

    axs[i, 0].set_xlim(100, -100)
    axs[i, 0].set_ylim(100, -100)
    axs[i, 1].tick_params(top=False, labeltop=False, bottom=True, labelbottom=True)
    axs[i, 1].set_xticks(list(range(5)))
    axs[i, 1].set_xticklabels([f"${z}$" for z in range(1, 6)])
    axs[i, 1].set_yticks(list(range(5)))
    axs[i, 1].set_yticklabels([f"${z}$" for z in range(1, 6)])
    for j, tick in enumerate(axs[i, 1].get_xticklabels()):
        tick.set_color(cols[j])
    for j, tick in enumerate(axs[i, 1].get_yticklabels()):
        tick.set_color(cols[j])
    axs[i, 0].set_xlabel("$F^{(2)}$ (Hz)")
    axs[i, 0].set_ylabel("$F^{(1)}$ (Hz)")
    axs[i, 0].text(0.01, 0.98, f"\\textbf{{{chr(97 + i)}1.}}", transform=axs[i, 0].transAxes, va="top")
    axs[i, 1].text(0.02, 0.98, f"\\textbf{{{chr(97 + i)}2.}}", transform=axs[i, 1].transAxes, va="top", color="w")



fig.savefig("figures/mmempm_poles/mmempm_poles.pdf")
