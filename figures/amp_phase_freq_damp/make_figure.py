# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Mon 15 May 2023 15:17:57 BST

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import nmrespy as ne


amps = [1, 0.5, 0.25]
phis = [0., np.pi / 4, np.pi]
init_freq = 0.1
freqs = [init_freq, init_freq / 2, 0.]
damps = [0.05, 0.1, 0.2]
params = [amps, phis, freqs, damps]
init_params = np.array([[amps[0], phis[0], freqs[0], damps[0]]])

expinfo = ne.ExpInfo(dim=1, sw=5., default_pts=256)
tp, = expinfo.get_timepoints()
shifts, = expinfo.get_shifts(pts=4096)

fig, axs = plt.subplots(
    nrows=4,
    ncols=2,
    gridspec_kw=dict(
        left=0.01,
        right=0.99,
        bottom=0.045,
        top=0.99,
        wspace=0.02,
        hspace=0.05,
    ),
    figsize=(5, 4),
)
for ax_row in axs:
    for ax in ax_row:
        ax.set_xticks([])
        ax.set_yticks([])

for i, variables in enumerate(params):
    ax_row = axs[i]
    for var in variables:
        p = np.copy(init_params)
        p[0, i] = var
        fid = expinfo.make_fid(p)
        spec = expinfo.make_fid(p, pts=4096)
        spec = ne.sig.ft(spec)
        line, = ax_row[0].plot(tp, fid.real)
        color = line.get_color()
        ax_row[0].plot(tp, fid.imag, color=color, alpha=0.4)
        line, = ax_row[1].plot(shifts, spec.real)
        color = line.get_color()
        ax_row[1].plot(shifts, spec.imag, color=color, alpha=0.4)


ylim0 = axs[1, 0].set_ylim()[1]
for ax in axs[:, 0]:
    ax.set_xlim(-tp[20], tp[-1])
    ax.set_ylim(-ylim0, ylim0)

ylim1 = axs[0, 1].get_ylim()[1]
for ax in axs[:, 1]:
    ax.axvline(0, color="k", ls=":", lw=0.7)
    ax.set_xlim(0.27, -0.07)
    ax.set_ylim(-ylim1, ylim1)

axs[3, 0].set_xlabel("$t$")
axs[3, 1].set_xlabel("$f$", labelpad=-8)
axs[3, 1].set_xticks([0])

for i in range(4):
    for j in range(2):
        letter = chr(97 + i)
        axs[i, j].text(0.01, 0.85, f"\\textbf{{{letter}{j + 1}.}}", transform=axs[i, j].transAxes)

fig.savefig("figures/amp_phase_freq_damp/amp_phase_freq_damp.pdf")
