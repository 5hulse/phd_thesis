# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Wed 21 Jun 2023 11:16:37 BST

from pathlib import Path
import nmrespy as ne
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


data_root = Path("/home/simon/Documents/DPhil/data/Single_Chirp_Experimental/")
data_dirs = [data_root / f"{i}/pdata/1001" for i in range(1, 4)]

spectra = []
for direc in data_dirs:
    spectra.append(np.fromfile(direc / "1r", "<i4"))

datapath = Path("~/Documents/DPhil/data/Single_Chirp_Experimental/1/pdata/1001").expanduser()
estimator = ne.BBQChili.new_bruker(datapath, 100.e-6, 400.e3)
estimator._default_pts = spectra[0].shape
estimator._pts = spectra[0].shape

top = 0.99
bottom = 0.065
left = 0.006
right = 0.985
sep = 0.01
seps = [sep, 0, sep, 0]
magn = 30.
limits = [
    (-3.8e7, 3.8e7),
    (-2.7e6, 3.8e7),
    (-4e5, 6.e5),
    (-2.7e6, 3.8e7),
    (-4e5, 6.e5),
]
heights = []
for i, lim in enumerate(limits):
    height = lim[1] - lim[0]
    if i in [2, 4]:
        height *= magn
    heights.append(height)

h_sum = sum(heights)
h_ratios = [h / h_sum * (top - bottom - sum(seps)) for h in heights]

fig = plt.figure(figsize=(6., 4.))
axs = []
for i in range(5):
    print([
        left,
        top - sum(h_ratios[:i + 1]) - sum(seps[:i]),
        right - left,
        h_ratios[i],
    ])
    axs.append(
        fig.add_axes(
            [
                left,
                top - sum(h_ratios[:i + 1]) - sum(seps[:i]),
                right - left,
                h_ratios[i],
            ]
        )
    )

shifts, = estimator.get_shifts()
axs[0].plot(shifts, spectra[0], color="k")
axs[1].plot(shifts, spectra[1], color="k")
axs[2].plot(shifts, spectra[1], color="k")
axs[3].plot(shifts, spectra[2], color="k")
axs[4].plot(shifts, spectra[2], color="k")

axs[1].spines["bottom"].set_visible(False)
axs[2].spines["top"].set_linestyle((0, (1, 2)))
axs[3].spines["bottom"].set_visible(False)
axs[4].spines["top"].set_ls((0, (1, 2)))

for ax in [axs[2], axs[4]]:
    ax.text(0.996, 0.95, f"$\\times {int(magn)}$", ha="right", va="top", transform=ax.transAxes, fontsize=7)

for i, (ax, lim) in enumerate(zip(axs, limits)):
    ax.set_xlim(shifts[0], shifts[-1])
    ax.set_ylim(lim)
    if i != 4:
        ax.set_xticks([])
    ax.set_yticks([])

axs[4].set_xlabel("\\unit{\\kilo\\hertz}")
axs[4].set_xticks([-i * 50000 for i in range(7)])
axs[4].set_xticklabels([f"{-i * 50}" for i in range(7)])

for i, ax in enumerate((axs[0], axs[1], axs[3])):
    ax.text(48e3, 3.e7, f"\\textbf{{{chr(97 + i)}.}}")

fig.savefig("figures/chirp_phase_vs_estimation_real_data/chirp_phase_vs_estimation_real_data.pdf")
