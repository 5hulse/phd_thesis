# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Thu 15 Jun 2023 19:43:54 BST

from pathlib import Path
import nmrespy as ne
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


folder = Path("/home/simon/Documents/DPhil/data/Single_Chirp_Experimental/1/pdata/1001")
real = np.fromfile(folder / "1r", "<i4")[::-1]
imag = np.fromfile(folder / "1i", "<i4")[::-1]
data = ne.sig.ift(real + 1j * imag)

datapath = Path("~/Documents/DPhil/data/Single_Chirp_Experimental/1/pdata/1001").expanduser()
estimator = ne.BBQChili.new_bruker(datapath, 100.e-6, 400.e3)
estimator._data = data
estimator._default_pts = data.shape
estimator._pts = data.shape
estimator.estimate(initial_guess=13, mpm_trim=4096, nlp_trim=16384)
estimator.prescan_delay = None

p0, p1 = (0.921,), (-1.552,)
init_spec = estimator.spectrum.real
phased_spec = ne.sig.phase(estimator.quadratic_phase(), p0=p0, p1=p1).real
bbqchili_fid = estimator.back_extrapolate()
bbqchili_fid[0] *= 0.5
bbqchili_spec = ne.sig.phase(ne.sig.ft(bbqchili_fid), p0=p0, p1=p1).real

top = 0.99
bottom = 0.05
left = 0.01
right = 0.99
sep = 0.02
seps = [sep, 0, sep, 0]
magn = 30.
limits = [
    (-3.8e7, 3.8e7),
    (-2.e6, 3.8e7),
    (-3.5e5, 1.e6),
    (-2.e6, 3.8e7),
    (-2.e5, 1.e6),
]
heights = []
for i, lim in enumerate(limits):
    height = lim[1] - lim[0]
    if i in [2, 4]:
        height *= magn
    heights.append(height)

h_sum = sum(heights)
h_ratios = [h / h_sum * (top - bottom - sum(seps)) for h in heights]

fig = plt.figure(figsize=(6., 6.))
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
axs[0].plot(shifts, init_spec)
axs[1].plot(shifts, phased_spec)
axs[2].plot(shifts, phased_spec)
axs[3].plot(shifts, bbqchili_spec)
axs[4].plot(shifts, bbqchili_spec)

axs[1].spines["bottom"].set_visible(False)
axs[2].spines["top"].set_ls(":")
axs[3].spines["bottom"].set_visible(False)
axs[4].spines["top"].set_ls(":")

for ax in [axs[2], axs[4]]:
    ax.text(0.02, 0.98, f"$\\times {int(magn)}$", va="top", transform=ax.transAxes)

for i, (ax, lim) in enumerate(zip(axs, limits)):
    ax.set_xlim(shifts[0], shifts[-1])
    ax.set_ylim(lim)
    if i != 4:
        ax.set_xticks([])
    ax.set_yticks([])

fig.savefig("figures/chirp_phase_vs_estimation_real_data/chirp_phase_vs_estimation_real_data.pdf")

