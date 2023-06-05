# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Sun 04 Jun 2023 19:35:42 BST

import nmrespy as ne
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


sw = 500.e3
pts = 2 ** 15
nspins =30
params = np.zeros((nspins, 4), dtype="float64")
params[:, 0] = 1.
freqs = np.linspace(-sw / 2, sw / 2, nspins + 2)
freqs = freqs[1:-1]
params[:, 2] = freqs
params[:, 3] = 1000.
estimator = ne.BBQChili.new_from_parameters(
    params,
    100.e-6,
    0.,
    500.e3,
    pts,
    sw,
    0.,
    snr=50.,
)
estimator.estimate(mpm_trim=2048, nlp_trim=4096)


shifts, = estimator.get_shifts()
spectra = []
spectra.append(estimator.spectrum.real)
phasecorr = estimator.quadratic_phase()
spectra.append(phasecorr.real)
estimated = estimator.back_extrapolate()
estimated[0] *= 0.5
spectra.append(ne.sig.ft(estimated).real)

padb, padt = (0.06, 0.02)
ylims = [(np.amin(x), np.amax(x)) for x in spectra]
heights = [ymax - ymin for (ymin, ymax) in ylims]
ylims = [(ymin - padb * h, ymax + padt * h) for (ymin, ymax), h in zip(ylims, heights)]
heights = [ymax - ymin for (ymin, ymax) in ylims]

fig, axs = plt.subplots(
    nrows=3,
    ncols=1,
    gridspec_kw={
        "left": 0.012,
        "right": 0.985,
        "bottom": 0.105,
        "top": 0.99,
        "hspace": 0.,
        "height_ratios": heights,
    },
    figsize=(6, 2.5),
)

for i, (ax, spectrum, ylim) in enumerate(zip(axs, spectra, ylims)):
    if i == 0:
        ax.spines["bottom"].set_visible(False)
        ax.set_xticks([])
    elif i == 1:
        ax.spines["bottom"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.set_xticks([])
    elif i == 2:
        ax.spines["top"].set_visible(False)

    ax.plot(shifts, spectrum, color="k")
    ax.set_xlim(shifts[0], shifts[-1])
    ax.set_ylim(ylim)
    ax.set_yticks([])

    ax.text(0.005, 0.97, f"\\textbf{{{chr(97 + i)}.}}", transform=ax.transAxes, va="top")


xticks = list(range(250, -300, -50))
axs[2].set_xticks([x * 1000 for x in xticks])
axs[2].set_xticklabels([str(i) for i in xticks])
axs[2].set_xlabel("\\unit{\\kilo\\hertz}")

fig.savefig("figures/chirp_phase_vs_estimation/chirp_phase_vs_estimation.pdf")
