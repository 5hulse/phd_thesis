# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Wed 07 Jun 2023 18:28:54 BST

import nmrespy as ne
import matplotlib as mpl
from matplotlib.markers import MarkerStyle
import matplotlib.pyplot as plt
import numpy as np


sw = 500.e3
pts = 2 ** 15
nspins = 31
params = np.zeros((nspins, 4), dtype="float64")
params[:, 0] = 1.
freqs = np.linspace(-sw / 2, sw / 2, nspins + 2)
freqs = freqs[1:-1]
params[:, 2] = freqs
params[:, 3] = 1000.
estimator = ne.BBQChili.new_from_parameters(
    params,
    100.e-6,
    6.5e-6,
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
estimated = ne.sig.ft(estimated)
spectra.append(estimated.real)

padb, padt = (0.06, 0.04)
ylims = [(np.amin(x), np.amax(x)) for x in spectra]
heights = [ymax - ymin for (ymin, ymax) in ylims]
ylims = [(ymin - padb * h, ymax + padt * h) for (ymin, ymax), h in zip(ylims, heights)]
heights = [ymax - ymin for (ymin, ymax) in ylims]

left = 0.012
right = 0.7
bottom = 0.105
top = 0.99
fig, axs = plt.subplots(
    nrows=len(spectra),
    ncols=1,
    gridspec_kw={
        "left": left,
        "right": right,
        "bottom": bottom,
        "top": top,
        "hspace": 0.,
        "height_ratios": heights,
    },
    figsize=(6, 3),
)

axs = axs.tolist()

for i, (ax, spectrum, ylim) in enumerate(zip(axs, spectra, ylims)):
    if i == 0:
        ax.spines["bottom"].set_visible(False)
        ax.set_xticks([])
    elif i == 2:
        ax.spines["top"].set_visible(False)
    else:
        ax.spines["bottom"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.set_xticks([])

    ax.plot(shifts, spectrum, color="k")
    ax.set_xlim(shifts[0], shifts[-1])
    ax.set_ylim(ylim)
    ax.set_yticks([])

    ax.text(0.005, 0.97, f"\\textbf{{{chr(97 + i)}.}}", transform=ax.transAxes, va="top")


xticks = list(range(250, -300, -50))
axs[2].set_xticks([x * 1000 for x in xticks])
axs[2].set_xticklabels([str(i) for i in xticks])
axs[2].set_xlabel("\\unit{\\kilo\\hertz}")


# Plot phases
params = estimator.get_params()
phi = params[:, 1]
phi_quad = phi.copy()

for i in range(16, -1, -1):
    while True:
        phase = phi_quad[i]
        if phase > phi_quad[i + 1]:
            break
        else:
            phi_quad[i] += 2 * np.pi

for i in range(18, 31):
    while True:
        phase = phi_quad[i]
        if phase > phi_quad[i - 1]:
            break
        else:
            phi_quad[i] += 2 * np.pi

freq = params[:, 2]

axs.append(
    fig.add_axes(
        [right + 0.02, bottom, 1 - left - right - 0.02, top - bottom]
    )
)
colors = mpl.rcParams["axes.prop_cycle"].by_key()["color"]
red, blue = colors[:2]
kwargs = dict(
    s=8,
    edgecolor="k",
    lw=0.4,
)
matching = np.where(phi - phi_quad == 0)[0]
left, right = matching[0], matching[-1]
# axs[-1].scatter(freq[:left], phi[:left], color=red, **kwargs)
# axs[-1].scatter(freq[right:], phi[right:], color=red, **kwargs)
# axs[-1].scatter(freq[:left], phi_quad[:left], color=blue, **kwargs)
# axs[-1].scatter(freq[right:], phi_quad[right:], color=blue, **kwargs)
# axs[-1].scatter(
#     freq[left:right + 1], phi[left:right + 1], color=red, marker=MarkerStyle("o", fillstyle="right"),
#     **kwargs,
# )
# axs[-1].scatter(
#     freq[left:right + 1], phi[left:right + 1], color=blue, marker=MarkerStyle("o", fillstyle="left"),
#     **kwargs,
# )

axs[-1].axhline(np.pi, color="k", ls=":")
axs[-1].axhline(-np.pi, color="k", ls=":")

phi2, phi1, phi0 = np.polyfit(freq, phi_quad, 2)
xs = np.linspace(-sw / 2, sw / 2, 100)
term1 = xs * (0.5 * estimator.pulse_length + estimator.prescan_delay)
term2 = 0.5 * (xs ** 2) * estimator.pulse_length / estimator.pulse_bandwidth
axs[-1].plot(xs, 2 * np.pi * (term1 + term2))
axs[-1].plot(xs, phi2 * xs ** 2 + phi1 * xs + phi0, color="k", zorder=-1)

fig.savefig("figures/chirp_phase_vs_estimation/chirp_phase_vs_estimation.pdf")
