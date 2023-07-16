# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Fri 14 Jul 2023 20:30:03 BST

import nmrespy as ne
import matplotlib as mpl
from matplotlib.markers import MarkerStyle
import matplotlib.pyplot as plt
import numpy as np


sw = 400.e3
pts = 2 ** 14
nspins = 30
params = np.zeros((nspins, 4), dtype="float64")
params[:, 0] = 1.
params[:, 1] = 0.
freqs = np.linspace(-sw / 2, sw / 2, nspins + 2)
freqs = freqs[1:-1]
params[:, 2] = freqs
params[:, 3] = 1000.
estimator = ne.BBQChili.new_from_parameters(
    params,
    100.e-6,
    0.,
    400.e3,
    pts,
    sw,
    0.,
    snr=25.,
)
estimator.estimate(mpm_trim=2048, nlp_trim=4096)
lines = []

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
right = 0.72
bottom = 0.095
top = 0.99
fig, axs = plt.subplots(
    nrows=len(spectra),
    ncols=1,
    gridspec_kw={
        "left": left,
        "right": right,
        "bottom": bottom,
        "top": top,
        "hspace": 0.04,
        "height_ratios": heights,
    },
    figsize=(6, 3),
)

axs = axs.tolist()

for i, (ax, spectrum, ylim) in enumerate(zip(axs, spectra, ylims)):
    if i != 2:
        ax.set_xticks([])

    ax.plot(shifts, spectrum, color="k")
    ax.set_xlim(shifts[0], shifts[-1])
    ax.set_ylim(ylim)
    ax.set_yticks([])

    ax.text(199e3, ylim[1] - 80., f"\\textbf{{{chr(97 + i)}.}}")

for osc in estimator.get_params():
    fid = estimator.make_fid(np.expand_dims(osc, axis=0))
    fid[0] *= 0.5
    line = ne.sig.ft(fid).real / 1.1
    axs[0].plot(shifts - 3000, line, lw=0.8, zorder=-1)

xticks = list(range(200, -250, -50))
axs[2].set_xticks([x * 1000 for x in xticks])
axs[2].set_xticklabels([f"${str(i)}$" for i in xticks])
axs[2].set_xlabel("\\unit{\\kilo\\hertz}")

# Plot phases
params = estimator.get_params()
phi = params[:, 1]
phi_quad = phi.copy()

pivot = 6

for i in range(pivot - 1, -1, -1):
    while True:
        phase = phi_quad[i]
        if phase < phi_quad[i + 1]:
            break
        else:
            phi_quad[i] -= 2 * np.pi

for i in range(pivot + 1, nspins):
    while True:
        phase = phi_quad[i]
        if phase < phi_quad[i - 1]:
            break
        else:
            phi_quad[i] -= 2 * np.pi

freq = params[:, 2]

axs.append(
    fig.add_axes(
        [right + 0.04, bottom, 1 - left - right - 0.04, top - bottom]
    )
)
colors = mpl.rcParams["axes.prop_cycle"].by_key()["color"]
red, blue = colors[:2]
kwargs = dict(
    s=10,
    edgecolor="k",
    lw=0.4,
)
matching = np.where(phi - phi_quad == 0)[0]
left, right = matching[0], matching[-1]
axs[-1].scatter(freq[:left], phi[:left], color=red, **kwargs)
axs[-1].scatter(freq[right:], phi[right:], color=red, **kwargs)
axs[-1].scatter(freq[:left], phi_quad[:left], color=blue, **kwargs)
axs[-1].scatter(freq[right:], phi_quad[right:], color=blue, **kwargs)
axs[-1].scatter(
    freq[left:right + 1], phi[left:right + 1], color=red, marker=MarkerStyle("o", fillstyle="right"),
    **kwargs,
)
axs[-1].scatter(
    freq[left:right + 1], phi[left:right + 1], color=blue, marker=MarkerStyle("o", fillstyle="left"),
    **kwargs,
)

axs[-1].axhline(np.pi, color="k", ls=(0, (1, 0.8)), zorder=-1)
axs[-1].axhline(-np.pi, color="k", ls=(0, (1, 0.8)), zorder=-1)

phi2, phi1, phi0 = np.polyfit(freq, phi_quad, 2)
xs = np.linspace(-sw / 2, sw / 2, 100)
axs[-1].plot(xs, phi2 * xs ** 2 + phi1 * xs + phi0, color="k", zorder=-1)
axs[-1].set_xticks([-2e5, 0, 2e5])
axs[-1].set_xticklabels(["$-200$", "$0$", "$200$"])
axs[-1].set_xlabel("$f$ (\\unit{\\kilo\\hertz})", labelpad=-1.5)
axs[-1].set_yticks([-np.pi, 0, np.pi])
axs[-1].set_yticklabels(["$-\\pi$", "$0$", "$\\pi$"], va="center")
axs[-1].set_ylabel("$\\phi$ (\\unit{\\radian})", labelpad=-10)
axs[-1].text(
    0.01, 0.99, "\\textbf{{d.}}", transform=axs[-1].transAxes, va="top",
    bbox={"facecolor": "w", "edgecolor": "none", "pad": 0.3},
)
axs[-1].set_xlim(reversed(axs[-1].get_xlim()))
axs[-1].set_ylim(top=3 * np.pi)

text_x = 0.18
text_y = 0.15
axs[-1].text(text_x, text_y, f"$\\phi=$", va="top", fontsize=6, transform=axs[-1].transAxes)
axs[-1].text(
    text_x + 0.1, text_y + 0.01,
    f"$\\num{{{phi2:.2e}}}f^2$\n$\\num{{{phi1:.2e}}}f$\n$\\num{{{phi0:.2f}}}$",
    va="top", fontsize=6, transform=axs[-1].transAxes,
)

rectangle = mpl.patches.Rectangle(
    (text_x - 0.015, text_y - 0.11),
    0.53,
    0.13,
    transform=axs[-1].transAxes,
    edgecolor="k",
    lw=0.6,
    facecolor="w",
)
axs[-1].add_patch(rectangle)

fig.savefig("figures/chirp_phase_vs_estimation/chirp_phase_vs_estimation.pdf")
