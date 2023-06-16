# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Wed 14 Jun 2023 21:19:11 BST

import pickle
import nmrespy as ne
import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, "figures/mdl")
from curlybrace import curlyBrace

colors = mpl.rcParams["axes.prop_cycle"].by_key()["color"]

def get_info(data):
    data /= np.linalg.norm(data)
    N = data.size
    L = int(np.floor(N / 3))
    HY = sp.linalg.hankel(data[: N - L], data[N - L - 1 :])
    _, sigma, _ = np.linalg.svd(HY)
    pdf = np.zeros(L)
    mdl = np.zeros(L)
    penalty = np.zeros(L)
    for k in range(L):
        penalty[k] = k * np.log(N) * (2 * L - k) / 2
        pdf[k] = (
            N * np.einsum("i->", np.log(sigma[k:])) -
            N * (L - k) * np.log(np.einsum("i->", sigma[k:]) / (L - k))
        )
        mdl[k] = -pdf[k] + penalty[k]

    return sigma, pdf, penalty, mdl


MAKE_FIDS = False

left = 0.01
right = 0.99
top = 0.95
bottom = 0.06
wspaces = [0.06, 0.045]
w = ((right - left) - sum(wspaces)) / 3
h = top - bottom
fig = plt.figure(figsize=(6., 3.))
axs = []
for i in range(3):
    l = left + (i * w) + sum(wspaces[:i])
    axs.append(fig.add_axes([l, bottom, w, h]))

plot_kwargs = {"lw": 1.}
scatter_kwargs = dict(
    edgecolor="k",
    lw=0.4,
    s=8,
)

expinfo = ne.ExpInfo(1, sw=1., default_pts=256)
quartet_center = -0.25
triplet_center = 0.3
splitting = 0.03
damping = 0.02
max_order = 14

q_min = quartet_center - 1.5 * splitting
q_max = quartet_center + 1.5 * splitting
t_min = triplet_center - splitting
t_max = triplet_center + splitting
q_freqs = np.linspace(q_min, q_max, 4)
t_freqs = np.linspace(t_min, t_max, 3)
q_amps = np.array([1, 3, 3, 1], dtype="float32")
t_amps = np.array([2, 4, 2], dtype="float32")

params = np.zeros((7, 4), dtype="float32")
params[:4, 0] = q_amps
params[4:, 0] = t_amps
params[:4, 2] = q_freqs
params[4:, 2] = t_freqs
params[:, 3] = damping

snrs = [7., 12., 20.]

shifts = expinfo.get_shifts()[0]
vshift = 0.
ns = np.arange(max_order)

for snr in snrs:
    if MAKE_FIDS:
        fid = expinfo.make_fid(params, snr=snr)
        with open(f"figures/mdl/snr_{snr}.pkl", "wb") as fh:
            pickle.dump(fid, fh)
    else:
        with open(f"figures/mdl/snr_{snr}.pkl", "rb") as fh:
            fid = pickle.load(fh)

    svs, pdf, penalty, mdl = [x [:max_order] for x in get_info(fid)]
    color = axs[1].plot(ns, svs, zorder=0, **plot_kwargs)[0].get_color()
    axs[1].scatter(ns, svs, color=color, **scatter_kwargs)
    axs[2].plot(ns, -pdf, zorder=0, color=color, **plot_kwargs, ls=(0, (1, 0.5)))
    axs[2].scatter(ns, -pdf, color=color, **scatter_kwargs, marker="s")
    axs[2].plot(ns, mdl, zorder=0, color=color, **plot_kwargs)
    mdl_min_idx = np.argmin(mdl)
    mdl_min = mdl[mdl_min_idx]
    mdl = np.delete(mdl, mdl_min_idx)
    ns_cp = ns.copy()
    ns_cp = np.delete(ns_cp, mdl_min_idx)
    axs[2].scatter(ns_cp, mdl, color=color, **scatter_kwargs)
    axs[2].scatter(mdl_min_idx, mdl_min, color=color, s=25, edgecolor="k", lw=0.4, marker="*", zorder=100)
    fid[0] *= 0.5
    spec = ne.sig.ft(fid).real
    specmin = np.amax(spec)
    spec += vshift
    vshift = np.amax(spec)
    axs[0].plot(shifts, spec)

for i, x in enumerate(params[:, 2], start=1):
    idx = expinfo.convert([float(x)], "hz->idx")[0]
    y = np.amax(spec[idx - 2 : idx + 2]) + 0.25
    if i == 1:
        shift = -0.01
    elif i == 2:
        shift = -0.01
    elif i == 3:
        shift = 0.01
    elif i == 4:
        shift = 0.01
    elif i == 5:
        shift = -0.01
    elif i == 7:
        shift = 0.01
    else:
        shift = 0
    axs[0].text(x + shift, y, f"({i})", color=color, ha="center", fontsize=6)

axs[1].text(0.4, 2.73, "(6)", color=color, fontsize=6)
curlyBrace(fig, axs[1], (1.3, 2.14), (2.3, 1.91), color=color, solid_capstyle="round")
axs[1].text(2.3, 2.05, "(2), (3)", color=color, fontsize=6)
curlyBrace(fig, axs[1], (3.3, 1.435), (4.3, 1.27), color=color, solid_capstyle="round")
axs[1].text(4.25, 1.4, "(5), (7)", color=color, fontsize=6)
curlyBrace(fig, axs[1], (5.7, 0.54), (4.7, 0.67), color=color, solid_capstyle="round")
axs[1].text(3., 0.5, "(1), (4)", color=color, fontsize=6)


axs[2].plot(ns, penalty, zorder=0, color="#808080", **plot_kwargs)
axs[2].scatter(ns, penalty, color="#808080", **scatter_kwargs)

axs[2].axvline(7., color="k", zorder=-1, lw=0.5)

axs[0].set_ylim(top=np.amax(spec) + 0.7)
axs[0].set_xticks([])
axs[0].set_yticks([])
axs[0].set_xlim(shifts[0], shifts[-1])
axs[0].set_xlabel("Hz", labelpad=1)

for label, ax in zip(("$r$", "$k$"), (axs[1], axs[2])):
    ax.set_xlabel(label, labelpad=-1, va="bottom")
    ax.set_xlim(left=-1.5)
    ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(1))
    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(5))
    ax.xaxis.get_minor_ticks()[1].draw = lambda *args:None

xlim1 = axs[1].get_xlim()
axs[1].axvspan(xlim1[0] , 6.5, facecolor=colors[4], ls=":", alpha=0.2, zorder=-1)
axs[1].axvspan(6.5, xlim1[1], facecolor=colors[5], ls=":", alpha=0.2, zorder=-1)
axs[1].set_xlim(xlim1)
axs[1].text(0.1, 0.1, "signal", color=colors[4], fontsize=8, transform=axs[1].transAxes)
axs[1].text(0.9, 0.9, "noise", color=colors[5], fontsize=8, transform=axs[1].transAxes, ha="right", va="top")
axs[1].set_ylabel("$\\symbf{\\sigma}[r]$", labelpad=2)
axs[1].set_yticks([0.] + list(axs[1].get_yticks())[:-1])

axs[2].ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)

for i, ax in enumerate(axs):
    ax.text(0.01, 0.95, f"\\textbf{{{chr(97 + i)}.}}", transform=ax.transAxes)
fig.savefig("figures/mdl/mdl.pdf")
