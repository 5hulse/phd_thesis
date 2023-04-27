# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Tue 18 Apr 2023 08:23:21 BST

# from numpy.linalg import norm
# from numpy.random import normal

import copy
from pathlib import Path

from matplotlib import cm, lines, pyplot as plt
import numpy as np
import nmrespy as ne

LEFT = 0.01
BOTTOM = 0.1
TOP = 0.98
RIGHT = 0.98
HSPACE = 0.07
VSPACE = 0.11

PTS = 64

width = (RIGHT - LEFT - HSPACE) / 2
height = (TOP - BOTTOM - VSPACE) / 2

fig = plt.figure(figsize=(6, 3))
axs = []
rectangles = [
    [LEFT, (BOTTOM + TOP + VSPACE) / 2, width, height],
    [LEFT, BOTTOM, width, height],
    [(LEFT + RIGHT + HSPACE) / 2, BOTTOM, width, TOP - BOTTOM],
]
for rectangle in rectangles:
    axs.append(fig.add_axes(rectangle))


if (file := Path("figures/optimisation_visualisation/estimator")).with_suffix(".pkl").is_file():
    estimator = ne.Estimator1D.from_pickle(file)
else:
    estimator = ne.Estimator1D.new_from_parameters(
        params=np.array([[1., 0., 1., .2]]),
        pts=PTS,
        sw=5.2,
        offset=0.,
        snr=10.,
    )
    for hessian in ("gauss-newton", "exact"):
        estimator.estimate(
            initial_guess=np.array([[1., 0., 1.1, .8]]),
            output_mode=1,
            save_trajectory=True,
            mode="fd",
            negative_amps="ignore",
            hessian=hessian,
        )
    estimator.to_pickle(file)

fid = estimator.data
tp, = estimator.get_timepoints()
shifts, = estimator.get_shifts(pts=(2 * PTS,))
trajectories = [res.trajectory for res in estimator._results]

freq_range = (0.95, 1.12)
damp_range = (0., 0.9)

axs[0].scatter(tp, np.real(fid), linewidths=0, color="k", zorder=100, s=8)
spec = copy.deepcopy(fid)
spec = ne.sig.exp_apodisation(spec, k=1.)
spec = ne.sig.zf(spec)
spec[0] *= 0.5
axs[1].scatter(shifts, ne.sig.ft(spec).real, color="k", linewidths=0, s=8, zorder=100)

for trajectory in trajectories:
    fs = [x[0, 2] for x in trajectory]
    ds = [x[0, 3] for x in trajectory]
    axs[2].plot(fs, ds, lw=0.8)

    for x, marker, ls in zip(
        (trajectory[0], trajectory[-1]),
        ("o", "s"),
        ("-", ":"),
    ):
        model = estimator.make_fid(x)
        axs[0].plot(tp, model.real, color="#808080", zorder=1000, ls=ls)
        model = ne.sig.exp_apodisation(model, k=1.)
        model = ne.sig.zf(model)
        model[0] *= 0.5
        spec = ne.sig.ft(model).real
        axs[1].plot(shifts, spec, color="#808080", ls=ls)
        axs[2].scatter(x[0, 2], x[0, 3], edgecolor="k", facecolor="#808080", s=6, zorder=100, linewidths=0.5, marker=marker)

samples = 101
cost_function_surface = np.zeros((samples, samples))
freqs = np.linspace(freq_range[0], freq_range[1], samples)
damps = np.linspace(damp_range[0], damp_range[1], samples)
params = copy.deepcopy(trajectory[0])
for i, f in enumerate(freqs):
    for j, d in enumerate(damps):
        params[0, 2] = f
        params[0, 3] = d
        model = estimator.make_fid(params)
        cost_function_surface[i, j] = np.linalg.norm(fid / np.linalg.norm(fid) - model / np.linalg.norm(fid)) ** 2

ff, dd = np.meshgrid(freqs, damps, indexing="ij")
base = 0.05
factor = 1.35
levels = [base * factor ** i for i in range(16)]
colors = cm.Greys(np.linspace(0.3, 1, len(levels)))
contour = axs[2].contour(ff, dd, cost_function_surface, levels=levels, colors=colors, zorder=-1)
lab = axs[2].clabel(contour, fontsize=4, fmt="%.2f")

axs[0].set_yticks([])
axs[1].set_yticks([])
axs[0].set_xlabel("$t (\\unit{\\second})$")
axs[0].set_ylim(bottom=1.05 * axs[0].get_ylim()[0])
axs[0].set_ylim(top=1.2 * axs[0].get_ylim()[1])
axs[1].set_xlim(1.5, 0.5)
axs[1].set_xlabel("$f (\\unit{\\hertz})$")
axs[2].set_xlabel("$f^{(1)} (\\unit{\\hertz})$")
axs[2].set_ylabel("$\\eta^{(1)} (\\unit{\\per\\second})$", labelpad=2)

for i, (x, y) in enumerate(zip((0.015, 0.015, 0.535), (0.95, 0.445, 0.95))):
    fig.text(x, y, f"\\textbf{{{chr(97 + i)}.}}")


fig.savefig("figures/optimisation_visualisation/optimisation_visualisation.pdf")
