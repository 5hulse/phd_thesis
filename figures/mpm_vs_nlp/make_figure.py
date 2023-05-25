# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Wed 24 May 2023 15:50:24 BST

from pathlib import Path
import pickle
import re
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import nmrespy as ne
from utils import transfer

COLORS = mpl.rcParams["axes.prop_cycle"].by_key()["color"]
COLORS[0], COLORS[1], COLORS[2], COLORS[4] = COLORS[1], COLORS[4], COLORS[0], COLORS[2]


def remove_noise_coomponents(result):
    damp_thold = 0.7
    damps = result[:, 3]
    result = np.delete(result, np.where(damps < damp_thold), axis=0)
    amp_thold = 0.1
    amps = result[:, 0]
    result = np.delete(result, np.where(amps < amp_thold), axis=0)
    return result


def get_colors(idxs):
    return [COLORS[i] for i in idxs]



offset = 250.
result_dir = Path("~/Documents/DPhil/results/onedim/simulated").expanduser()
estimators, thetas = [], []
for i in range(1, 6):
    estimators.append(
        ne.Estimator1D.from_pickle(result_dir / f"estimators/estimator_{i}")
    )
    with open(result_dir / f"params/params_{i}.pkl", "rb") as fh:
        thetas.append(pickle.load(fh))

data_ylim = (10, 305)
osc_ylim = (-30, 305)

data_h = data_ylim[1] - data_ylim[0]
osc_h = osc_ylim[1] - osc_ylim[0]
ax_heights = [data_h, osc_h, osc_h, osc_h]
ax_height_sum = sum(ax_heights)
height_ratios = [h / ax_height_sum for h in ax_heights]

fig, axs = plt.subplots(
    nrows=4,
    ncols=5,
    gridspec_kw={
        "left": 0.005,
        "right": 0.995,
        "bottom": 0.05,
        "top": 0.99,
        "hspace": 0.0,
        "wspace": 0.04,
        "height_ratios": height_ratios
    },
    figsize=(9, 5),
)

kwargs = {
    "xaxis_unit": "hz",
    "spectrum_line_kwargs": {"linewidth": 0.6},
    "oscillator_line_kwargs": {"linewidth": 0.6},
    "residual_line_kwargs": {"linewidth": 0.6},
    "plot_model": False,
}

residual_shifts = 5 * [30.]  #[30., 26., 26., 34., 50.]
ytops = [245, 330, 280, 340, 310]
oscillator_cols_mpm = [
    get_colors([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1]),
    get_colors([0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0]),
    get_colors([2, 0, 0, 3, 3, 3, 0, 4, 4, 1, 1, 1, 1, 1, 0, 0, 0, 4, 4, 0]),
    get_colors([0, 1, 5, 5, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 4, 0, 0, 0, 0, 2, 0]),
    get_colors([0, 0, 0, 0, 0, 1, 1, 5, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0]),
]
oscillator_cols_nlp = [
    get_colors([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1]),
    get_colors([0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0]),
    get_colors([2, 0, 0, 3, 3, 3, 0, 4, 1, 1, 1, 1, 1, 0, 0, 0, 4, 0, 0, 0]),
    get_colors([0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 4, 4, 0, 0, 0, 0, 2, 0]),
    get_colors([0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0]),
]
oscillator_cols_true = [
    get_colors([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1]),
    get_colors([0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0]),
    get_colors([2, 2, 0, 0, 3, 3, 0, 4, 4, 1, 1, 1, 1, 1, 0, 0, 0, 4, 4, 0]),
    get_colors([0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 4, 4, 4, 0, 0, 0, 0, 2, 2, 0]),
    get_colors([0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0]),
]

for i, (
    ax_col,
    estimator,
    theta,
    residual_shift,
    ytop,
    osc_cols_mpm,
    osc_cols_nlp,
    osc_cols_true,
) in enumerate(zip(
    axs.T,
    estimators,
    thetas,
    residual_shifts,
    ytops,
    oscillator_cols_mpm,
    oscillator_cols_nlp,
    oscillator_cols_true,
)):
    kwargs["residual_shift"] = residual_shift
    nlp_params = estimator.get_params()
    mpm_params = estimator.get_results()[0].trajectory[0]
    remove_noise_coomponents(mpm_params)
    remove_noise_coomponents(nlp_params)

    # Data
    _ax = estimator.plot_result(oscillator_colors=osc_cols_mpm, **kwargs)[1][0, 0]
    for line in _ax.lines:
        if line.get_color() != "#000000":
            line.remove()
    transfer(_ax, axs[0, i], fig)

    # MPM result
    estimator._results[0].params = remove_noise_coomponents(mpm_params)
    _ax = estimator.plot_result(oscillator_colors=osc_cols_mpm, **kwargs)[1][0, 0]
    _ax.lines[0].remove()
    transfer(_ax, axs[1, i], fig)

    # NLP result
    estimator._results[0].params = remove_noise_coomponents(nlp_params)
    _ax = estimator.plot_result(oscillator_colors=osc_cols_nlp, **kwargs)[1][0, 0]
    _ax.lines[0].remove()
    transfer(_ax, axs[2, i], fig)

    # True parameters
    estimator._results[0].params = theta
    _ax = estimator.plot_result(oscillator_colors=osc_cols_true, **kwargs)[1][0, 0]
    _ax.lines[0].remove()
    transfer(_ax, axs[3, i], fig)

    for j, ax in enumerate(ax_col):
        if j == 0:
            ax.set_ylim(data_ylim)
        else:
            ax.set_ylim(osc_ylim)
        if j == 3:
            ax.set_xticklabels(
                [
                    str(int(
                        re.search(r"(\d+)", tick.get_text()).group(0)
                    ) - 250)
                    for tick in ax.get_xticklabels()
                ]
            )
        else:
            ax.set_xticklabels([])


for i, ax in enumerate(axs.T[0]):
    ax.text(0.02, 0.88, f"\\textbf{{{chr(97 + i)}.}}", transform=ax.transAxes)
for i, ax in enumerate(axs[0]):
    ax.text(
        0.98, 0.88, f"\\textbf{{Run {i + 1}}}", transform=ax.transAxes, ha="right",
        bbox={"facecolor": "w", "pad": 0., "edgecolor": "none"},
    )
fig.text(0.5, 0.005, "Hz", transform=fig.transFigure, ha="center", fontsize=8)

fig.savefig("figures/mpm_vs_nlp/mpm_vs_nlp.pdf")
