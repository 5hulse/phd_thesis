# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Wed 03 Jan 2024 18:18:54 GMT

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
FORMAT = "portrait"

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

# New format: 3 runs on top, 2 runs on bottom, for portrait view
if FORMAT == "portrait":
    figsize = (6, 6.7)
    left = 0.01
    right = 0.99
    top = 0.995
    bottom = 0.04
    wspace = 0.01
    hspace = 0.03

    width = (1. - left - (1. - right) - 2 * wspace) / 3.
    height = (1. - bottom - (1. - top) - hspace) / 2.
    run_ax_heights = [height * ratio for ratio in height_ratios]

    top_left = left
    top_right = right
    top_bottom = bottom + height + hspace
    top_top = top

    bottom_extra = (width + wspace) / 2.
    bottom_left = left + bottom_extra
    bottom_right = right - bottom_extra
    bottom_bottom = bottom
    bottom_top = bottom + height

    rectangles = np.zeros((4, 5, 4))

    # Create top axes rectangles
    for i in range(3):
        l = top_left + i * (width + wspace)
        for j in range(4):
            b = top_bottom + sum(run_ax_heights[j + 1:])
            rectangles[j, i] = np.array([l, b, width, run_ax_heights[j]])

    # Create bottom axes rectangles
    for i in range(2):
        l = bottom_left + i * (width + wspace)
        for j in range(4):
            b = bottom_bottom + sum(run_ax_heights[j + 1:])
            rectangles[j, 3 + i] = np.array([l, b, width, run_ax_heights[j]])

    axs = np.empty((4, 5), dtype=mpl.axes.Axes)
    fig = plt.figure(figsize=figsize)
    for j in range(3, -1, -1):
        for i in range(5):
            axs[j, i] = fig.add_axes(rectangles[j, i])


# Old format: all 5 runs alongside each other, for landscape view
elif FORMAT == "landscape":
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

else:
    print("FORMAT must be \"landscape\" or \"portrait\"")
    exit()

kwargs = {
    "xaxis_unit": "hz",
    "spectrum_line_kwargs": {"linewidth": 0.6},
    "oscillator_line_kwargs": {"linewidth": 0.6},
    "residual_line_kwargs": {"linewidth": 0.6},
    "plot_model": False,
    "indices": [0],
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
# oscillator_cols_true = [
#     get_colors([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1]),
#     get_colors([0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0]),
#     get_colors([2, 2, 0, 0, 3, 3, 0, 4, 4, 1, 1, 1, 1, 1, 0, 0, 0, 4, 4, 0]),
#     get_colors([0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 4, 4, 4, 0, 0, 0, 0, 2, 2, 0]),
#     get_colors([0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0]),
# ]

for i, (
    ax_col,
    estimator,
    theta,
    residual_shift,
    ytop,
    osc_cols_mpm,
    osc_cols_nlp,
    # osc_cols_true,
) in enumerate(zip(
    axs.T,
    estimators,
    thetas,
    residual_shifts,
    ytops,
    oscillator_cols_mpm,
    oscillator_cols_nlp,
    # oscillator_cols_true,
)):
    kwargs["residual_shift"] = residual_shift
    nlp_params = estimator.get_params(indices=[0])
    mpm_params = estimator.get_results()[0].trajectory[0]
    remove_noise_coomponents(mpm_params)
    remove_noise_coomponents(nlp_params)

    # Data
    _ax = estimator.plot_result(oscillator_colors=osc_cols_mpm, **kwargs)[1][0, 0]
    for line in _ax.lines:
        if line.get_color() != "#000000":
            line.remove()
    transfer(_ax, axs[0, i], fig)

    # True parameters
    estimator._results[0].params = theta
    _ax = estimator.plot_result(oscillator_colors="k", **kwargs)[1][0, 0]
    _ax.lines[0].remove()
    transfer(_ax, axs[1, i], fig)

    # MPM result
    estimator._results[0].params = remove_noise_coomponents(mpm_params)
    _ax = estimator.plot_result(oscillator_colors=osc_cols_mpm, **kwargs)[1][0, 0]
    _ax.lines[0].remove()
    transfer(_ax, axs[2, i], fig)

    # NLP result
    estimator._results[0].params = remove_noise_coomponents(nlp_params)
    _ax = estimator.plot_result(oscillator_colors=osc_cols_nlp, **kwargs)[1][0, 0]
    _ax.lines[0].remove()
    transfer(_ax, axs[3, i], fig)

    for j, ax in enumerate(ax_col):
        if j in [0, 1]:
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


y = 0.88 if FORMAT == "landscape" else 0.84
for i, ax in enumerate(axs.T[0]):
    ax.text(0.02, y, f"\\textbf{{{chr(97 + i)}.}}", transform=ax.transAxes)
for i, ax in enumerate(axs[0]):
    ax.text(
        0.98, y, f"\\textbf{{Run {i + 1}}}", transform=ax.transAxes, ha="right",
        bbox={"facecolor": "w", "pad": 0., "edgecolor": "none"},
    )
fig.text(0.5, 0.005, "Hz", transform=fig.transFigure, ha="center", fontsize=8)

axs[0, 2].text(0.715, 0.27, "*", transform=axs[0, 2].transAxes)

for ax in axs[1]:
    for line in ax.lines:
        if line.get_color() == "#808080":
            line.remove()
fig.savefig(f"figures/mpm_vs_nlp/mpm_vs_nlp_{FORMAT}.pdf")
