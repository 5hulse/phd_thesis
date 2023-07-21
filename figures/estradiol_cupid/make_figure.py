# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Fri 21 Jul 2023 17:53:12 BST

from bruker_utils import BrukerDataset
from pathlib import Path
import matplotlib as mpl
from matplotlib.patches import ConnectionPatch
import numpy as np

import nmrespy as ne

from utils import (
    add_pure_shift_labels,
    get_pure_shift_labels,
    fix_linewidths,
    panel_labels,
    raise_axes,
)


# === CONFIGURATION ===
data_dir = Path("~/Documents/DPhil/data/estradiol").expanduser()
psyche_dir = data_dir / "6300/pdata/1"
estimator_pe_path = Path(
    "~/Documents/DPhil/results/cupid/estradiol_low_snr/estimator_postedit"
).expanduser()
estimator_path = Path(
    "~/Documents/DPhil/results/cupid/estradiol_low_snr/estimator"
).expanduser()


cupid_shift = 3.5e4
residual_shift = cupid_shift + 1.65e5
multiplet_shift = residual_shift + 1.7e4
onedim_shift = multiplet_shift + 7e3
ax_top = onedim_shift + 6e4

xticks = ([2.25], [2.05], [1.9, 1.8], [1.6], [1.4, 1.3, 1.2, 1.1])
thold_scale = 6.
label_x = 0.055
label_ys = (0.96, 0.7, 0.47, 0.395, 0.35)
psyche_shift = 1.e8
# =====================

estimator = ne.Estimator2DJ.from_pickle(estimator_pe_path)
estimator_tmp = ne.Estimator2DJ.from_pickle(estimator_path)
estimator_tmp.predict_multiplets(indices=[0], rm_spurious=True, thold=2., max_iterations=1, check_neg_amps_every=1)
estimator._results.pop(0)
estimator._results.insert(0, estimator_tmp._results[0])

colors = mpl.rcParams["axes.prop_cycle"].by_key()["color"]
colors.append("#808080")
mp_cols = [
    colors[i] for i in [
        0, 1, 2, 3, 4, 5, -1, 0, 1, 2, -1, 3, 4, -1, 0, 1, -1, 2, -1,
    ]
]
# estimator_tmp =
thold = thold_scale * estimator.default_multiplet_thold
fig, axs = estimator.plot_result(
    multiplet_thold=thold,
    region_unit="ppm",
    multiplet_colors=mp_cols,
    axes_bottom=0.065,
    axes_left=0.05,
    axes_right=0.995,
    axes_top=0.98,
    xaxis_label_height=0.01,
    axes_region_separation=0.035,
    marker_size=6.,
    contour_base=2.e8,
    contour_factor=1.2,
    contour_nlevels=10,
    multiplet_show_center_freq=True,
    multiplet_show_45=False,
    xaxis_ticks=[
        (0, (2.25, 2.2)),
        (1, (2.1, 2.05)),
        (2, (1.9, 1.85, 1.8, 1.75)),
        (3, (1.6, 1.55)),
        (4, (1.4, 1.35, 1.3, 1.25, 1.2, 1.15, 1.1, 1.05)),
    ],
    figsize=(6, 4.5),
)

fix_linewidths(axs, mpl.rcParams["lines.linewidth"])
raise_axes(axs, 1.98)

cupid_max = 0.
res = estimator.sw(unit="ppm")[1] / estimator.default_pts[1]
psyche_estimator = ne.Estimator1D.new_bruker("/home/simon/Documents/DPhil/data/estradiol/6300/pdata/1")
psyche_spectra = []
psyche_shifts = []
n_axs = len(axs[0])
for i, ax in enumerate(axs[0]):
    if i in [0]:
        cupid_idx = -2
        xs = [1]
    elif i == n_axs - 1:
        cupid_idx = -2
        xs = [0]
    else:
        cupid_idx = -3
        xs = [0, 1]

    ys = len(xs) * [0]

    cupid_line = ax.get_lines()[cupid_idx]
    # Number was found by printing relavent shift out in source code for
    # Estimator2DJ.plot_result
    line_max = np.amax(cupid_line.get_ydata()) - 56291.79
    if line_max > cupid_max:
        cupid_max = line_max
    cupid_line.set_ydata(cupid_line.get_ydata())
    # hack: data I generate seems to be off by resolution / 2
    cupid_line.set_xdata(cupid_line.get_xdata() + res / 2)

    onedim_idx = cupid_idx - 1
    onedim_line = ax.get_lines()[onedim_idx]
    # hack: data I generate seems to be off by resolution / 2
    onedim_line.set_xdata(onedim_line.get_xdata() + res / 2)

    multiplet_slice = slice(0, onedim_idx)
    shifts = onedim_line.get_xdata()
    slice_ = estimator.convert([None, [shifts[0], shifts[-1]]], "ppm->idx")[1]

    mp_lines = ax.get_lines()[multiplet_slice]
    if i == 0:
        mp_colors = iter([line.get_color() for line in mp_lines])
    for mp_line in mp_lines:
        mp_ydata = mp_line.get_ydata()
        # hack: data I generate seems to be off by resolution / 2
        mp_line.set_xdata(mp_line.get_xdata() + res / 2)

    left, right = ax.get_xlim()
    slc = psyche_estimator.convert([[left, right]], "ppm->idx")[0]
    psyche_shifts.append(psyche_estimator.get_shifts(unit="ppm")[0][slc[0] : slc[1] + 1])
    psyche_spectra.append(psyche_estimator.spectrum.real[slc[0] : slc[1] + 1])

psyche_max = max([np.amax(psyche_spectrum) for psyche_spectrum in psyche_spectra])
psyche_scale = cupid_max / psyche_max

for ax, psyche_shift, psyche_spectrum in zip(axs[0], psyche_shifts, psyche_spectra):
    ax.plot(psyche_shift, 2.e5 + psyche_scale * psyche_spectrum, color="k", lw=0.5)

panel_labels(fig, label_x, label_ys)

xs, ys, ss = get_pure_shift_labels(estimator, yshift=6.5e4, thold=3e4)
add_pure_shift_labels(axs[0], xs, ys, ss, fs=7)

fig.text(
    0.99,
    0.825,
    f"\\times {psyche_scale:.1f}",
    transform=fig.transFigure,
    fontsize=7,
    ha="right",
)

jres_dataset = BrukerDataset("/home/simon/Documents/DPhil/data/estradiol/32/pdata/1")
jres_spectrum = jres_dataset.data[::-1]
params = jres_dataset.get_parameters()
acqus, acqu2s = [params[file] for file in ("acqus", "acqu2s")]
jres_expinfo = ne.ExpInfo(
    dim=2,
    sw=(float(acqu2s["SW_h"]), float(acqus["SW_h"])),
    offset=(0., float(acqus["O1"])),
    sfo=(None, float(acqus["SFO1"])),
    nuclei=(None, "1H"),
    default_pts=jres_spectrum.shape,
)
f1, f2 = jres_expinfo.get_shifts(unit="ppm")
levels = [3.e3 * 1.8 ** i for i in range(10)]

for i, ax in enumerate(axs[1]):
    left, right = [jres_expinfo.convert((None, x), "ppm->idx")[1] for x in ax.get_xlim()]
    f1_ = f1[:, left:right]
    f2_ = f2[:, left:right]
    spec = jres_spectrum[:, left:right]
    ax.contour(f2_, f1_, spec, levels=levels, colors="k", linewidths=0.2)

for ax in axs[1]:
    ax.set_ylim(25, -25)

n_ax = axs[0].size
for i, (ax0, ax1) in enumerate(zip(axs[0], axs[1])):
    if i in [0, n_ax - 1]:
        lines = ax1.get_lines()[:-1]
    else:
        lines = ax1.get_lines()[:-2]
    for line in lines:
        line.remove()
        x = line.get_xdata()[0]
        con = ConnectionPatch(
            xyA=(x, ax1.get_ylim()[0]),
            xyB=(x, ax0.get_ylim()[0]),
            coordsA="data",
            coordsB="data",
            axesA=ax1,
            axesB=ax0,
            color=line.get_color(),
            lw=0.5,
            zorder=-1,
        )
        ax1.add_patch(con)
        ax0.axvline(x, color=line.get_color(), lw=0.5, zorder=-1)


fig.savefig("figures/estradiol_cupid/estradiol_cupid.pdf")
