# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Tue 02 May 2023 21:58:27 BST

from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import nmrespy as ne

contour_lw = 0.3
onedim_shift = 7.
cupid_shift = 23.
ps_yshift = 20.
ps_xshift = 3.
ytop = 52.
ax1_ylim = (15, -15)
plots_left = 0.055
plots_right = 0.99

colors = mpl.rcParams["axes.prop_cycle"].by_key()["color"]

estimator_dir = Path("~/Documents/DPhil/results/cupid/four_multiplets").expanduser()

fig, plot_axes = plt.subplots(
    nrows=2,
    ncols=5,
    gridspec_kw={
        "top": 0.99,
        "bottom": 0.09,
        "left": plots_left,
        "right": plots_right,
        "hspace": 0.,
        "wspace": 0.04,
        "height_ratios": (2, 1),
    },
    figsize=(6, 3),
)

estimator_dir = estimator_dir / "estimators"
for i, axs in enumerate(plot_axes.T):
    path_2dj = estimator_dir / f"estimator_{i}"
    estimator = ne.Estimator2DJ.from_pickle(path_2dj)
    _, axs_ = estimator.plot_result(
        multiplet_colors=colors,
        contour_base=0.7,
        contour_factor=1.6,
        contour_nlevels=10,
        multiplet_show_center_freq=False,
        multiplet_show_45=False,
    )

    path_ps = estimator_dir / f"estimator_pure_shift_{i}"
    estimator_ps = ne.Estimator2DJ.from_pickle(path_ps)
    ps_spectrum = estimator_ps.spectrum
    ps_shifts, = estimator_ps.get_shifts(unit="hz")

    axs_ = axs_[:, 0]
    k_count = 0
    for line in axs_[0].get_lines():
        ydata = line.get_ydata()
        ymin = np.amin(ydata)
        color = line.get_color()
        if color in ["k", "#000000"]:
            if k_count == 0:
                shift = onedim_shift
            elif k_count == 1:
                shift = cupid_shift
            k_count += 1
        else:
            shift = 0.

        axs[0].plot(
            ps_shifts + ps_xshift,
            ps_spectrum.real + cupid_shift,
            color="#c0c0c0",
        )
        axs[0].plot(
            line.get_xdata(),
            ydata - ymin + shift,
            color=line.get_color(),
        )

    for collection in axs_[1].collections:
        if collection.get_edgecolor().size == 0:
            pass
        else:
            for path in collection.get_paths():
                vertices = path.vertices
                x, y = vertices.T
                axs[1].plot(x, y, color="k", lw=contour_lw, zorder=0)

    for ax in axs:
        ax.set_xlim(axs_[0].get_xlim())

    # Have to hack scatter in.
    # Don't know how to extract from the axes
    multiplets = estimator.predict_multiplets()
    params = estimator.get_params()
    for (cf, mp), color in zip(reversed(multiplets.items()), colors):
        f1 = params[mp, 2]
        f2 = params[mp, 3]
        axs[1].scatter(f2, f1, s=4, edgecolor="none", facecolor=color, alpha=1)
        for ax in axs:
            ax.axvline(cf, color=color, ls="-", zorder=-1, lw=0.5)

    axs[0].text(0.985, 0.945, f"\\textbf{{Run {i + 1}}}", transform=axs[0].transAxes, fontsize=7, ha="right")

    axs[0].set_xticks([])
    axs[1].set_xticks([20, 0, -20])
    axs[0].set_yticks([])
    axs[0].spines["bottom"].set_visible(False)
    axs[1].spines["top"].set_visible(False)
    if i == 0:
        axs[1].set_ylabel("Hz")
    else:
        axs[1].set_yticks([])
    axs[0].set_ylim(top=ytop)
    axs[1].set_ylim(ax1_ylim)

plot_axes[1, 0].set_yticks([-10, -5, 0, 5, 10])

fig.text((plots_left + plots_right) / 2, 0.01, "\\textsuperscript{1}H (Hz)", ha="center", fontsize=8)  # noqa: E501

label_x = 0.06
label_ys = (0.96, 0.615, 0.44, 0.355)
for i, label_y in enumerate(label_ys, start=97):
    char = chr(i)
    fig.text(label_x, label_y, f"\\textbf{{{char}.}}")


fig.savefig("figures/four_multiplets/four_multiplets.pdf")
