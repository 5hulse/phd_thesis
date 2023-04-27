# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Wed 26 Apr 2023 16:02:01 BST

# Had to expand TeX's memory for this.
# See Leo Liu's comment at:
# https://tex.stackexchange.com/questions/7953/how-to-expand-texs-main-memory-size-pgfplots-memory-overload

import colorsys
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt, colors
import nmrespy as ne


def lighten_color(color, amount=0.5):
    try:
        c = colors.cnames[color]
    except:
        c = color
    c = np.array(colorsys.rgb_to_hls(*colors.to_rgb(c)))
    return colorsys.hls_to_rgb(c[0],1-amount * (1-c[1]),c[2])


mpl.rcParams["xtick.major.pad"] = 1
mpl.rcParams["ytick.major.pad"] = 1
mpl.rcParams["axes.ymargin"] = 0.070

COLORS = mpl.rcParams["axes.prop_cycle"].by_key()["color"]
red, blue = [lighten_color(c, 1.2) for c in COLORS[:2]]

full_cmap = colors.LinearSegmentedColormap.from_list("bwr_custom", [blue, "#ffffff", red])
half_cmap = colors.LinearSegmentedColormap.from_list("wr_custom", ["#ffffff", red], N=128)

pts = (24, 48)
fpts = (pts[0], pts[1] * 2 - 1)
expinfo = ne.ExpInfo(
    dim=2,
    sw=(50., 200.),
    nuclei=(None, "1H"),
    default_pts=pts,
)
fid = ne.sig.add_noise(
    expinfo.make_fid(
        np.array(
            [
                [2, 0, 10, 45, 3, 7],
                [4, 0, 0, 55, 3, 7],
                [2, 0, -10, 65, 3, 7],
                [2, 0, 10, -40, 3, 7],
                [4, 0, 0, -50, 3, 7],
                [2, 0, -10, -60, 3, 7],
            ]
        )
    ),
    8.,
)

estimator = ne.Estimator2DJ(fid, expinfo)
_, f2 = estimator.get_shifts(pts=fpts)
t1, _ = estimator.get_timepoints(pts=fpts)
t1_1d = t1[:, 0]
f2_1d = f2[0]
filter_ = ne.freqfilter.Filter(
    estimator.data,
    estimator.expinfo,
    region=(None, (70., 20.)),
    noise_region=(None, (90., 85.)),
    twodim_dtype="hyper",
    sg_power=10.,
)
spectrum = filter_.spectrum.real
spectrum_max = np.amax(spectrum)
sg = filter_.sg.real
noise = filter_.sg_noise.real
filtered_spectrum = filter_._filtered_unfixed_spectrum.real

fig = plt.figure(figsize=(6, 2.5))

l_pad = 0.05
r_pad = 0.01
t_pad = 0.075
b_pad = 0.08
ax_space = 0.05
cb_space = 0.02
h_ratios = [2.5, 1, 0.15]
ax_w = (1 - (l_pad + r_pad + 3 * ax_space)) / 4
total_h = 1 - (t_pad + b_pad + cb_space)
hr_sum = sum(h_ratios)
twodim_h, onedim_h, cb_h = [(x * total_h) / hr_sum for x in h_ratios]

lefts = [l_pad + i * (ax_w + ax_space) for i in range(4)]
twodim_bottom = b_pad + cb_h + cb_space + onedim_h
onedim_bottom = b_pad + cb_h + cb_space

axs = []
cb_axs = []
for i, left in enumerate(lefts):
    axs.append(
        [
            fig.add_axes(
                [left, twodim_bottom, ax_w, twodim_h]
            ),
            fig.add_axes(
                [left, onedim_bottom, ax_w, onedim_h]
            ),
        ]
    )
    if i < 2:
        cb_axs.append(
            axs[i][0].inset_axes(
                [left, b_pad, ax_w, cb_h],
                transform=fig.transFigure,
            )
        )

for ax_group in axs:
    for ax in ax_group:
        ax.set_xticks([])
        ax.set_yticks([])

ax_geom = axs[0][0].get_position()
x_loc = ((ax_geom.x0 + ax_geom.x1) / 2, ax_geom.y1)
y_loc = (ax_geom.x0, (ax_geom.y0 + ax_geom.y1) / 2)
fig.text(
    x_loc[0],
    x_loc[1] + 0.01,
    "$F^{(2)}$",
    transform=fig.transFigure,
    ha="center",
    va="bottom",
)
fig.text(
    y_loc[0] - 0.005,
    y_loc[1],
    "$t^{(1)}$",
    rotation="horizontal",
    transform=fig.transFigure,
    ha="right",
    va="center",
)

spec_norm = colors.CenteredNorm(halfrange=np.abs(np.amax(spectrum)))
for i, obj in enumerate((spectrum, sg, noise, filtered_spectrum)):
    norm = spec_norm if i != 1 else None
    cmap = full_cmap if i != 1 else half_cmap
    pcm = axs[i][0].pcolormesh(obj, norm=norm, cmap=cmap)
    axs[i][1].plot(obj[0], color="k")
    if i < 2:
        fig.colorbar(pcm, ax=axs[i][0], cax=cb_axs[i], orientation="horizontal")

ax_min = 1000.
ax_max = -1000.
for i in (0, 2, 3):
    mn, mx = axs[i][1].get_ylim()
    ax_min = mn if mn < ax_min else ax_min
    ax_max = mx if mx > ax_max else ax_max

high_tick = 150.
for i in (0, 2, 3):
    axs[i][1].set_ylim(ax_min, ax_max)
    axs[i][1].set_yticks([0, high_tick])
cb_axs[0].set_xticks([-high_tick, 0, high_tick])
cb_axs[1].set_xticks([0, 0.5, 1])

axs[1][1].set_yticks([0, 1])

for i in range(4):
    pos = axs[i][0].get_position()
    x = pos.x0 + 0.005
    y = pos.y1 - 0.05
    fig.text(x, y, f"\\textbf{{{chr(97 + i)}.}}", ha="left", va="bottom")
    if i != 1:
        fig.text(
            pos.x1 - 0.005,
            pos.y1 - 0.008,
            "*",
            ha="right",
            va="top",
        )
    else:
        fig.text(
            pos.x1 - 0.005,
            pos.y1 - 0.008,
            "\\textsuperscript{\\dag}",
            ha="right",
            va="top",
        )

cb0_pos = cb_axs[0].get_position()
cb1_pos = cb_axs[1].get_position()
fig.text(
    x=cb0_pos.x1 + 0.005,
    y=(cb0_pos.y0 + cb0_pos.y1) / 2,
    s="*",
    va="center",
)
fig.text(
    x=cb1_pos.x1 + 0.005,
    y=(cb1_pos.y0 + cb1_pos.y1) / 2,
    s="\\textsuperscript{\\dag}",
    va="center",
)

fig.savefig("figures/jres_filtering/jres_filtering.pdf")
