# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Tue 11 Jul 2023 23:52:02 BST

#!/usr/bin/python3

# Illustration of the steps taken in going from an FID
# to a frequency filtered signal

import copy
import os
import pathlib
import pickle
import subprocess
from typing import Dict

import colorsys
import numpy as np
from numpy.fft import fft, fftshift, ifft, ifftshift
from numpy.random import normal

import matplotlib as mpl
import matplotlib.colors as mc
from matplotlib.patches import ConnectionPatch, FancyArrowPatch, Rectangle
import matplotlib.pyplot as plt

import nmrespy as ne

COLORS = mpl.rcParams["axes.prop_cycle"].by_key()["color"]
mpl.rcParams["axes.ymargin"] = 0.070


def lighten_color(color, amount=0.5):
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = np.array(colorsys.rgb_to_hls(*mc.to_rgb(c)))
    return colorsys.hls_to_rgb(c[0],1-amount * (1-c[1]),c[2])

pts = 248
# CONSTRUCT SIGNALS
signals = {}
params = np.array(
    [
        [1, 0, 115, 7],
        [2, 0, 125, 7],
        [2, 0, 135, 7],
        [1, 0, 145, 7],
        [1.5, 0, -20, 7],
        [3, 0, -30, 7],
        [1.5, 0, -40, 7],
        [2, 0, -180, 7],
        [2, 0, -195, 7],
    ],
    dtype="float64",
)
expinfo = ne.ExpInfo(dim=1, default_pts=pts, sw=500., offset=0.)
fid = expinfo.make_fid(params, snr=22.)
filter_obj = ne.freqfilter.Filter(
    fid,
    expinfo,
    region=(165., 95.),
    noise_region=(50., 20.),
    sg_power=10.,
)
left, right = filter_obj._region[0]

signals["ve"] = ne.sig.ift(filter_obj._spectrum)
signals["spectrum"] = filter_obj._spectrum
signals["sg"] = filter_obj.sg
signals["sg_noise"] = filter_obj.sg_noise
signals["filtered_spectrum"], _ = filter_obj.get_filtered_spectrum(cut_ratio=None)
signals["filtered_ve"] = ne.sig.ift(signals["filtered_spectrum"])

region = filter_obj._region
noise_region = filter_obj._noise_region

x = np.arange(2 * pts)
lpad = x[20]
rpad = x[5]
fig, axs = plt.subplots(
    nrows=6,
    ncols=1,
    gridspec_kw={
        "hspace": 0.,
        "left": 0.05,
        "bottom": 0.01,
        "right": 0.99,
        "top": 0.99,
    },
    figsize=(5, 7),
)
axs = list(axs)
insetw = 0.25
axs.append(axs[0].inset_axes([0.518 - (insetw / 2), 0.6, insetw, 0.4]))
axs.append(axs[5].inset_axes([0.518 - (insetw / 2), 0.6, insetw, 0.4]))
axs = np.array(axs)

for ax in axs[:6]:
    ax.set_xlim(x[0] - lpad, x[-1] + rpad)
for ax in axs[6:]:
    ax.set_xlim(x[pts - 20], x[pts + 20])
    ax.set_ylim(-0.5, 0.5)

for ax in axs:
    ax.set_xticks([])
    ax.set_yticks([])

for ax in axs[[0, 6]]:
    ax.plot(x[:pts + 1], signals["ve"][:pts + 1].imag, color=lighten_color(COLORS[0]))
    ax.plot(x[:pts + 1], signals["ve"][:pts + 1].real, color=COLORS[0])
    ax.plot(x[pts:], signals["ve"][pts:].imag, color=lighten_color(COLORS[1]))
    ax.plot(x[pts:], signals["ve"][pts:].real, color=COLORS[1])

for ax in axs[[0, 5]]:
    rectangle = Rectangle(
        (x[pts - 20], -0.5),
        x[pts + 20] - x[pts - 20],
        1,
        facecolor="none",
        edgecolor="k",
        lw=0.8,
        zorder=100,
    )
    ax.add_patch(rectangle)

axs[0].plot(
    (x[pts - 20], 184.4),
    (0.5, 4.82),
    color="k",
    lw=0.8,
)
axs[0].plot(
    (x[pts + 20], 313.6),
    (0.5, 4.82),
    color="k",
    lw=0.8,
)

axs[5].plot(
    (x[pts - 20], 184.4),
    (0.5, 1.368),
    color="k",
    lw=0.8,
)
axs[5].plot(
    (x[pts + 20], 313.6),
    (0.5, 1.368),
    color="k",
    lw=0.8,
)


axs[1].plot(signals["spectrum"], color=COLORS[0])
axs[2].plot(signals["sg"], color=COLORS[0])
axs[3].plot(signals["sg_noise"], color=COLORS[0])
axs[4].plot(signals["filtered_spectrum"], color=COLORS[0])

center = (left + right) / 2
bottom = axs[2].get_ylim()[0]
for xpos in (left, right, center):
    axs[2].plot([xpos, xpos], [bottom, 1.1], ls=":", color="k")
axs[2].set_ylim(bottom=bottom)

axs[2].text(left - 2, 0.7, "$l_{\\text{idx}}$", ha="right", fontsize=8)
axs[2].text(right + 2, 0.7, "$r_{\\text{idx}}$", ha="left", fontsize=8)
axs[2].text(center + 2, 0.7, "$c$", ha="left", fontsize=8)
axs[2].text(center, 1.13, "$b$", ha="center", fontsize=8)
arrow = FancyArrowPatch(posA=(left, 1.1), posB=(right, 1.1), arrowstyle="<->", shrinkA=0, shrinkB=0, mutation_scale=10, lw=1)
axs[2].add_patch(arrow)
axs[2].set_ylim(top=1.29)

for ax in axs[[5, 7]]:
    ax.plot(x[:pts + 1], signals["filtered_ve"][:pts + 1].imag, color=lighten_color(COLORS[0]))
    ax.plot(x[:pts + 1], signals["filtered_ve"][:pts + 1].real, color=COLORS[0])
    ax.plot(x[pts:], signals["filtered_ve"][pts:].imag, color=lighten_color(COLORS[1]))
    ax.plot(x[pts:], signals["filtered_ve"][pts:].real, color=COLORS[1])

for ax in axs[6:]:
    ax.axvline(x[pts], color="k", ls=":")
for i, ax in enumerate(axs[:6]):
    ax.text(0.008, 0.88, f"\\textbf{{{chr(97 + i)}.}}", transform=ax.transAxes)

axs[4].set_ylim(axs[1].get_ylim())

for r, color in zip((region, noise_region), (lighten_color(COLORS[5], 0.2), "#e0e0e0")):
    r = r[0]
    axs[1].axvspan(r[0], r[1], zorder=0, color=color)

axs[0].set_yticks([-10, 0, 10])
axs[1].set_yticks([0, 150, 300])
axs[2].set_yticks([0, 1])
axs[3].set_yticks([-8, 0, 8])
axs[4].set_yticks([0, 150, 300])
axs[5].set_yticks([-5, 0, 5])

fig.savefig("figures/filtering/filtering.pdf")
