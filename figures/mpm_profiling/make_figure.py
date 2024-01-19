# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Fri 19 Jan 2024 16:23:07 EST

import os
import re
import pickle
import subprocess
import sys
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerTuple
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit

scatter_kwargs = dict(
    s=12, lw=0.5, edgecolor="k",
)

TOTAL_TIME_REGEX = re.compile(r"Total time: (\d+\.\d+) s")
MEMPROF_FILE_REGEX = re.compile(r".*_(\d+)pts\.memprof")
MEM_REGEX = re.compile(r"^MEM (\d+\.\d+)")

def linear(x, a, b):
    return a * x + b

def quadratic(x, a, b):
    return a * x ** 2 + b

def cubic(x, a, b):
    return a * x ** 3 + b

def get_var(directory, var):
    var_set = set()
    for file in directory.iterdir():
        match = re.search(r"(\d+)" + var, file.name)
        if match is not None:
            var_set.add(int(match.group(1)))
    var_list = np.array(sorted(list(var_set)))
    return var_list


colors = mpl.rcParams["axes.prop_cycle"].by_key()["color"]

pkl = False
pkldir = Path("figures/mpm_profiling/pickled").resolve()
if not pkldir.is_dir():
    os.mkdir(pkldir)
    pkl = True


fig, axs = plt.subplots(
    ncols=3,
    nrows=2,
    gridspec_kw={
        "left": 0.065,
        "top": 0.995,
        "bottom": 0.06,
        "right": 0.995,
        "wspace": 0.2,
        "hspace": 0.15,
    },
    figsize=(6, 5),
)
axs[1, 2].remove()

#### MPM Times
if pkl:
    directory = Path("~/Documents/DPhil/results/profiling/mpm_profiles").expanduser()
    mpm_pts = get_var(directory, "pts")
    mpm_times = []
    for pts in mpm_pts:
        mpm_time = 0.
        for repeat in range(5):
            file = directory / f"mpm_{pts}pts_{repeat}.timeprof"
            result = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "line_profiler",
                    file,
                ],
                stdout=subprocess.PIPE,
                encoding="utf-8",
            ).stdout
            mpm_time += float(re.search(TOTAL_TIME_REGEX, result).group(1)) / 5.
        mpm_times.append(mpm_time)
    mpm_times = np.array(mpm_times)
    with open(pkldir / f"mpm_pts", "wb") as fh:
        pickle.dump(mpm_pts, fh)
    with open(pkldir / f"mpm_times", "wb") as fh:
        pickle.dump(mpm_times, fh)

else:
    with open(pkldir / f"mpm_pts", "rb") as fh:
        mpm_pts = pickle.load(fh)
    with open(pkldir / f"mpm_times", "rb") as fh:
        mpm_times = pickle.load(fh)

cutoff = 6
for p, t, m in zip((mpm_pts[:cutoff], mpm_pts[cutoff:]), (mpm_times[:cutoff], mpm_times[cutoff:]), ("s", "o")):
    axs[0, 0].scatter(
        p,
        t ** (1 / 3),
        marker=m,
        facecolor=colors[0],
        zorder=100,
        **scatter_kwargs,
    )
axs[0, 0].set_ylabel("wall clock time (\\unit{\\second})")

(a, b), _ = curve_fit(linear, mpm_pts[cutoff:], mpm_times[cutoff:] ** (1 / 3))
x = np.array([mpm_pts[0], mpm_pts[-1]])
axs[0, 0].plot(x, a * x + b, color=colors[0])

ticks = np.linspace(0, 15 ** (1 / 3), 10)
ticklabels = [f"{x ** 3:.2f}" for x in ticks]
axs[0, 0].set_yticks(ticks)
axs[0, 0].set_yticklabels(ticklabels)

##### MPM max memory
if pkl:
    max_mem = []
    func_regex = re.compile(r"^FUNC nmrespy\.mpm\._mpm_1d (\d+\.\d+)")
    for pts in mpm_pts:
        file = directory / f"mpm_{pts}pts.memprof"
        with open(file, "r") as fh:
            lines = fh.readlines()
        max_ = 0.
        for line in lines:
            mem_match = re.search(MEM_REGEX, line)
            if mem_match is not None:
                mem = float(mem_match.group(1))
                if mem > max_:
                    max_ = mem
                continue
            func_match = re.search(func_regex, line)
            if func_match is not None:
                start_mem = float(func_match.group(1))

        max_mem.append(max_ - start_mem)
    max_mem = np.array(max_mem) / 1024
    with open(pkldir / "mpm_mem", "wb") as fh:
        pickle.dump(max_mem, fh)

else:
    with open(pkldir / "mpm_mem", "rb") as fh:
        max_mem = pickle.load(fh)

axs[1, 0].scatter(
    mpm_pts,
    max_mem ** (1 / 2),
    zorder=100,
    **scatter_kwargs,
)

(a, b), _ = curve_fit(linear, mpm_pts, max_mem ** (1 / 2))
x = np.array([mpm_pts[0], mpm_pts[-1]])
axs[1, 0].plot(x, a * x + b, color=colors[0], zorder=0)
max_gibi_pines = 0.09375
n_thold_pines = ((max_gibi_pines ** 0.5) - b) / a
lineopts = {"ls": ":", "color": "k"}
axs[1, 0].plot((0, n_thold_pines), 2 * [max_gibi_pines ** 0.5], **lineopts, zorder=-1)
axs[1, 0].plot((n_thold_pines, n_thold_pines), (0., max_gibi_pines ** 0.5), **lineopts, zorder=-1)
axs[1, 0].set_xlim(left=axs[0, 0].get_xlim()[0])
axs[1, 0].set_ylim(bottom=0)
axs[1, 0].text(n_thold_pines + 64, 0.04, int(n_thold_pines), fontsize=6)
axs[1, 0].text(400, 0.09375 ** 0.47, "\\qty{96}{\\mebi\\byte}", fontsize=6)

ticks = np.linspace(0, 2. ** (1 / 2), 10)
ticklabels = [f"{x ** 2:.2f}" for x in ticks]
axs[1, 0].set_yticks(ticks)
axs[1, 0].set_yticklabels(ticklabels)
axs[1, 0].set_ylabel("peak memory usage (\\unit{\\gibi\\byte})", labelpad=2)

### MMEMPM timings
if pkl:
    directory = Path("~/Documents/DPhil/results/profiling/mmempm_profiles").expanduser()
    mmempm_pts = get_var(directory, "pts")
    mmempm_oscs = get_var(directory, "oscs")
    mmempm_times = {osc: [] for osc in mmempm_oscs}
    for pts in mmempm_pts:
        for oscs in mmempm_oscs:
            mmempm_time = 0.
            for repeat in range(5):
                file = directory / f"mmempm_{pts}pts_{oscs}oscs_{repeat}.timeprof"
                result = subprocess.run(
                    [
                        sys.executable,
                        "-m",
                        "line_profiler",
                        file,
                    ],
                    stdout=subprocess.PIPE,
                    encoding="utf-8",
                ).stdout
                mmempm_time += float(re.search(TOTAL_TIME_REGEX, result).group(1)) / 5.
                if oscs == 80 and pts <= 160:
                    print(f"{repeat}, {pts}, {mmempm_time}")
            mmempm_times[oscs].append(mmempm_time)

    with open(pkldir / "mmempm_pts", "wb") as fh:
        pickle.dump(mmempm_pts, fh)
    with open(pkldir / "mmempm_times", "wb") as fh:
        pickle.dump(mmempm_times, fh)

else:
    with open(pkldir / "mmempm_pts", "rb") as fh:
        mmempm_pts = pickle.load(fh)
    with open(pkldir / "mmempm_times", "rb") as fh:
        mmempm_times = pickle.load(fh)
        marker="s",

starts = {10: 0, 20: 0, 40: 1, 80: 4}
scatters = []
for i, (oscs, times) in enumerate(mmempm_times.items()):
    scatter_pair = []
    pts = mmempm_pts[starts[oscs]:]
    times = np.array(times)[starts[oscs]:] ** 0.5
    scatter_pair.append(
        axs[0, 1].scatter(
            pts[:-10],
            times[:-10],
            color=colors[i],
            marker="s",
            **scatter_kwargs,
        )
    )
    scatter_pair.append(
        axs[0, 1].scatter(
            pts[-10:],
            times[-10:],
            color=colors[i],
            **scatter_kwargs,
        )
    )
    scatters.append(tuple(scatter_pair))
    (a, b), _ = curve_fit(linear, pts[-10:], times[-10:])
    x = np.array([pts[0], pts[-1]])
    axs[0, 1].plot(x, a * x + b, color=colors[i], zorder=-1)

ticks = np.linspace(0, 50 ** 0.5, 10)
ticklabels = [f"{x ** 2:.2f}" for x in ticks]
axs[0, 1].set_yticks(ticks)
axs[0, 1].set_yticklabels(ticklabels)

legend = axs[0, 1].legend(
    scatters,
    [f"$M = {x}$" for x in mmempm_times.keys()],
    loc=(0.16, 0.7),
    fontsize=7,
    # ncol=2,
    columnspacing=0.5,
    handletextpad=0.5,
    handler_map={tuple: HandlerTuple(ndivide=None)},
)
legend.get_frame().set_linewidth(0.8)

##### MMEMPM SVD time
if pkl:
    svd_times = []
    for pts in mmempm_pts:
        svd_time = 0.
        for repeat in range(5):
            file = directory / f"mmempm_{pts}pts_20oscs_{repeat}.timeprof"
            result = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "line_profiler",
                    file,
                ],
                stdout=subprocess.PIPE,
                encoding="utf-8",
            ).stdout
            svd_time += float(re.search("(\d+\.\d+).*UM", result).group(1)) / 5.
        svd_times.append(svd_time)

    with open(pkldir / "mmempm_svd_times", "wb") as fh:
        pickle.dump(svd_times, fh)

else:
    with open(pkldir / "mmempm_svd_times", "rb") as fh:
        svd_times = 1e-6 * np.array(pickle.load(fh))

pts1 = mmempm_pts[:6]
pts2 = mmempm_pts[6:]

axs[0, 2].scatter(
    pts1,
    svd_times[:6],
    color=colors[2],
    marker="s",
    **scatter_kwargs,
)
axs[0, 2].scatter(
    pts2,
    svd_times[6:],
    color=colors[2],
    **scatter_kwargs,
)

(a1, b1), _ = curve_fit(quadratic, pts1, svd_times[:6])
(a2, b2), _ = curve_fit(quadratic, pts2, svd_times[6:])
x = np.linspace(pts1[0], pts2[-1], 100)
line1 = a1 * x ** 2 + b1
line2 = a2 * x ** 2 + b2
axs[0, 2].plot(x, line2, color=colors[2], zorder=-1)
ylim = axs[0, 2].get_ylim()
axs[0, 2].plot(x, line1, color=colors[2], ls=":", zorder=-1)
axs[0, 2].set_ylim(ylim)

##### MMEMPM max memory (M=20)
if pkl:
    max_mem = []
    func_regex = re.compile(r"^FUNC nmrespy\.mpm\._mpm_2d (\d+\.\d+)")
    for pts in mmempm_pts:
        file = directory / f"mmempm_{pts}pts_20oscs.memprof"
        with open(file, "r") as fh:
            lines = fh.readlines()
        max_ = 0.
        for line in lines:
            mem_match = re.search(MEM_REGEX, line)
            if mem_match is not None:
                mem = float(mem_match.group(1))
                if mem > max_:
                    max_ = mem
                continue
            func_match = re.search(func_regex, line)
            if func_match is not None:
                start_mem = float(func_match.group(1))

        max_mem.append(max_ - start_mem)

    # To gibi bytes
    max_mem = np.array(max_mem) / 1024
    with open(pkldir / "mmempm_mem", "wb") as fh:
        pickle.dump(max_mem, fh)

else:
    with open(pkldir / "mmempm_mem", "rb") as fh:
        max_mem = pickle.load(fh)

axs[1, 1].scatter(
    mmempm_pts,
    np.array(max_mem) ** 0.5,
    **scatter_kwargs,
)

(a, b), _ = curve_fit(linear, mmempm_pts, max_mem ** (1 / 2))
x = np.array([mmempm_pts[0], mmempm_pts[-1]])
axs[1, 1].plot(x, a * x + b, color=colors[0], zorder=0)

ticks = np.linspace(0, 4 ** 0.5, 10)
ticklabels = [f"{x ** 2:.2f}" for x in ticks]
axs[1, 1].set_yticks(ticks)
axs[1, 1].set_yticklabels(ticklabels)

ticks = [x * 1024 for x in range(1, 9)]
ticklabels = [str(x // 1024) for x in ticks]
for ax in axs[:, 0]:
    ax.set_xlabel("$\\nicefrac{N}{1024}$")
    ax.set_xticks(ticks)
    ax.set_xticklabels(ticklabels)

ticks = [x * 64 for x in range(1, 9)]
ticklabels = [str(x // 64) for x in ticks]
for ax in [*axs[:, 1]] + [axs[0, 2]]:
    ax.set_xlabel("$\\nicefrac{N^{(2)}}{64}$")
    ax.set_xticks(ticks)
    ax.set_xticklabels(ticklabels)

sqrts = [3, 2, 1, 2, 2]
for ax, sqrt in zip(axs.flatten(), sqrts):
    if sqrt == 1:
        txt = "linear"
    else:
        txt = f"\\hspace*{{2pt}}$\\sqrt[\\leftroot{{3}}\\uproot{{2}}{sqrt}]{{\\hspace*{{3pt}}}}$"
    txt = f"({txt} scale)"
    ax.text(
        0.02, 0.92, txt, va="top", rotation=90, transform=ax.transAxes, fontsize=8,
    )

labels = ("a1", "b1", "b2", "a2", "b3")
for label, ax in zip(labels, axs.flatten()):
    ax.text(0.02, 0.94, f"\\textbf{{{label}.}}", transform=ax.transAxes)


fig.savefig("figures/mpm_profiling/mpm_profiling.pdf")
