# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Fri 01 Sep 2023 13:01:15 BST

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
pkldir = Path("figures/nlp_profiling/pickled").resolve()
if not pkldir.is_dir():
    os.mkdir(pkldir)
    pkl = True


fig, axs = plt.subplots(
    ncols=2,
    nrows=2,
    gridspec_kw={
        "left": 0.1,
        "top": 0.99,
        "bottom": 0.06,
        "right": 0.99,
        "wspace": 0.2,
        "hspace": 0.15,
    },
    figsize=(4., 5),
)

### 1D timings
if pkl:
    directory = Path("~/Documents/DPhil/results/profiling/nlp1d_profiles").expanduser()
    nlp1d_pts = get_var(directory, "pts")
    nlp1d_oscs = get_var(directory, "oscs")
    nlp1d_times = {pts: [] for pts in nlp1d_pts}
    for oscs in nlp1d_oscs:
        for pts in nlp1d_pts:
            nlp1d_time = 0.
            for repeat in range(5):
                file = directory / f"nlp1d_{oscs}oscs_{pts}pts_{repeat}.timeprof"
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
                nlp1d_time += float(re.search(TOTAL_TIME_REGEX, result).group(1)) / 5.
            nlp1d_times[pts].append(nlp1d_time)

    with open(pkldir / "nlp1d_pts", "wb") as fh:
        pickle.dump(nlp1d_pts, fh)
    with open(pkldir / "nlp1d_oscs", "wb") as fh:
        pickle.dump(nlp1d_oscs, fh)
    with open(pkldir / "nlp1d_times", "wb") as fh:
        pickle.dump(nlp1d_times, fh)

else:
    with open(pkldir / "nlp1d_pts", "rb") as fh:
        nlp1d_pts = pickle.load(fh)
    with open(pkldir / "nlp1d_oscs", "rb") as fh:
        nlp1d_oscs = pickle.load(fh)
    with open(pkldir / "nlp1d_times", "rb") as fh:
        nlp1d_times = pickle.load(fh)
        marker="s",

scatters = []
for i, (oscs, times) in enumerate(nlp1d_times.items()):
    if i % 2 == 1:
        scatters.append(
            axs[0, 0].scatter(
                nlp1d_oscs,
                np.array(times) ** 0.5,
                **scatter_kwargs,
            )
        )

legend = axs[0, 0].legend(
    scatters,
    [f"$N = {x}$" for i, x in enumerate(nlp1d_times.keys()) if i % 2 == 1],
    loc=(0.15, 0.70),
    fontsize=7,
    columnspacing=0.,
    handletextpad=0.,
)
legend.get_frame().set_linewidth(0.8)

ticks = np.linspace(0, 2. ** 0.5, 10)
ticklabels = [f"{x ** 2:.2f}" for x in ticks]
axs[0, 0].set_yticks(ticks)
axs[0, 0].set_yticklabels(ticklabels)
axs[0, 0].set_xlabel("$M$")
axs[0, 0].set_xticks([i * 10 for i in range(1, 9)])

##### 1D Hessian max memory
if pkl:
    max_mem = {}
    func_regex = re.compile(r"^FUNC __main__\.hess_1d (\d+\.\d+)")
    for oscs in nlp1d_oscs:
        mem_for_osc = []
        for pts in nlp1d_pts:
            file = directory / f"nlp1d_{oscs}oscs_{pts}pts.memprof"
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

            mem_for_osc.append(max_ - start_mem)

        max_mem[oscs] = np.array(mem_for_osc) / 1024

    with open(pkldir / "hess1d_mem", "wb") as fh:
        pickle.dump(max_mem, fh)

else:
    with open(pkldir / "hess1d_mem", "rb") as fh:
        max_mem = pickle.load(fh)

scatters = []
for i, (oscs, mem) in enumerate(max_mem.items()):
    if i % 4 == 3:
        scatters.append(
            axs[1, 0].scatter(
                nlp1d_pts,
                mem,
                **scatter_kwargs,
            )
        )

legend = axs[1, 0].legend(
    scatters,
    [f"$M= {x}$" for i, x in enumerate(max_mem.keys()) if i % 4 == 3],
    loc=(0.15, 0.7),
    fontsize=7,
    columnspacing=0.,
    handletextpad=0.,
)
legend.get_frame().set_linewidth(0.8)


axs[1, 0].set_xticks([i * 1024 for i in range(1, 9)])
axs[1, 0].set_xticklabels([f"${i}$" for i in range(1, 9)])
axs[1, 0].set_xlabel("$\\nicefrac{N}{1024}$")

### 2D timings
if pkl:
    directory = Path("~/Documents/DPhil/results/profiling/nlp2d_profiles").expanduser()
    nlp2d_pts = get_var(directory, "pts")
    nlp2d_oscs = get_var(directory, "oscs")
    nlp2d_times = {pts: [] for pts in nlp2d_pts}
    for oscs in nlp2d_oscs:
        for pts in nlp2d_pts:
            nlp2d_time = 0.
            for repeat in range(5):
                file = directory / f"nlp2d_{oscs}oscs_{pts}pts_{repeat}.timeprof"
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
                nlp2d_time += float(re.search(TOTAL_TIME_REGEX, result).group(1)) / 5.
            nlp2d_times[pts].append(nlp2d_time)

    with open(pkldir / "nlp2d_pts", "wb") as fh:
        pickle.dump(nlp2d_pts, fh)
    with open(pkldir / "nlp2d_oscs", "wb") as fh:
        pickle.dump(nlp2d_oscs, fh)
    with open(pkldir / "nlp2d_times", "wb") as fh:
        pickle.dump(nlp2d_times, fh)

else:
    with open(pkldir / "nlp2d_pts", "rb") as fh:
        nlp2d_pts = pickle.load(fh)
    with open(pkldir / "nlp2d_oscs", "rb") as fh:
        nlp2d_oscs = pickle.load(fh)
    with open(pkldir / "nlp2d_times", "rb") as fh:
        nlp2d_times = pickle.load(fh)
        marker="s",

scatters = []
for i, (oscs, times) in enumerate(nlp2d_times.items()):
    if i % 4 == 3:
        scatters.append(
            axs[0, 1].scatter(
                nlp2d_oscs,
                np.array(times),
                **scatter_kwargs,
            )
        )

legend = axs[0, 1].legend(
    scatters,
    [f"$N^{{(2)}} = {x}$" for i, x in enumerate(nlp2d_times.keys()) if i % 4 == 3],
    loc=(0.15, 0.70),
    fontsize=7,
    columnspacing=0.,
    handletextpad=0.,
)
legend.get_frame().set_linewidth(0.8)

axs[0, 1].set_xticks([i * 10 for i in range(1, 9)])
axs[0, 1].set_xlabel("$M$")

##### 2D Hessian max memory
if pkl:
    max_mem = {}
    func_regex = re.compile(r"^FUNC __main__\.hess_2d (\d+\.\d+)")
    for oscs in nlp2d_oscs:
        mem_for_osc = []
        for pts in nlp2d_pts:
            file = directory / f"nlp2d_{oscs}oscs_{pts}pts.memprof"
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

            mem_for_osc.append(max_ - start_mem)

        max_mem[oscs] = np.array(mem_for_osc) / 1024

    with open(pkldir / "hess2d_mem", "wb") as fh:
        pickle.dump(max_mem, fh)

else:
    with open(pkldir / "hess2d_mem", "rb") as fh:
        max_mem = pickle.load(fh)

scatters = []
for i, (oscs, mem) in enumerate(max_mem.items()):
    if i % 4 == 3:
        scatters.append(
            axs[1, 1].scatter(
                nlp2d_pts[1::2],
                mem[1::2],
                **scatter_kwargs,
            )
        )

legend = axs[1, 1].legend(
    scatters,
    [f"$M= {x}$" for i, x in enumerate(max_mem.keys()) if i % 4 == 3],
    loc=(0.15, 0.7),
    fontsize=7,
    columnspacing=0.,
    handletextpad=0.,
)
legend.get_frame().set_linewidth(0.8)


axs[1, 1].set_xticks([i * 64 for i in range(1, 9)])
axs[1, 1].set_xticklabels([f"${i}$" for i in range(1, 9)])
axs[1, 1].set_xlabel("$\\nicefrac{N^{(2)}}{64}$")


sqrts = [2, 1, 1, 1]
for ax, sqrt in zip(axs.flatten(), sqrts):
    if sqrt == 1:
        txt = "linear"
    else:
        txt = f"\\hspace*{{2pt}}$\\sqrt[\\leftroot{{3}}\\uproot{{2}}{sqrt}]{{\\hspace*{{3pt}}}}$"
    txt = f"({txt} scale)"
    ax.text(
        0.02, 0.92, txt, va="top", rotation=90, transform=ax.transAxes, fontsize=8,
    )

axs[0, 0].set_ylabel("wall clock time (\\unit{\\second})", labelpad=2)
axs[1, 0].set_ylabel("peak memory usage (\\unit{\\gibi\\byte})", labelpad=2)

labels = ("a1", "b1", "a2", "b2")
for label, ax in zip(labels, axs.flatten()):
    ax.text(0.02, 0.94, f"\\textbf{{{label}.}}", transform=ax.transAxes)

fig.savefig("figures/nlp_profiling/nlp_profiling.pdf")
