# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Thu 27 Apr 2023 14:04:36 BST

import re
import subprocess
import sys
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

colors = mpl.rcParams["axes.prop_cycle"].by_key()["color"]

fig, axs = plt.subplots(
    nrows=2,
    ncols=4,
)
axs = axs.flatten()

directory = Path("~/Documents/DPhil/results/profiling/mpm_profiles").resolve()

pts_list = []
regex = re.compile(r"mpm1d_(\d+)_pts\.memprof")
for file in directory.iterdir():
    match = re.search(regex, file.name)
    if match is not None:
        pts_list.append(int(match.group(1)))
pts_list = np.array(sorted(pts_list))


#### MPM Times
total_regex = re.compile(r"Total time: (\d+\.\d+) s")

mpm_times = []
for pts in pts_list:
    mpm_time = 0.
    for repeat in range(5):
        file = directory / f"mpm1d_{pts}_pts_{repeat}.timeprof"
        result = subprocess.run(
            [
                "python",
                "-m",
                "line_profiler",
                file,
            ],
            stdout=subprocess.PIPE,
            encoding="utf-8",
        ).stdout
        mpm_time += float(re.search(total_regex, result).group(1)) / 5.
    mpm_times.append(mpm_time)

mpm_times = np.array(mpm_times)
cutoff = 6
for p, t, m in zip((pts_list[:cutoff], pts_list[cutoff:]), (mpm_times[:cutoff], mpm_times[cutoff:]), ("s", "o")):
    axs[0].scatter(
        p,
        t ** (1 / 3),
        edgecolor="k",
        marker=m,
        s=12,
        lw=0.5,
        facecolor=colors[0],
        zorder=100,
    )
linear = stats.linregress(pts_list[cutoff:], mpm_times[cutoff:] ** (1 / 3))
c, m = linear.intercept, linear.slope
x = np.array([pts_list[0], pts_list[-1]])
axs[0].plot(x, m * x + c, color=colors[0])
x = np.array([pts_list[-1], 2 * pts_list[-1]])
axs[0].plot(x, m * x + c, color=colors[0], ls=":")

for ax in axs[:2]:
    ax.set_xlim(0, 2 * 8192 + 512)
    ax.set_xlabel("$N^{(1)}")
    ax.set_xticks(list(range(2048, 16384 + 2024, 2048)))

axs[0].set_ylabel("wall clock time (\\unit{\\second}) (\\hspace*{2pt}$\\sqrt[\\leftroot{3}\\uproot{2}3]{\\hspace*{3pt}}$ scale)")
ticks = np.linspace(0, 140 ** (1 / 3), 10)
ticklabels = [f"{x ** 3:.2f}" for x in ticks]
axs[0].set_yticks(ticks)
axs[0].set_yticklabels(ticklabels)

##### MPM max memory
max_mem = []
mem_regex = re.compile(r"^MEM (\d+\.\d+)")
func_regex = re.compile(r"^FUNC nmrespy\.mpm\._mpm_1d (\d+\.\d+)")
for pts in pts_list:
    file = directory / f"mpm1d_{pts}_pts.memprof"
    with open(file, "r") as fh:
        lines = fh.readlines()
    max_ = 0.
    for line in lines:
        mem_match = re.search(mem_regex, line)
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
axs[1].scatter(
    pts_list,
    max_mem ** (1 / 2),
    s=12,
    edgecolor="k",
    lw=0.5,
    color=colors[1],
    zorder=100,
)
linear = stats.linregress(pts_list, max_mem ** (1 / 2))
c, m = linear.intercept, linear.slope
x = np.array([pts_list[0], pts_list[-1]])
axs[1].plot(x, m * x + c, color=colors[1])
x = np.array([pts_list[-1], 2 * pts_list[-1]])
axs[1].plot(x, m * x + c, color=colors[1], ls=":")
ticks = np.linspace(0, 7.5 ** (1 / 2), 10)
ticklabels = [f"{x ** 2:.2f}" for x in ticks]
axs[1].set_yticks(ticks)
axs[1].set_yticklabels(ticklabels)
axs[1].set_ylabel("peak memory usage (\\unit{\\gibi\\byte}) (\\hspace*{2pt}$\\sqrt[\\leftroot{3}\\uproot{2}2]{\\hspace*{3pt}}$ scale)")

### NLP timings
directory = Path("~/Documents/DPhil/results/profiling/nlp_profiles").resolve()

pts_list = []
M_list = []
regex = re.compile(r"(\d+)_oscs_(\d+)_pts\.memprof")
for file in directory.iterdir():
    match = re.search(regex, file.name)
    if match is not None:
        M, pts = int(match.group(1)), int(match.group(2))
        if pts not in pts_list:
            pts_list.append(pts)
        if M not in M_list:
            M_list.append(M)

pts_list = np.array(sorted(pts_list))
M_list = np.array(sorted(M_list))
print(pts_list, M_list)

txt = ""
ogh_times = np.zeros((pts_list.size, M_list.size), dtype="float64")
for i, pts in enumerate(pts_list):
    for j, M in enumerate(M_list):
        ogh_time = 0.
        for repeat in range(5):
            file = directory / f"{M}_oscs_{pts}_pts_{repeat}.timeprof"
            result = subprocess.run(
                [
                    "python",
                    "-m",
                    "line_profiler",
                    file,
                ],
                stdout=subprocess.PIPE,
                encoding="utf-8",
            ).stdout
            txt +=f"{result}\n"
            ogh_time += float(re.search(total_regex, result).group(1)) / 5.
        ogh_times[i, j] = ogh_time

for times in ogh_times.T:
    axs[2].scatter(pts_list, times)
for times in ogh_times:
    axs[3].scatter(M_list, times)

fig.savefig("figures/timings/timings.pdf")
