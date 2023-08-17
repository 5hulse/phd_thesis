# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Tue 15 Aug 2023 23:19:45 BST

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
    nrows=1,
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

print(nlp1d_times)
for i, times in enumerate(nlp1d_times.values()):
    if i % 2 == 1:
        axs[0].scatter(
            nlp1d_oscs,
            np.array(times) ** 0.5,
            # color=colors[i],
            **scatter_kwargs,
        )
fig.savefig("figures/nlp_profiling/nlp_profiling.pdf")
