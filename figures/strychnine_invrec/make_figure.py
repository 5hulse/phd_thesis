# make_figure.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Thu 04 May 2023 19:51:54 BST

from pathlib import Path
import pickle
import numpy as np
import nmrespy as ne

RESULT_DIR = Path("~/Documents/DPhil/results/invrec/strychinine").expanduser()

estimator = ne.EstimatorInvRec.from_pickle(RESULT_DIR/ "estimator_new")
sfo, = estimator.sfo

with open(RESULT_DIR / "shifts.pkl", "rb") as fh:
    shifts = pickle.load(fh)
with open(RESULT_DIR / "couplings.pkl", "rb") as fh:
    couplings = pickle.load(fh)
with open(RESULT_DIR / "t1s.pkl", "rb") as fh:
    t1s = pickle.load(fh)

nspins = len(shifts)
coupling_array = np.zeros((nspins, nspins))
for (s1, s2, coupling) in couplings:
    coupling_array[s1 - 1, s2 - 1] = np.abs(coupling)
    coupling_array[s2 - 1, s1 - 1] = np.abs(coupling)

spans = [np.abs(np.sum(row)) / (2 * sfo) for row in coupling_array]

n_regions = len(estimator._results)
sep = 0.005
indices = list(range(n_regions))

fig, axs = estimator.plot_oscs_vs_fits(
    y_range=(0.3, 1.3),
    y_pts=128,
    indices=indices,
    contour_base=2.e-4,
    contour_factor=1.6,
    scale=6.,
    xaxis_unit="ppm",
    region_separation=sep,
    contour_kwargs={"linewidths": 0.2},
    oscillator_line_kwargs={"linewidth": 0.4},
    figsize=(6, 4),
)
merge_regions = estimator._plot_regions(indices, "ppm")[1]
merge_region_spans = estimator._get_3d_xaxis_spans(merge_regions, sep)

max_ = max([r[0] for r in merge_regions])
min_ = min([r[1] for r in merge_regions])

for shift, span, t1 in zip(shifts, spans, t1s):
    if max_ > shift > min_:
        xs_ppm = [shift - span, shift + span]
        xs_ax = estimator._transform_freq_to_xaxis(xs_ppm, merge_regions, merge_region_spans)
        print(xs_ax)
        axs[1].plot(xs_ax, [t1, t1], color="k")

fig.savefig("figures/strychnine_invrec/strychnine_invrec.pdf")
