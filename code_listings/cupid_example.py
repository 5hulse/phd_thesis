from pathlib import Path
import nmrespy as ne

estimate_kwargs = {
    "noise_region": (5.25, 5.21),
    "region_unit": "ppm",
    "max_iterations": 100,
    "nlp_trim": 2048,
    "check_neg_amps_every": 25,
}
estimator = ne.Estimator2DJ.new_bruker(Path("~/Documents/data/camphor/1").expanduser())
estimator.phase_data(p0=5.238, p1=-6.262)
estimator.baseline_correction()
regions = [
    (2.55, 2.475),
    (2.35, 2.23),
    (2.09, 2.025),
    (1.95, 1.75),
    (1.7, 1.61),
    (1.375, 1.215),
]
for region in regions:
    estimator.estimate(region=region, **estimate_kwargs)
    estimator.to_pickle(estimator_path, force_overwrite=True)

thold = 2. * (estimator.sw()[1] / estimator.default_pts[1])
fig, _ = estimator.plot_result(multiplet_thold=thold)
