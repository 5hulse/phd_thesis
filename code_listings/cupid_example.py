from pathlib import Path
import nmrespy as ne

datapath = Path("~/data/camphor/").expanduser()
estimator = ne.Estimator2DJ.new_bruker(datapath / "1")  |\label{ln:newbruker}|
estimator.phase_data(p0=5.238, p1=-6.262)  |\label{ln:phase}|
estimator.baseline_correction()  |\label{ln:baseline}|

regions = [  |\label{ln:regionstart}|
    (2.55, 2.475), (2.35, 2.23), (2.09, 2.025),
    (1.95, 1.75), (1.7, 1.61), (1.375, 1.215),
]  |\label{ln:regionend}|
for region in regions:  |\label{ln:estimatestart}|
    estimator.estimate(
        region=region, noise_region=(5.25, 5.21), region_unit="ppm",
        nlp_trim=2048,
    )  |\label{ln:estimateend}|

estimator.to_pickle(datapath / "camphor_2dj_espy")  |\label{ln:pickle}|
estimator.write_cupid_to_bruker(datapath, expno=2)  |\label{ln:writecupid}|
thold = 2. * estimator.sw()[1] / estimator.default_pts[1]  # $\epsilon = \nicefrac{2\fswtwo}{\Ntwo}$
fig, _ = estimator.plot_result(multiplet_thold=thold)  |\label{ln:makefig}|
fig.savefig("camphor_2dj.pdf")  |\label{ln:savefig}|
