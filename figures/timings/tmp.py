# get_timings_hess1d.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Thu 20 Apr 2023 20:10:29 BST

import nmrespy as ne
import numpy as np

M = 64
pts = 8192
expinfo = ne.ExpInfo(
    dim=1,
    sw=500.,
    default_pts=pts,
)
params = np.zeros((64, 4), dtype="float64")
params[:, 0] = 1.
params[:, 2] = np.linspace(-200., 200, M)
params[:, 3] = 2
fid = expinfo.make_fid(params=params, snr=20.)
norm = np.linalg.norm(fid)
fid /= norm
params[:, 0] /= norm
active = params.reshape((4 * M,), order="F")
args = (
    fid,
    expinfo.get_timepoints(),
    M,
    np.array([]),
    [0, 1, 2, 3],
    False,
)
obj, grad, hess = ne.nlp._funcs.obj_grad_true_hess_1d(
    active,
    *args,
)
