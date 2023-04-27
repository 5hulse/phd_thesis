from typing import Tuple
import numpy as np

def obj_grad_hess(
    theta: np.ndarray, *args,
) -> Tuple[float, np.ndarray, np.ndarray]:
    """Compute the objective, gradient and hessian for 2D data.

    Parameters
    ----------
    theta
        1D array of parameters, of length 2(1+D)M

    args
        Contains elements in the following order:

        * data: Array of the original FID data.
        * tp: The time-points the signal was sampled at
        * phase_variance: If ``True``, include the oscillator phase variance to
          the cost function.

    Returns
    -------
    obj: float
        Value of the objective.

    grad: numpy.ndarray
        Gradient of the objective.

    hess: numpy.ndarray
        Gauss-Newton Hessian of the objective.
    """
    # Unpack args
    data, tp, phasevar = args
    # Determine the data dimensionality and number of oscillators
    D = data.ndim
    m = theta.size // (2 * (1 + D))
