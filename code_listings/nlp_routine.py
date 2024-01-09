# nlp_routine.py
# Simon Hulse
# simonhulse@protonmail.com
# Last Edited: Tue 09 Jan 2024 17:39:49 EST

def nlp(
    y: np.ndarray,
    sw: float,
    offset: float,
    theta0: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Nonlinear programming routine for 1D FID estimation.

    Parameters
    ----------
    y
        FID.
    sw
        Spectral width (Hz).
    offset
        Tranmitter offset (Hz).
    theta0
        Initial guess, of shape (M, 4)

    Returns
    -------
    theta
        Parameter estimate.
    errors
        Errors associated with theta.
    """
    norm = np.linalg.norm(y)
    y /= norm
    N = y.shape[0]
    M = theta0.shape[0]
    # Flatten parameter array: Fortran (column-wise) ordering
    theta0_vec = theta0.flatten(order="F")
    theta0_vec[:M] /= norm  # Normalise amplitudes
    # Remove transmitter offset from frequencies
    theta0_vec[2 * M : 3 * M] -= offset

    # Extra arguments needed to compute the objective, grad, and Hessian:
    # FID and timepoints sampled
    opt_args = [y, np.linspace(0, float(N - 1) / sw, N)]

    while True:
        theta_vec, errors_vec, negative_amps = trust_ncg(
            theta0=theta0_vec,
            function_factory=FunctionFactoryGaussNewton1D,
            args=opt_args,
        )

        if negative_amps:
            # Negative amps exist: remove these
            negative_idx = list(np.where(theta_vec[:M] <= 0.)[0])
            slice_ = []
            for idx in range(negative_idx):
                slice_.extend([i * M + idx for i in range(4)])
            theta_vec = np.delete(theta_vec, slice_)
            M -= len(negative_idx)  # New model order

        else:
            # No negative amps: routine complete
            break

    # Reshape parameter array back to (M, 4)
    theta = theta_vec.reshape((M, 4), order="F")
    errors = errors_vec.reshape((M, 4), order="F")
    errors_vec /= N - 1
    theta[:, 2] += offset  # Re-add transmitter offset to frequencies
    theta[:, 0] *= norm  # Re-scale amplitudes
    errors[:, 0] *= norm
    theta[:, 1] = (theta[:, 1] + np.pi) % (2 * np.pi) - np.pi  # Wrap phases

    return theta, errors
