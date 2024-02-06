def mpm(y: np.ndarray, sw: float, offset: float, M: int = 0) -> np.array:
    """Matrix Pencil Method for estimating a 1D FID.

    Parameters
    ----------
    y
        FID.
    sw
        Spectral width (Hz).
    offset
        Transmitter offset (Hz).
    M
        Number of oscillators. If 0, this is estimated using the MDL.
    """
    N = y.size  # $N$
    L = N // 3  # $L$: pencil parameter
    norm = np.linalg.norm(y)  # $\left\lVert\symbf{y}\right\rVert$
    y /= norm  # $\nicefrac{\symbf{Y}}{\left\lVert\symbf{y}\right\rVert}$

    # SVD of $\Hy$
    col = Y[: N - L]  # First column of $\Hy$
    row = Y[N - L - 1 :]  # Last row of $\Hy$
    Hy = sp.linalg.hankel(col, row)
    U, sigma, _ = np.linalg.svd(Hy)  # Determine $\symbf{\sigma}$ and $\symbf{U}$

    # Optional model order selection
    if M == 0:
        M = mdl(sigma, N, L)

    # Eigendecomposition of $\symbf{U}_{M1}^{+}\symbf{U}_{M2}^{\vphantom{+}}$ to get signal poles, $\symbf{z}$
    UM = U[:, :M]
    UM1 = Um[:-1, :]
    UM2 = Um[1:, :]
    UM1inv = np.linalg.pinv(UM1)
    UM1invUM2 = UM1inv @ UM2
    z, _ = np.linalg.eig(UM1invUM2)

    # Determine complex amplitudes, $\symbf{\alpha}$
    Z = np.power.outer(z, np.arange(N)).T  # Vandermonde pole matrix
    alpha = np.linalg.pinv(Z) @ y

    # Determine $\bda$, $\bdphi$, $\bdf$, and $\bdeta$, and store in $M \times 4$ array
    M = z.size
    theta = np.zeros((M, 4), dtype="float64")
    theta[:, 0] = np.abs(alpha) * norm
    theta[:, 1] = np.arctan2(np.imag(alpha), np.real(alpha))
    theta[:, 2] = (sw / (2 * np.pi)) * np.imag(np.log(z)) + offset
    theta[:, 3] = -sw * np.real(np.log(z))
    theta = theta[np.argsort(theta[:, 2])]  # order by $\bdf$

    # Remove oscillators with negative damping factors
    neg_damp_idx = np.nonzero(damp < 0.)[0]
    theta = np.delete(theta, neg_damp_idx, axis=0)

    return theta
