def mpm(Y: np.ndarray, sw: float, offset: float, M: int = 0) -> np.array:
    """Matrix Pencil Method for estimating a 1D FID.

    Parameters
    ----------
    Y
        FID.

    sw
        Sweep width (Hz).

    offset
        Transmitter offset (Hz).

    M
        Number of oscillators. If 0, this is estimated using the MDL.
    """
    N = Y  # $N^{(1)}$
    L = N // 3  # $L^{(1)}$: pencil parameter
    norm = np.linalg.norm(Y)  # $\left\lVert\symbf{Y}\right\rVert$
    Y /= norm  # $\nicefrac{\symbf{Y}}{\left\lVert\symbf{Y}\right\rVert}$
    col = Y[: N - L]  # First column of $\symbf{H}_{\symbf{Y}}$
    row = Y[N - L - 1 :]  # Last row of $\symbf{H}_{\symbf{Y}}$
    HY = sp.linalg.hankel(col, row)  # $\symbf{H}_{\symbf{Y}}$
    _, S, Vh = np.linalg.svd(HY)  # $\symbf{\sigma}$ and $\symbf{V}^{\dagger}$
    V = Vh.T  # $\symbf{V}$

    if M == 0:
        M = mdl(sigma, N, L)

    VM = V[:, :M]  # $\symbf{V}_M$
    VM1 = Vm[:-1, :]  # $\symbf{V}_{M1}$
    VM2 = Vm[1:, :]  # $\symbf{V}_{M2}$
    VM1inv = np.linalg.pinv(VM1)  # $\symbf{V}_{M1}^{+}$
    VM1invVM2 = VM1inv @ VM2  # $\symbf{V}_{M1}^{+}\symbf{V}_{M2}^{\vphantom{+}}$
    z, _ = np.linalg.eig(VM1invVM2)  # $\symbf{z}^{(1)}$: signal poles
    Z = np.power.outer(z, np.arange(N)).T  # $\symbf{Z}^{(1)}$
    alpha = np.linalg.pinv(Z) @ Y  # $\symbf{\alpha}$: complex amplitudes

    # Extract amplitude, phase, frequency and damping factor
    amp = np.abs(alpha) * norm  # $\symbf{a}$
    phase = np.arctan2(np.imag(alpha), np.real(alpha))  # $\symbf{\phi}$
    freq = (sw / (2 * np.pi)) * np.imag(np.log(z)) + offset  # $\symbf{f}^{(1)}$
    damp = -sw * np.real(np.log(z))  # $\symbf{\eta}^{(1)}$
    theta = np.vstack((amp, phase, freq, damp)).T  # $\symbf{\theta}$, as a $M \times 4$ array

    # Remove negative damping factors
    neg_damp_idx = np.nonzero(damp < 0.)[0]
    theta = np.delete(theta, neg_damp_idx, axis=0)

    theta = theta[np.argsort(theta[:, 2])]  # order by $\symbf{f}^{(1)}$
    return theta
