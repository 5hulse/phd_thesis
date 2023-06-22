def mdl(sigma: np.ndarray, N: int, L: int) -> int:
    """Compute the Minimum Description Length

    Parameters
    ----------
    sigma
        Vector of singular values associated with the data matrix.

    N
        Number of points the signal comprises.

    L
        Pencil parameter.
    """
    mdl_vec = np.zeros(L)
    for k in range(L):
        mdl_vec[k] = (
            -N * np.sum(np.log(sigma[k:])) +
            N * (L - k) * np.log(np.sum(sigma[k:]) / (L - k)) +
            k * np.log(N) * (2 * L - k) / 2
        )
    M = sp.signal.argrelextrema(mdl_vec, np.less)[0][0]
    return M
