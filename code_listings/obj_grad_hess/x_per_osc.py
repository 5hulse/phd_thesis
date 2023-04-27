exp_Z = np.exp(
    np.outer(
        tp,  # $\symbf{\tau}$
        2j * np.pi * theta[2 * M : 3 * M] - theta[3 * M:],  # $\symbf{z}$
    )
)
alpha = theta[: M] * np.exp(1j * theta[M : 2 * M]
X_m = exp_Z * alpha
