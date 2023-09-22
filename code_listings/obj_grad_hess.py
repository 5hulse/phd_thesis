class FunctionFactory:
    """Object which computes and memoises the objective, gradient and
    hessian for a given set of parameters."""
    def __init__(self, theta: np.ndarray, fun: callable, *args) -> None:
        self.theta = theta
        self._new_theta = True
        self.fun = fun
        self.obj = None
        self.grad = None
        self.hess = None
        self.args = args

    def _compute_if_needed(self):
        """Compute the obj, grad and Hess if they haven't been yet"""
        if self.obj is None:
            self.obj, self.grad, self.hess = self.fun(self.theta, *self.args)

    def model(self, p) -> float:
        return self.objective + self.gradient @ p + 0.5 * (p.T @ self.hessian @ p)

    @property
    def objective(self) -> float:
        self._compute_if_needed()
        return self.obj

    @property
    def gradient(self) -> np.ndarray:
        self._compute_if_needed()
        return self.grad

    @property
    def gradient_norm(self) -> float:
        return sp.linalg.norm(self.gradient)

    @property
    def hessian(self) -> np.ndarray:
        self._compute_if_needed()
        return self.hess


class FunctionFactory1D(FunctionFactory):  !\label{ln:ff1}!
    def __init__(self, theta: np.ndarray, *args) -> None:
        super().__init__(theta, obj_grad_hess_1d, *args)  !\label{ln:ff2}!


def obj_grad_hess_1d(  # $\label{ln:obj-grad-hess-start}$
    theta: np.ndarray, *args: Tuple[int, int, np.ndarray],
) -> Tuple[float, np.ndarray, np.ndarray]:
    Y, tp = args  # Unpack args: FID, timepoints
    N = Y.shape[0]
    M = theta.shape[0] // 4
    X_per_osc = np.exp(
        np.outer(tp, (2j * np.pi * theta[2 * M : 3 * M] - theta[3 * M :]))
    ) * (theta[:M] * np.exp(1j * theta[M : 2 * M]))

    # === First partial derivatives ===
    d1 = np.zeros((N, 4 * M), dtype="complex")
    d1[:, :M] = X_per_osc / theta[:M]  # $\nicefrac{\partial x}{\partial a_m}$
    d1[:, M : 2 * M] = 1j * X_per_osc  # $\nicefrac{\partial x}{\partial \phi_m}$
    d1[:, 2 * M : 3 * M] = \  # $\nicefrac{\partial x}{\partial f_m}$
        np.einsum("ij,i->ij", X_per_osc, 2j * np.pi * tp)
    d1[:, 3 * M :] = np.einsum("ij,i->ij", X_per_osc, -tp)  # $\nicefrac{\partial x}{\partial \eta_m}$

    # === Second partial derivatives ===
    d2 = np.zeros((N, 10 * M), dtype="complex")
    # Note $\nicefrac{\partial^2 x}{\partial a_m^2} = 0$ always.
    d2[:, M : 2 * M] = 1j * d1[:, M : 2 * M]  # $\nicefrac{\partial^2 x}{\partial \phi_m^2}$
    d2[:, 2 * M : 3 * M] = \  # $\nicefrac{\partial^2 x}{\partial f_m^2}$
        np.einsum("ij,i->ij", d1[:, 2 * M : 3 * M], 2j * np.pi * tp)
    d2[:, 3 * M : 4 * M] = np.einsum("ij,i->ij", d2[: 3 * M :], -tp)  # $\nicefrac{\partial^2 x}{\partial \eta_m^2}$
    d2[:, 4 * M : 5 * M] = 1j * d1[:, :M]  # $\nicefrac{\partial^2 x}{\partial a_m \partial \phi_m}$
    d2[:, 5 * M : 6 * M] = 1j * d1[:, 2 * M : 3 * M]  # $\nicefrac{\partial^2 x}{\partial \phi_m \partial f_m}$
    d2[:, 6 * M : 7 * M] = \  # $\nicefrac{\partial^2 x}{\partial f_m \partial \eta_m}$
        np.einsum("ij,i->ij", d2[: 2 * M : 3 * M], -tp)
    d2[:, 7 * M : 8 * M] = d1[:, 2 * M : 3 * M] / theta[:M]  # $\nicefrac{\partial^2 x}{\partial a_m \partial f_m}$
    d2[:, 8 * M : 9 * M] = 1j * d1[:, 3 * M : 4 * M]  # $\nicefrac{\partial^2 x}{\partial \phi_m \partial \eta_m}$
    d2[:, 9 * M :] = d1[:, 3 * M : 4 * M] / theta[:M]  # $\nicefrac{\partial^2 x}{\partial a_m \partial \eta_m}$

    X = np.einsum("ij->i", X_per_osc)
    Y_minus_X = Y - X

    # === Objective ===
    obj = np.real(Y_minus_X.conj().T @ Y_minus_X)

    # === Gradient ===
    grad = -2 * np.real(Y_minus_X.conj().T @ d1)

    # === Hessian ===
    hess = np.zeros((4 * M, 4 * M))
    # Compute non-zero elements of term $\text{\circled{2}}$ of the Hessian,
    # and assign to upper right triangle.
    hess_term_2 = -2 * np.real(np.einsum("ji,j->i", d2.conj(), Y_minus_X))
    term_2_row_idxs = []
    term_2_col_idxs = []
    for i in range(4):
        rows, cols = _diag_idx(4 * M, i * M)
        term_2_row_idxs.extend(rows)
        term_2_col_idxs.extend(cols)
    term_2_idxs = (term_2_row_idxs, term_2_col_idxs)
    hess[term_2_idxs] = hess_term_2
    # Add the transpose to fill left botton triangle
    # Halve the main diagonal before transposing to avoid doubling
    main_diag = _diag_idx(4 * M, 0)
    hess[main_diag] *= 0.5
    hess += hess.T
    # Compute term $\text{\circled{1}}$ of the Hessian
    hess_term_1 = 2 * np.real(np.einsum("ki,kj->ij", d1.conj(), d1))
    hess += hess_term_1

    # === Phase variance and derivatives ===
    phi = theta[M : 2 * M]
    cos_phi = np.cos(phi)
    cos_sum = np.sum(cos_phi)
    sin_phi = np.sin(phi)
    sin_sum = np.sum(sin_phi)
    R = np.sqrt(cos_sum ** 2 + sin_sum ** 2)
    pv_obj = 1 - (R / M)
    pv_grad = (sin_phi * cos_sum - cos_phi * sin_sum) / (M * R)
    x = (sin_phi * cos_sum) - (cos_phi * sin_sum)
    term_1 = np.outer(x, x) / (R ** 2)
    phi_array = np.zeros((M, M))
    phi_array[:] = phi
    term_2 = -np.cos(phi_array.T - phi_array)
    pv_hess = term_1 + term_2
    pv_hess[np.diag_indices(M)] += cos_phi * cos_sum + sin_phi * sin_sum
    pv_hess /= M * R
    obj += pv_obj  # $\mathcal{F}_{\phi}\left(\symbf{\theta}\right)$
    grad[M : 2 * M] += pv_grad  # $\nabla\mathcal{F}_{\phi}\left(\symbf{\theta}\right)$
    hess[M : 2 * M, M : 2 * M] += pv_hess  # $\nabla^2\mathcal{F}_{\phi}\left(\symbf{\theta}\right)$

    return obj, grad, hess  # $\label{ln:obj-grad-hess-end}$


def _diag_idx(size: int, disp: int) -> Tuple[np.ndarray, np.ndarray]:
    """Return the indices of an array's kth diagonal.

    Parameters
    ----------
    size
        The size of the sqaure matrix.

    disp
        Displacement from the main diagonal.

    Returns
    -------
    rows
        0-axis coordinates of indices.

    cols
        1-axis coordinates of indices.
    """
    rows, cols = np.diag_indices(size)
    if disp < 0:
        return rows[-disp:], cols[:disp]
    elif disp > 0:
        return rows[:-disp], cols[disp:]
    else:
        return rows, cols
