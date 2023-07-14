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


class FunctionFactoryGaussNewton1D(FunctionFactory):  !\label{ln:ff1}!
    def __init__(self, theta: np.ndarray, *args) -> None:
        super().__init__(theta, obj_grad_gauss_newton_hess_1d, *args)  !\label{ln:ff2}!


def obj_grad_gauss_newton_hess_1d(
    theta: np.ndarray, *args: Tuple[int, int, np.ndarray],
) -> Tuple[float, np.ndarray, np.ndarray]:
    Y, tp = args  # Unpack args: FID, timepoints
    N = Y.shape[0]
    M = theta.shape[0] // 4
    X_per_osc = np.exp(
        np.outer(tp, (2j * np.pi * theta[2 * M : 3 * M] - theta[3 * M :]))
    ) * (theta[:M] * np.exp(1j * theta[M : 2 * M]))

    # Jacobian: all first partial derivatives, $\nicefrac{\partial\symbf{X}}{\partial\symbf{\theta}}$
    jac = np.zeros((N, 4 * M), dtype="complex")
    jac[:, :M] = X_per_osc / theta[:M]  # $\nicefrac{\partial\symbf{X}}{\partial\symbf{a}}$
    jac[:, M : 2 * M] = 1j * X_per_osc  # $\nicefrac{\partial\symbf{X}}{\partial\symbf{\phi}}$
    jac[:, 2 * M : 3 * M] = \
        np.einsum("ij,i->ij", X_per_osc, 2j * np.pi * tp)  # $\nicefrac{\partial\symbf{X}}{\partial\symbf{f}^{(1)}}$
    jac[:, 3 * M :] = np.einsum("ij,i->ij", X_per_osc, -tp)  # $\nicefrac{\partial\symbf{X}}{\partial\symbf{\eta}^{(1)}}$

    X = np.einsum("ij->i", X_per_osc)
    Y_minus_X = Y - X

    obj = (Y_minus_X.conj().T @ Y_minus_X).real  # $\mathcal{F}\left(\symbf{\theta}\right)$
    grad = -2 * (Y_minus_X.conj().T @ jac).real  # $\nabla\mathcal{F}\left(\symbf{\theta}\right)$
    hess = 2 * (jac.conj().T @ jac).real  # $\nabla^2\mathcal{F}\left(\symbf{\theta}\right)$

    # === Determine phase variance and derivatives ===
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

    return obj, grad, hess
