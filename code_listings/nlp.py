from typing import Any, Iterable, Optional, Tuple, Union
import numpy as np
import scipy as sp


class FunctionFactory:
    """Object which computes and memoises the objective, gradient and
    hessian for a given set of parameters.
    """
    def __init__(self, theta: np.ndarray, fun: callable, *args) -> None:
        self.theta = theta
        self._new_theta = True
        self.fun = fun
        self.obj = None
        self.grad = None
        self.hess = None
        self.args = args

    def _compute_if_needed(self):
        """Compute the objective, gradient and Hessian if they haven't been yet"""
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


class FunctionFactoryGaussNewton1D(FunctionFactory):
    def __init__(self, theta: np.ndarray, *args) -> None:
        super().__init__(theta, obj_grad_gauss_newton_hess_1d, *args)


def obj_grad_gauss_newton_hess_1d(
    theta: np.ndarray, *args: Tuple[int, int, np.ndarray],
) -> Tuple[float, np.ndarray, np.ndarray]:
    Y, M, tp = args  # Unpack args
    N = Y.shape[0]
    X_per_osc = np.exp(
        np.outer(tp, (2j * np.pi * theta[2 * M : 3 * M] - theta[3 * M :]))
    ) * (theta[:M] * np.exp(1j * theta[M : 2 * M]))

    # Jacobian: all first partial derivatives, ∂x/∂θᵢ
    jac = np.zeros((N, 4 * M), dtype="complex")
    jac[:, :M] = X_per_osc / theta[:M]  # a
    jac[:, M : 2 * M] = 1j * X_per_osc  # φ
    jac[:, 2 * M : 3 * M] = np.einsum("ij,i->ij", X_per_osc, 2j * np.pi * tp)  # f
    jac[:, 3 * M :] = np.einsum("ij,i->ij", X_per_osc, -tp) # η

    X = np.einsum("ij->i", X_per_osc)
    Y_minus_X = Y - X

    obj = (Y_minus_X.conj().T @ Y_minus_X).real  # ℱ(θ)
    grad = -2 * (Y_minus_X.conj().T @ jac).real  # ∇ℱ(θ)
    hess = 2 * (jac.conj().T @ jac).real  # ∇²ℱ(θ)

    # Incorporate phase variance
    phi = theta[M : 2 * M]
    cos_phi = np.cos(phi)
    cos_sum = np.sum(cos_phi)
    sin_phi = np.sin(phi)
    sin_sum = np.sum(sin_phi)
    R = np.sqrt(cos_sum ** 2 + sin_sum ** 2)
    pv_obj = 1 - (R / M)  # Var(φ)
    pv_grad = (sin_phi * cos_sum - cos_phi * sin_sum) / (M * R)  # ∂Var(φ)/∂φᵢ
    x = (sin_phi * cos_sum) - (cos_phi * sin_sum)  # ∂²Var(φ)/∂φᵢ∂φⱼ
    term_1 = np.outer(x, x) / (R ** 2)
    phi_array = np.zeros((M, M))
    phi_array[:] = phi
    term_2 = -np.cos(phi_array.T - phi_array)
    pv_hess = term_1 + term_2
    pv_hess[np.diag_indices(M)] += cos_phi * cos_sum + sin_phi * sin_sum
    pv_hess /= M * R
    obj += pv_obj
    grad[M : 2 * M] += pv_grad
    hess[M : 2 * M, M : 2 * M] += pv_hess

    return obj, grad, hess


def trust_ncg(
    theta0: np.ndarray,
    function_factory: FunctionFactory,
    args: Iterable[Any] = (),
    eta: float = 0.15,
    initial_trust_radius: Optional[float] = None,
    max_trust_radius: Optional[float] = None,
    epsilon: float = 1.e-8,
    max_iterations: int = 100,
    check_neg_amps_every: int = 10,
) -> Tuple[np.ndarray, np.ndarray]:
    theta = theta0
    M = theta.shape[0] // 4
    factory = function_factory(theta, *args)
    if initial_trust_radius is None:
        initial_trust_radius = 0.1 * factory.gradient_norm
    if max_trust_radius is None:
        max_trust_radius = 16 * initial_trust_radius
    trust_radius = min(initial_trust_radius, max_trust_radius)

    k = 0  # Current iterate
    amp_slice = slice(0, M)  # Slice to extract amplitudes

    while True:
        # Required accuracy of the computed solution
        epsi = min(0.5, np.sqrt(factory.gradient_norm)) * factory.gradient_norm
        z = np.zeros_like(theta)
        r = factory.gradient
        d = -r
        while True:
            Bd = factory.hessian @ d
            dBd = d.T @ Bd
            if dBd <= 0:
                ta, tb = get_boundaries(z, d, trust_radius)
                pa = z + ta * d
                pb = z + tb * d
                p = min(factory.model(pa), factory.model(pb))
                hits_boundary = True
                break
            r_sq = r.T @ r
            alpha = r_sq / dBd
            z_next = z + alpha * d
            if sp.linalg.norm(z_next) >= trust_radius:
                _, tb = get_boundaries(z, d, trust_radius)
                p = z + tb * d
                hits_boundary = True
                break
            r_next = r + alpha * Bd
            r_next_sq = r_next.T @ r_next
            if np.sqrt(r_next_sq) < epsi:
                hits_boundary = False
                p = z_next
                break
            beta_next = r_next_sq / r_sq
            d_next = -r_next + beta_next * d
            z = z_next
            r = r_next
            d = d_next

        predicted_value = factory.model(p)
        theta_proposed = theta + p
        factory_proposed = function_factory(theta_proposed, *args)
        actual_reduction = factory.objective - factory_proposed.objective
        predicted_reduction = factory.objective - predicted_value
        if predicted_reduction <= 0:
            # No improvement could be found: terminate
            negative_amps = False
            break
        rho = actual_reduction / predicted_reduction
        if rho < 0.25:
            # Quadratic model perfomring poorly: reduce TR
            trust_radius *= 0.25
        elif rho > 0.75 and hits_boundary:
            # Quadratic model perfomring well: increase TR
            trust_radius = min(2 * trust_radius, max_trust_radius)
        if rho > eta:
            # Accept update
            theta = theta_proposed
            factory = factory_proposed

        k += 1
        if (k % check_neg_amps_every == 0):
            neg_amps = np.where(theta[amp_slice] <= 0)[0]
            print(neg_amps)
            if neg_amps.size > 0:
                # Negative amps found: this run in order to purge
                negative_amps = True
                break
        if factory.gradient_norm < epsilon:
            # Convergence
            negative_amps = False
            break
        if k == max_iterations:
            # Maximum allowed iterations reached
            negative_amps = False
            break

    errors = np.sqrt(
        factory.objective * np.abs(np.diag(np.linalg.inv(factory.hessian)))
    )
    return theta, errors, negative_amps


def get_boundaries(z, d, trust_radius):
    a = d.T @ d
    b = 2 * z.T @ d
    c = (z.T @ z) - (trust_radius ** 2)
    aux = b + np.copysign(
        np.sqrt(b * b - 4 * a * c),
        b,
    )
    return sorted([-aux / (2 * a), -(2 * c) / aux])


def nlp(
    Y: np.ndarray,
    sw: float,
    offset: float,
    theta0: np.ndarray,
    max_iterations: int = 100,
    epsilon: float = 1.0e-8,
    eta: float = 0.15,
    initial_trust_radius: float = None,
    max_trust_radius: float = None,
    check_neg_amps_every: int = 25,
):
    norm = np.linalg.norm(Y)
    Y /= norm
    N = Y.shape[0]
    M = theta0.shape[0]
    theta0_vec = theta0.flatten(order="F")  # Flatten parameter array: (m, p) -> (p * m,)
    theta0_vec[:M] /= norm  # Normalise amplitudes
    theta0_vec[2 * M : 3 * M] -= offset  # Center frequencies

    # Extra arguments needed to compute the objective, grad, and Hessian:
    # FID, model order, timepoints sampled
    opt_args = [Y, M, np.linspace(0, float(N - 1) / sw, N)]

    while True:
        theta_vec, errors_vec, negative_amps = trust_ncg(
            theta0=theta0_vec,
            function_factory=FunctionFactoryGaussNewton1D,
            args=opt_args,
            eta=eta,
            initial_trust_radius=initial_trust_radius,
            max_trust_radius=max_trust_radius,
            epsilon=epsilon,
            max_iterations=max_iterations,
            check_neg_amps_every=check_neg_amps_every,
        )

        if negative_amps:
            negative_idx = list(np.where(theta_vec[:M] <= 0.)[0])
            # Determine indices of parameters to remove
            slice_ = []
            for i in negative_idx:
                slice_.extend([i * M + j for j in negative_idx])

            theta_vec = np.delete(theta_vec, slice_)
            M -= len(negative_idx)  # New model order
            opt_args[1] = M

        else:
            break

    theta = theta_vec.reshape((M, 4), order="F")
    errors = errors_vec.reshape((M, 4), order="F")
    errors_vec /= N - 1
    theta[:, 2] += offset  # Correct frequencies
    theta[:, 0] *= norm  # Re-scale amplitudes
    errors[:, 0] *= norm
    theta[:, 1] = (theta[:, 1] + np.pi) % (2 * np.pi) - np.pi  # Wrap phases

    return theta, errors


if __name__ == "__main__":
    import nmrespy as ne
    sw = 200.
    offset = 0.
    pts = 1024
    expinfo = ne.ExpInfo(1, 200., 0., default_pts=pts)
    fid = expinfo.make_fid(
        np.array([
            [1., 0., 20., 0.2],
        ]),
        snr=30.,
    )
    theta0 = np.array([
        [1.1, 0.1, 20.3, 0.25],
    ])
    print(ne.nlp.nonlinear_programming(expinfo, fid.copy(), theta0))
    print(nlp(fid.copy(), sw, offset, theta0, max_iterations=1000, initial_trust_radius=1.))
