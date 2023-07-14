def trust_steihaug_toint(
    theta0: np.ndarray,
    function_factory: FunctionFactory,
    args: Iterable[Any] = (),
) -> Tuple[np.ndarray, np.ndarray, bool]:
    """Trust Region algorithm with Steihaug-Toint subroutine.

    Parameters
    ----------
    theta0
        Initial guess, with shape (M, 4).

    function_factory
        Object for computing the objective, gradient and Hessian.

    args
        Extra arguments required for computing the objective and its
        derivatives.

    Returns
    -------
    theta
        Parameter vector at termination.

    errors
        Errors associated with parameter vector.

    negative_amps
        Flag indicating whether or not termination occurred because
        negative amplitudes were detected.
    """
    theta = theta0
    M = theta.shape[0] // 4
    factory = function_factory(theta, *args)

    # === Define relevent parameters ===
    # These have been hard-coded, though in NMR-EsPy they are all
    # configurable.
    eta = 0.15,
    initial_trust_radius = 0.1 * factory.gradient_norm
    max_trust_radius = 16 * initial_trust_radius
    epsilon: float = 1.e-8,
    max_iterations: int = 200,
    check_neg_amps_every: int = 25,

    k = 0
    while True:
        # === Steihaug-Toint ===
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

        # === Assess effectiveness of update ===
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
            # Accept update: new iteration
            theta = theta_proposed
            factory = factory_proposed
            k += 1

        # === Check for termination criteria ===
        if (k % check_neg_amps_every == 0):  !\label{ln:negamp1}!
            neg_amps = np.where(theta[amp_slice] <= 0)[0]
            print(neg_amps)
            if neg_amps.size > 0:
                # Negative amps found: this run in order to purge
                negative_amps = True
                break  !\label{ln:negamp2}!

        if factory.gradient_norm < epsilon:
            # Convergence
            negative_amps = False
            break

        if k == max_iterations:
            # Maximum allowed iterations reached
            negative_amps = False
            break

    # Routine terminated: compute errors and return parameter array
    errors = np.sqrt(
        factory.objective *
        np.abs(np.diag(np.linalg.inv(factory.hessian)))
    )
    return theta, errors, negative_amps


def get_boundaries(
    z: np.ndarray, d: np.ndarray, trust_radius: float
) -> Tuple[float, float]:
    """Determine the intersections of the search direction and the trust
    region."""
    a = d.T @ d
    b = 2 * z.T @ d
    c = (z.T @ z) - (trust_radius ** 2)
    aux = b + np.copysign(
        np.sqrt(b * b - 4 * a * c),
        b,
    )
    return sorted([-aux / (2 * a), -(2 * c) / aux])
