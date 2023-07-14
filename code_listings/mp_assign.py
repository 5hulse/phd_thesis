def predict_multiplets(
    params: np.ndarray,
    thold: float,
) -> Dict[float, Iterable[int]]:
    """
    Parameters
    ----------
    params
        Estimated parameter array with shape (M, 6) such that each row
        provides the parameters of a particular oscillator, in the order
        [a, φ, f1, f2, η1, η2].

    thold
        Frequency threshold ε > 0.

    Returns
    -------
    A dictionary with the multiplet central frequencies as keys and
    oscillator indices as values
    """
    multiplets = {}
    for m, osc in enumerate(params):
        assigned = False
        f1, f2 = osc[2], osc[3]  # Extract $f^{(1)}$ and $f^{(2)}$
        fc = osc[3] - osc[2]  # Central frequency for oscillator: $f^{(2)}_m - f^{(1)}_m$
        # Check whether central frequency agrees with any
        # alread-established multiplet grouping
        for fmp in multiplets:
            if fmp - thold < fc < fmp + thold:
                # Update grouping:
                # Add m to list of oscillators
                # Update central frequency: mean of all oscillators
                ms = multiplets.pop(fmp)
                ms.append(m)
                mp_size = len(ms)
                new_fmp = ((mp_size - 1) * fmp + fc) / mp_size
                multiplets[new_fmp] = ms
                assigned = True
                break
        # No match: create new multiplet group
        if not assigned:
            multiplets[fc] = [m]

    return multiplets
