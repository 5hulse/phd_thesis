# mmempm.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Fri 04 Aug 2023 14:43:44 BST

def mmempm(
    Y: np.ndarray, sw1: float, sw2: float,
    offset1: float, offset2: float, M: int = 0,
) -> np.ndarray:
    """Modified Matrix Enhancment Matrix Pencil Method for estimating a
   2D FID.

    Parameters
    ----------
    Y
        FID.

    sw1
        Sweep width in first (indirect) dimension (Hz).

    sw2
        Sweep width in second (direct) dimension (Hz).

    offset1
        Transmitter offset in first dimension (Hz).

    offset2
        Transmitter offset in second dimension (Hz).

    M
        Number of oscillators. If 0, this is estimated by applying the MDL
        to the first direct dimension FID (i.e. Y[0])
    """
    N1, N2 = Y.shape  # $N^{(1)}, N^{(2)}$
    L1, L2 = N1 // 2, N2 // 2  # $L^{(1)}, L^{(2)}$: pencil parameters
    norm = np.linalg.norm(Y)  # $\left\lVert\symbf{Y}\right\rVert$
    Y /= norm  # $\nicefrac{\symbf{Y}}{\left\lVert\symbf{Y}\right\rVert}$

    if M == 0:
        # Extract first direct dimension signal, and estimate model order
        # with the MDL
        L_mdl = N2 // 3
        Y_mdl = Y[0]
        col_mdl = Y_mdl[: N2 - L_mdl]
        row_mdl = Y_mdl[N2 - L_mdl - 1 :]
        HY_mdl = sp.linalg.hankel(col_mdl, row_mdl)
        _, sigma_mdl, _ = np.linalg.svd(HY_mdl)
        M = mdl(sigma_mdl, N2, L_mdl)

    # === Construct block Hankel $\symbf{E}_{\symbf{Y}}$ ===
    row_size = L2
    col_size = N2 - L2 + 1
    EY = np.zeros(
        (L1 * L2, (N1 - L1 + 1) * (N2 - L2 + 1)),
        dtype="complex",
    )
    for n1 in range(N1):
        # Construct $\symbf{H}_{\symbf{Y},n^{(1)}}$, and assign to appropriate positions in $\symbf{E}_{\symbf{Y}}$
        col = Y[n1, :L2]
        row = Y[n1, L2 - 1:]
        HYn1 = sp.linalg.hankel(col, row)
        for n in range(n1 + 1):
            r, c = n, n1 - n
            if r < L1 and c < N1 - L1 + 1:
                EY[
                    r * row_size : (r + 1) * row_size,
                    c * col_size : (c + 1) * col_size
                ] = HYn1

    EY = sp.sparse.csr_matrix(EY)  # Make $\symbf{E}_{\symbf{Y}}$ sparse $\label{ln:sparse1}$
    UM, *_ = sp.sparse.linalg.svds(EY, k=M)  # $\symbf{U}_M$ $\label{ln:sparse2}$

    # === Construct permutation matrix $\symbf{P}$ ===
    P = np.zeros((L1 * L2, L1 * L2))
    r = 0
    for l2 in range(L2):
        for l1 in range(L1):
            c = l1 * L2 + l2
            P[r, c] = 1
            r += 1

    UM1 = UM[: L1 * (L2 - 1)]  # Last $L^{(2)}$ rows deleted: $\symbf{U}_{M1}$
    UM2 = UM[L2:]  # First $L^{(2)}$ rows deleted: $\symbf{U}_{M2}$
    z1, W1 = np.linalg.eig(np.linalg.pinv(UM1) @ UM2)  # $\symbf{z}^{(1)}$, $\symbf{W}^{(1)}$

    UMP = P @ UM  # $\symbf{U}_{MP}$
    UMP1 = UMP[: (L1 - 1) * L2]  # Last $L^{(1)}$ rows deleted: $\symbf{U}_{MP1}$
    UMP2 = UMP[L1:]  # First $L^{(1)}$ rows deleted: $\symbf{U}_{MP2}$
    G = np.linalg.inv(W1) @ np.linalg.pinv(UMP1) @ UMP2 @ W1  # $\symbf{G}$
    z2 = np.diag(G).copy()  # $\symbf{z}^{(2)}$: copy needed as slice is readonly

    # === Check for and deal with similar frequencies in $\symbf{f}^{(1)}$ ===
    freq1 = (0.5 * sw1 / np.pi) * np.imag(np.log(z1)) + offset1  # $\symbf{f}^{(1)} \label{ln:similar-f-start}$
    threshold = sw1 / N1  # $\nicefrac{f_{\text{sw}}^{(1)}}{N^{(1)}}$
    groupings = {}
    # Iterate through values in $\symbf{f}^{(1)}$ and group any with
    # similar frequencies together
    for idx, f1 in enumerate(freq1):
        assigned = False
        for group_f1, indices in groupings.items():
            if np.abs(f1 - group_f1) < threshold:
                indices.append(idx)
                n = len(indices)
                indices = sorted(indices)
                # Get new mean freq of the group
                new_group_f1 = (n * group_f1 + f1) / (n + 1)
                groupings[new_group_f1] = groupings.pop(group_f1)
                assigned = True
                break
        if not assigned:
            groupings[f1] = [idx]

    for indices in groupings.values():
        n = len(indices)
        if n != 1:
            Gr_slice = tuple(zip(*product(indices, repeat=2)))
            Gr = G[Gr_slice].reshape(n, n)
            new_group_z2, _ = np.linalg.eig(Gr)
            z2[indices] = new_group_z2  # $\label{ln:similar-f-end}$

    # === Construct $\symbf{E}_{\text{L}}$ and $\symbf{E}_{\text{R}}$ ===
    ZL2 = np.power.outer(z2, np.arange(L2)).T  # $\symbf{Z}_{\text{L}}^{(2)}$
    ZR2 = np.power.outer(z2, np.arange(N2 - L2 + 1))  # $\symbf{Z}_{\text{R}}^{(2)}$
    Z1D = np.diag(z1)  # $\symbf{Z}_{\text{D}}^{(1)}$

    EL = np.zeros((L1 * L2, M), dtype="complex")
    Z2LZ1D = ZL2
    for i in range(L1):
        EL[i * row_size : (i + 1) * row_size] = Z2LZ1D
        Z2LZ1D = Z2LZ1D @ Z1D
    ER = np.zeros((M, (N1 - L1 + 1) * (N2 - L2 + 1)), dtype="complex")
    Z1DZ2R = ZR2
    for i in range(N1 - L1 + 1):
        ER[:, i * col_size : (i + 1) * col_size] = Z1DZ2R
        Z1DZ2R = Z1D @ Z1DZ2R
    # $\symbf{\alpha} = \operatorname{diag}\left(\symbf{E}_{\text{L}}^+ \symbf{E}_{\symbf{Y}} \symbf{E}_{\text{R}}^+\right)$
    alpha = np.diag(np.linalg.pinv(EL) @ EY @ np.linalg.pinv(ER))

    amp = np.abs(alpha) * norm  # $\symbf{a}$
    phase = np.arctan2(np.imag(alpha), np.real(alpha))  # $\symbf{\phi}$
    freq2 = (0.5 * sw2 / np.pi) * np.imag(np.log(z2)) + offset2  # $\symbf{f}^{(2)}$
    damp1 = -sw1 * np.real(np.log(z1))  # $\symbf{\eta}^{(1)}$
    damp2 = -sw2 * np.real(np.log(z2))  # $\symbf{\eta}^{(2)}$
    # $\symbf{\theta}$, as a $M \times 6$ array
    theta = np.vstack((amp, phase, freq1, freq2, damp1, damp2)).T

    # Remove negative damping factors
    neg_damp_idx1 = list(np.nonzero(damp1 < 0.)[0])
    neg_damp_idx2 = list(np.nonzero(damp2 < 0.)[0])
    neg_damp_idx = list(set(neg_damp_idx1 + neg_damp_idx2))
    theta = np.delete(theta, neg_damp_idx, axis=0)

    theta = theta[np.argsort(theta[:, 3])]  # order by $\symbf{f}^{(2)}$
    return theta
