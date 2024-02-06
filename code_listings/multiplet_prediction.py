from typing import Dict, Iterable
import numpy as np
np.set_printoptions(precision=3)


def multiplet_prediction(theta, epsilon) -> Iterable[Iterable[int]]:
    # Dictionary with central frequencies as keys, and lists of indices as values.
    multiplets = {}
    # Reshape for easier processing $\rightarrow \symbf{\theta} \in \mathbb{R}^{M \times 6}$
    theta = theta.reshape((theta.size // 6, 6))
    # Iterate over each oscillator index and parameters
    for idx, params in enumerate(theta):
        center = params[3] - params[2]  # $f_{\text{c}} = \ftwo - \fone$
        assigned = False
        for freq in multiplets.keys():
            # Test whether $\left\lvert \left(\ftwo_i - \fone_i\right) - \left(\ftwo_j - \fone_j\right) \right\rvert < \epsilon$
            if abs(center - freq) < epsilon:
                assigned = True
                multiplets[freq].append(idx)
                break
        if not assigned:
            multiplets[center] = [idx]

    return list(multiplets.values())

theta = np.array(
    [
        [4, 0, 0, 20, 3, 3],
        [2, 0, -7, -2, 3, 3],
        [2, 0, 7, 12, 3, 3],
        [1, 0, -10, -25, 3, 3],
        [2, 0, 0, -15, 3, 3],
        [1, 0, 10, -5, 3, 3],
    ],
    dtype="float64",
)
theta[:, 2:4] += np.random.uniform(-0.05, 0.05, (theta.shape[0], 2))
print("θ =\n{}".format(theta))
epsilon = 0.1
multiplets = multiplet_prediction(theta, epsilon)
print("Predicted multiplets:\n{}".format(multiplets))

# Output
# ======
# θ =
# [[ 4.000e+00  0.000e+00 -3.978e-02  1.997e+01  3.000e+00  3.000e+00]
#   [ 2.000e+00  0.000e+00 -7.044e+00 -2.044e+00  3.000e+00  3.000e+00]
#   [ 2.000e+00  0.000e+00  6.969e+00  1.198e+01  3.000e+00  3.000e+00]
#   [ 1.000e+00  0.000e+00 -9.983e+00 -2.502e+01  3.000e+00  3.000e+00]
#   [ 2.000e+00  0.000e+00 -9.125e-04 -1.495e+01  3.000e+00  3.000e+00]
#   [ 1.000e+00  0.000e+00  9.982e+00 -5.027e+00  3.000e+00  3.000e+00]]
# Predicted multiplets:
# [[0], [1, 2], [3, 4, 5]]
