from __future__ import annotations

import numpy as np

NA, NB, DELAY = 3, 3, 1


def regressors(y: np.ndarray, U: np.ndarray, na: int = NA, nb: int = NB, delay: int = DELAY):
    """Rows x_t = [-y_{t-1..t-na}, u_j(t-delay-0..nb-1) for every input j, 1] and targets y_t."""
    t0 = max(na, nb + delay - 1 + 1)
    X, T = [], []
    for t in range(t0, len(y)):
        row = [-y[t - i] for i in range(1, na + 1)]
        for j in range(U.shape[1]):
            row += [U[t - delay - i, j] for i in range(nb)]
        row.append(1.0)
        X.append(row)
        T.append(y[t])
    return np.array(X), np.array(T)


def fit_arx(X: np.ndarray, T: np.ndarray, ridge: float = 1e-8) -> np.ndarray:
    A = X.T @ X + ridge * np.eye(X.shape[1])
    return np.linalg.solve(A, X.T @ T)


def simulate(
    theta: np.ndarray,
    y_init: np.ndarray,
    U: np.ndarray,
    n_steps: int,
    na: int = NA,
    nb: int = NB,
    delay: int = DELAY,
):
    """Free-run simulation from ``y_init`` (last na outputs, oldest first) driven by the inputs U (rows 0.. aligned with the steps)."""
    y = list(y_init)
    for k in range(n_steps):
        row = [-y[-i] for i in range(1, na + 1)]
        for j in range(U.shape[1]):
            row += [U[k - delay - i, j] if k - delay - i >= 0 else 0.0 for i in range(nb)]
        row.append(1.0)
        y.append(float(np.dot(theta, row)))
    return np.array(y[len(y_init) :])


def siso_polynomials(
    theta: np.ndarray, input_index: int, na: int = NA, nb: int = NB, delay: int = DELAY
):
    """Denominator a (1 + a1 z^-1 + ...) and numerator b (delay included) of the transfer function from input ``input_index``."""
    a = np.r_[1.0, theta[:na]]
    b = np.r_[np.zeros(delay), theta[na + input_index * nb : na + (input_index + 1) * nb]]
    return a, b
