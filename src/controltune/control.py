from __future__ import annotations

import numpy as np

UNSTABLE = 1e3


def closed_loop(
    a: np.ndarray, b: np.ndarray, kp: np.ndarray, ki: np.ndarray, n: int = 400, bound: float = 50.0
):
    """Unit-step setpoint response of the loop PI (velocity form) around the plant a(z^-1) y = b(z^-1) u, vectorised.

    a, b: (M, la) and (M, lb) polynomial coefficients (a[:, 0] = 1) for M models.
    kp, ki: (G,) controller gains. Returns y of shape (G, M, n) and u of shape (G, M, n).
    Plant output y is a deviation from the operating point; setpoint r = 1 from t = 0.
    """
    G, M = len(kp), len(a)
    la, lb = a.shape[1], b.shape[1]
    y = np.zeros((G, M, n + 1))
    u = np.zeros((G, M, n + 1))
    e_prev = np.zeros((G, M))
    for t in range(1, n + 1):
        acc = np.zeros((G, M))
        for i in range(1, la):
            if t - i >= 0:
                acc -= a[None, :, i] * y[:, :, t - i]
        for i in range(lb):
            if t - 1 - i >= 0:
                acc += b[None, :, i] * u[:, :, t - 1 - i]
        y[:, :, t] = np.clip(acc, -1e6, 1e6)
        e = 1.0 - y[:, :, t]
        u[:, :, t] = u[:, :, t - 1] + kp[:, None] * (e - e_prev) + ki[:, None] * e
        u[:, :, t] = np.clip(u[:, :, t], -1e6, 1e6)
        e_prev = e
    return y[:, :, 1:], u[:, :, 1:]


def cost(y: np.ndarray, u: np.ndarray, bound: float = 50.0, w_tv: float = 0.01):
    """IAE of the step response plus a small control-variation term; unstable / divergent loops get UNSTABLE."""
    iae = np.abs(1.0 - y).sum(axis=-1) / y.shape[-1]
    tv = np.abs(np.diff(u, axis=-1)).sum(axis=-1) / y.shape[-1]
    bad = (np.abs(y).max(axis=-1) > bound) | ~np.isfinite(iae)
    c = iae + w_tv * tv
    return np.where(bad, UNSTABLE, c), bad


def overshoot(y: np.ndarray) -> np.ndarray:
    return np.maximum(y.max(axis=-1) - 1.0, 0.0)
