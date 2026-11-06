import numpy as np

from controltune.control import closed_loop, cost, overshoot
from controltune.ident import fit_arx, regressors, simulate, siso_polynomials


def _plant(n=3000, seed=0):
    rng = np.random.default_rng(seed)
    u = rng.normal(0, 1, (n, 1))
    y = np.zeros(n)
    for t in range(3, n):  # y_t = 0.7 y_{t-1} + 0.5 u_{t-1} + 0.2 u_{t-2} + noise + 1 (offset)
        y[t] = (
            0.7 * y[t - 1]
            + 0.5 * u[t - 1, 0]
            + 0.2 * u[t - 2, 0]
            + 1.0
            + 0.01 * rng.standard_normal()
        )
    return y, u


def test_arx_recovers_known_plant_and_polynomials():
    y, u = _plant()
    X, T = regressors(y, u)
    th = fit_arx(X, T)
    a, b = siso_polynomials(th, 0)
    assert abs(a[1] + 0.7) < 0.02 and abs(b[1] - 0.5) < 0.02 and abs(b[2] - 0.2) < 0.02
    ysim = simulate(th, y[2000:2003], u[2000:2400], 300)
    assert np.isfinite(ysim).all()


def test_closed_loop_pi_reaches_setpoint_and_flags_unstable_gains():
    a = np.array([[1.0, -0.7, 0.0, 0.0]])
    b = np.array([[0.0, 0.5, 0.2, 0.0]])
    y, u = closed_loop(a, b, np.array([0.5, 50.0]), np.array([0.2, 20.0]), n=300)
    c, bad = cost(y, u)
    assert not bad[0, 0] and abs(y[0, 0, -1] - 1.0) < 0.05
    assert bad[1, 0] and c[1, 0] == 1e3
    assert overshoot(y)[0, 0] >= 0.0
