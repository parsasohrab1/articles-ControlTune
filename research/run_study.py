"""DaISy steam generator (Abbott Power Plant, real data, 9600 samples at 3 s): identification and uncertainty-aware PI tuning.

usage: python research/run_study.py <steamgen.dat>
Plant: drum pressure (output 1) from fuel (input 1); the other three inputs are part of the identified MISO ARX(3,3) model and are
held at their operating values in the control simulation. Identification on the first 70 % of the record, validation (free-run
simulation fit) on the last 30 %. Uncertainty: block bootstrap of the regression rows (block 100).
Tuning (grid search over normalised PI gains) on models identified from the first 70 %:
  nominal - best cost on the nominal model; robust - best 90th-percentile cost over a 100-model bootstrap ensemble.
Evaluation: (i) the tuning ensemble, (ii) an independent ensemble identified from the LAST 30 % of the record (plant drift).
All closed-loop numbers are SIMULATIONS on identified models: no controller was run on the real plant.
"""

import json
import sys
from pathlib import Path

import numpy as np

from controltune.control import closed_loop, cost, overshoot
from controltune.ident import NA, fit_arx, regressors, simulate, siso_polynomials

M = 100
BLOCK = 100
STEPS = 400


def boot_models(X, T, rng, m=M):
    n = len(T)
    thetas = []
    for _ in range(m):
        idx = np.concatenate(
            [
                np.arange(s, min(s + BLOCK, n))
                for s in rng.integers(0, n - BLOCK, int(np.ceil(n / BLOCK)))
            ]
        )[:n]
        thetas.append(fit_arx(X[idx], T[idx]))
    return np.array(thetas)


def polys(thetas):
    a, b = zip(*[siso_polynomials(t, 0) for t in thetas], strict=True)
    return np.array(a), np.array(b)


def main(path: str):
    d = np.loadtxt(path)
    U, y = d[:, 1:5], d[:, 5]
    n = len(y)
    ntr = int(0.7 * n)
    X, T = regressors(y, U)
    t0 = len(y) - len(T)
    split = ntr - t0
    th = fit_arx(X[:split], T[:split])
    # validation: free-run simulation on the last 30 %
    yv = y[ntr - NA : ntr]
    ysim = simulate(
        th, yv, U[ntr - 10 :], n - ntr
    )  # U shifted by 10 rows so that delays refer to real inputs
    # simpler and exact: re-simulate with aligned inputs
    ysim = np.array(
        simulate(th, y[ntr - NA : ntr], np.vstack([U[ntr - 3 : ntr], U[ntr:]]), n - ntr)
    )
    ytrue = y[ntr:]
    fit = 100 * (1 - np.linalg.norm(ytrue - ysim) / np.linalg.norm(ytrue - ytrue.mean()))
    a0, b0 = siso_polynomials(th, 0)
    g0 = b0.sum() / a0.sum()
    out = {
        "note": "real DaISy steam-generator data; closed-loop results are simulations on identified models",
        "n": n,
        "validation_simulation_fit_percent": float(fit),
        "dc_gain_fuel_to_pressure": float(g0),
        "nominal_poles_abs": np.sort(np.abs(np.roots(a0))).tolist(),
    }
    rng = np.random.default_rng(0)
    tune_th = boot_models(X[:split], T[:split], rng)
    late_th = boot_models(X[split:], T[split:], rng)
    # structural / operating-condition variation: models identified separately on 1200-sample windows of the record
    win = 1200
    wins = [fit_arx(X[i : i + win], T[i : i + win]) for i in range(0, len(T) - win + 1, win)]
    half = len(wins) // 2
    early_w, late_w = np.array(wins[:half]), np.array(wins[half:])
    out["window_models"] = {
        "n_early": len(early_w),
        "n_late": len(late_w),
        "dc_gain_early": [
            float(siso_polynomials(t, 0)[1].sum() / siso_polynomials(t, 0)[0].sum())
            for t in early_w
        ],
        "dc_gain_late": [
            float(siso_polynomials(t, 0)[1].sum() / siso_polynomials(t, 0)[0].sum()) for t in late_w
        ],
    }
    a_t, b_t = polys(
        np.vstack([tune_th, early_w])
    )  # tuning ensemble = bootstrap models + early windows
    a_l, b_l = polys(
        np.vstack([late_th, late_w])
    )  # evaluation ensemble = late bootstrap models + late windows
    a_w, b_w = polys(late_w)
    a_n, b_n = a0[None], b0[None]
    kp_n = np.logspace(-2, 3, 41)
    ki_n = np.logspace(-3, 1.5, 41)
    KP, KI = np.meshgrid(kp_n, ki_n, indexing="ij")
    kp, ki = KP.ravel() / abs(g0), KI.ravel() / abs(g0)
    # grid evaluation on the nominal model and on the tuning ensemble
    yn, un = closed_loop(a_n, b_n, kp, ki, STEPS)
    cn, _ = cost(yn, un)
    ye, ue = closed_loop(a_t, b_t, kp, ki, STEPS)
    ce, bad_e = cost(ye, ue)
    best_nom = int(np.argmin(cn[:, 0]))
    robust_obj = np.quantile(ce, 0.9, axis=1)
    best_rob = int(np.argmin(robust_obj))
    out["tuned_gains_normalised"] = {
        "nominal": [float(KP.ravel()[best_nom]), float(KI.ravel()[best_nom])],
        "robust": [float(KP.ravel()[best_rob]), float(KI.ravel()[best_rob])],
    }
    # evaluation sets
    sets = {
        "tuning_ensemble": (a_t, b_t),
        "late_ensemble_plant_drift": (a_l, b_l),
        "late_windows_only": (a_w, b_w),
        "nominal_model": (a_n, b_n),
    }
    out["evaluation"] = {}
    for label, idx in (("nominal_tuning", best_nom), ("robust_tuning", best_rob)):
        out["evaluation"][label] = {}
        for sname, (a, b) in sets.items():
            yy, uu = closed_loop(a, b, kp[[idx]], ki[[idx]], STEPS)
            c, bad = cost(yy, uu)
            good = ~bad[0]
            out["evaluation"][label][sname] = {
                "share_unstable": float(bad[0].mean()),
                "median_cost": float(np.median(c[0])),
                "p90_cost": float(np.quantile(c[0], 0.9)),
                "median_overshoot_stable": float(np.median(overshoot(yy)[0][good]))
                if good.any()
                else None,
            }
    # how many grid points are stable for all late models
    yl, ul = closed_loop(a_l, b_l, kp, ki, STEPS)
    _, bad_l = cost(yl, ul)
    out["grid"] = {
        "share_gain_pairs_stable_on_nominal": float((cn[:, 0] < 1e3).mean()),
        "share_gain_pairs_stable_on_all_tuning_models": float((~bad_e).all(axis=1).mean()),
        "share_gain_pairs_stable_on_all_late_models": float((~bad_l).all(axis=1).mean()),
    }
    Path(__file__).with_name("results.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main(sys.argv[1])
