# Software Requirements Specification – ControlTune

- SR-1 `ident.regressors`/`fit_arx`: MISO ARX(3,3) with one-sample delay and an offset, least squares; `simulate`: free-run simulation; `siso_polynomials`: denominator and numerator for one input.
- SR-2 `control.closed_loop`: unit-step setpoint response of a velocity-form PI loop around the identified polynomials, vectorised over models and gain pairs; `cost`: mean absolute error plus 0.01 x control variation, 1000 for divergent loops (|y| > 50); `overshoot`.
- SR-3 `research/run_study.py`: identification on the first 70 % of the record, free-run validation on the last 30 %; block-bootstrap ensembles (100 models, block 100) of the first 70 % and the last 30 %; additional models from 1200-sample windows; grid search of normalised gains (41 x 41); nominal tuning = best cost on the nominal model; robust tuning = best 90th-percentile cost over the tuning ensemble (bootstrap models and early windows).
- SR-4 Two unit tests (ARX recovers a known plant; PI loop converges and divergent gains are flagged).
