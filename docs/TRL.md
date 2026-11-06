# TRL assessment – ControlTune

**Claimed level: TRL 4** (simulation-based, on models identified from public real plant data). Not TRL 5: the controller was never applied to the plant or a hardware-in-the-loop setup.

Data: DaISy steam generator, 9600 samples (3 s). Model of drum pressure from the four inputs (fuel, air, reference level, disturbance), ARX(3,3). Free-run simulation fit on the last 30 % of the record: **73 %**. Fuel-to-pressure DC gain 318; poles at 0.25, 0.60, 0.98 (one slow pole close to 1).

Window models (1200 samples each): DC gain early 292, 320, 302; late 315, 378, 346, 309.

Tuned normalised gains (Kp x gain, Ki x gain): nominal (10.0, 0.23); robust (10.0, 0.23).

| Tuning | Evaluated on | Share unstable | Median cost | 90th-percentile cost | Median overshoot (stable) |
|---|---|---|---|---|---|
| nominal tuning | tuning ensemble | 0.00 | 0.0168 | 0.0176 | 0.20 |
| nominal tuning | late ensemble plant drift | 0.00 | 0.0184 | 0.0196 | 0.28 |
| nominal tuning | late windows only | 0.00 | 0.0176 | 0.0183 | 0.21 |
| nominal tuning | nominal model | 0.00 | 0.0168 | 0.0168 | 0.20 |
| robust tuning | tuning ensemble | 0.00 | 0.0168 | 0.0176 | 0.20 |
| robust tuning | late ensemble plant drift | 0.00 | 0.0184 | 0.0196 | 0.28 |
| robust tuning | late windows only | 0.00 | 0.0176 | 0.0183 | 0.21 |
| robust tuning | nominal model | 0.00 | 0.0168 | 0.0168 | 0.20 |

Share of the 41 x 41 gain pairs that are stable on the nominal model 0.46; on all tuning-ensemble models 0.44; on all late models 0.43.

## Findings and caveats
- **The uncertainty-aware tuning is identical to the nominal tuning** on this plant: the cost-optimal gains lie well inside the stable region of every model, and the identified models differ little (gain varies by about ±15 %). No benefit of the robust variant is demonstrated; it would matter only for plants whose uncertainty reaches the stability boundary.
- The cost-optimal gains give a median overshoot of about 20–28 % in simulation: the criterion (absolute error) favours speed over damping; no overshoot constraint was used.
- The identified model is only moderately accurate (simulation fit 73 %); the single-loop design ignores the coupling with the other outputs and the disturbance input, which were held constant in the simulation. These unmodelled effects are the main risk and are not covered by the bootstrap or window ensembles.
- About 55 % of the gain grid is unstable on every model — so tuning by trial on the plant would need care; the stable region itself is stable across models.
- Stationarity: the late-window models have a higher gain than the early ones (up to 378 against 320), but the closed loop remained stable.
