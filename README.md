# articles-ControlTune

Identification of a real plant from input–output data (DaISy steam generator, Abbott Power Plant: fuel → drum pressure, MISO ARX model)
and tuning of a discrete PI controller with and without model uncertainty (bootstrap and window-to-window variation).
**All closed-loop results are simulations on identified models; no controller was run on a real plant.** TRL 4 (simulation-based).

```
pip install -e ".[dev]" && python -m pytest
# data: DaISy steamgen.dat (https://ftp.esat.kuleuven.be/pub/SISTA/data/process_industry/steamgen.dat.gz; not committed)
python research/run_study.py <steamgen.dat>
```
