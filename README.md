# ODE

A small, model-agnostic ODE simulation pipeline: define a system of ODEs once, then integrate it, sweep its parameters across many runs, and collect the results into a dataset for downstream analysis or modeling.

## Requirements

- Python 3
- `numpy`, `scipy`, `pandas`, `matplotlib`

```bash
pip install numpy scipy pandas matplotlib
```

## Usage

Run a single simulation and see the resulting trajectories:

```bash
python3 ODE_solver.py
```

Run a full perturbation study — 200 simulations with the base parameters each randomly jittered by ±10%, saved to CSV and plotted as an overlay:

```bash
python3 perturbation_study.py
```

This writes `perturbation_results.csv` (one row per run per time point, with a `success` column flagging failed runs) and `perturbation_results.png` (all trajectories overlaid, per state).

To run the study against a different registered ODE system, edit `MODEL_NAME` at the top of `perturbation_study.py`.

## Structure

```
ODE_solver.py       generic validation, integrator (scipy.integrate.solve_ivp) and plotting, no system-specific knowledge
data_generator.py   runs an ODE model many times over sampled parameter sets, returns a long-format DataFrame
perturbation_study.py   driver script: random-perturbation study, saves CSV + PNG
models/
    base.py          ODESystemSpec, SimulationRequest, SimulationResult -- the interface every system implements
    bioprocess.py    5-state cell-culture kinetic model (VCD, Glc, Gln, Lac, Amm)
    sird.py          4-state SIRD epidemic model (S, I, R, D)
    __init__.py      MODELS registry mapping name -> ODESystemSpec
tests/
    test_model_contract.py   runs every registered model through the full pipeline
```

## Tests

```bash
pip install pytest
python3 -m pytest
```

## Adding a new ODE system

Add a new file to `models/` defining an `ODESystemSpec` (RHS function, state names, parameter names, defaults) and register it in `models/__init__.py`. No changes are needed to the solver or dataset generator. See `models/CLAUDE.md` for the exact steps and conventions.
