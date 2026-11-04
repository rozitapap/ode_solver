# ODE

A model-agnostic ODE engine and a catalog of mathematical models. Define a system of ODEs once, then integrate it, sweep its parameters, and collect the results for downstream analysis or modeling. Catalog entries also record each model's biological context, sources, parameter provenance, runnable scenarios, and validation status.

## Requirements

- Python 3.10+
- `numpy`, `scipy`, `pandas`, `matplotlib`

```bash
pip install -e .            # from this folder; add ".[test]" to include pytest
```

## Usage

```python
from ode_catalog import CATALOG
from ode_solver import Pipeline

entry = CATALOG["brain_heart"]
run = Pipeline().run(entry.spec, entry.request("setpoint_step"), plot_path="brain_heart.png")
```

Run a single simulation and see the resulting trajectories:

```bash
python3 single_run_demo.py
```

Run a full perturbation study — 200 simulations with the base parameters each randomly jittered by ±10%, saved to CSV and plotted as an overlay:

```bash
python3 perturbation_study.py
```

This writes `perturbation_results.csv` (one row per run per time point, with a `success` column flagging failed runs) and `perturbation_results.png` (all trajectories overlaid, per state).

To run the study against a different registered ODE system, edit `MODEL_NAME` at the top of `perturbation_study.py`.

## Structure

```
ode_solver/             the engine, with no knowledge of any specific model
    spec.py             ODESystemSpec, SimulationRequest, SimulationResult
    solve.py            validate_request, solve_ode (scipy.integrate.solve_ivp)
    sweep.py            generate_dataset: many parameter sets, long-format DataFrame
    analysis.py         per-state summary statistics, declared conservation checks
    plotting.py         plot_results
    pipeline.py         validate -> solve -> analyze -> plot, with logging
ode_catalog/            the model collection
    schema.py           ModelEntry, Parameter, Reference, Scenario, ModelStatus
    ADDING_A_MODEL.md   steps and checklist for a new model
    _template/          copy this to add a model
    brain_heart/        equations and sources, README, reference data, validation tests
models/                 models not yet migrated to the catalog
tests/                  contract tests that run every model, and the regression baseline
single_run_demo.py, perturbation_study.py   driver scripts
```

The catalog imports the engine, but the engine never imports the catalog.

## Validation status

Catalog entries use the BioModels empirical quality control terms:

| Status | Definition |
|---|---|
| executable | Simulation runs without crashing |
| numerically verified | Results agree across independent implementations or against a trusted reference within defined tolerances |
| scientifically reproduced | Simulation reproduces a documented, published result |

## Tests

```bash
python3 -m pytest
```

## Adding a model

Follow [ode_catalog/ADDING_A_MODEL.md](ode_catalog/ADDING_A_MODEL.md). The tests cover a new model automatically once it's in `CATALOG`.
