# ODE

A model-agnostic ODE engine and a catalog of mathematical models. Define a system of ODEs once, then integrate it, sweep its parameters, and collect the results for downstream analysis or modeling. Catalog entries also record each model's biological context, sources, parameter provenance, runnable scenarios, and validation status.

See the [model catalog](ode_catalog/INDEX.md) for what's included.

## Requirements

Python 3.10+. Install from this folder:

```bash
pip install -e ".[test]"
```

This installs `numpy`, `scipy`, `pandas` and `matplotlib`, plus `pytest` for the tests.

To use it from another project, once the `v0.1.0` tag has been pushed:

```bash
pip install "ode-solver @ git+https://github.com/rozitapap/ode_solver.git@v0.1.0"
```

## Usage

```python
from ode_catalog import CATALOG
from ode_solver import Pipeline

entry = CATALOG["brain_heart"]
run = Pipeline().run(entry.spec, entry.request("setpoint_step"), plot_path="brain_heart.png", **entry.units)
```

Or run a scenario from the command line:

```bash
python3 experiments/single_run.py --model brain_heart --scenario setpoint_step
```

Run a perturbation study, with 200 simulations and each parameter randomly scaled by up to ±10%:

```bash
python3 experiments/perturbation_study.py --model sird
```

Results go to `outputs/<experiment>/<model>/`. In the sweep's `trajectories.csv`, the `success` column flags failed runs, so filter on it before using the data.

## Structure

```
ode_solver/             the engine, with no knowledge of any specific model
    spec.py             ODESystemSpec, SimulationRequest, SimulationResult
    solve.py            validate_request, solve_ode (scipy.integrate.solve_ivp)
    sweep.py            generate_dataset: many parameter sets, long-format DataFrame
    analysis.py         per-state summary statistics, declared conservation checks
    plotting.py         plot_results, plot_sweep
    pipeline.py         validate -> solve -> analyze -> plot, with logging
ode_catalog/            the model collection
    INDEX.md            generated list of models, sources and statuses
    schema.py           ModelEntry, Parameter, Reference, Scenario, ModelStatus
    validation.py       tools for verification and figure comparison
    ADDING_A_MODEL.md   steps and checklist for a new model
    _template/          copy this to add a model
    <model>/            equations and sources, README, validation tests
experiments/            scripts that run catalog models
tests/                  contract tests that run every model, and the regression baseline
outputs/                generated results, not tracked by git
```

The catalog imports the engine, but the engine never imports the catalog.

## Validation status

Catalog entries use the BioModels empirical quality control terms:

| Status | Definition |
|---|---|
| executable | Simulation runs without crashing |
| numerically verified | Results agree across independent implementations or against a trusted reference within defined tolerances |
| scientifically reproduced | Simulation reproduces a documented, published result |

No model is scientifically reproduced yet.

## Tests

```bash
python3 -m pytest
```

## Adding a model

Follow [ode_catalog/ADDING_A_MODEL.md](ode_catalog/ADDING_A_MODEL.md). The tests cover a new model automatically once it's in `CATALOG`.
