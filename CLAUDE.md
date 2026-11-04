# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

No linter is configured. `pyproject.toml` defines the `ode-solver` distribution; install it with `pip install -e ".[test]"`. Runtime dependencies are `numpy`, `scipy`, `pandas`, `matplotlib`, and the tests need `pytest`.

- Run the tests: `python3 -m pytest`. This covers `tests/` and each catalog model's `test_validation.py`.
- Run a single-simulation demo (bioprocess model, prints final state, shows a plot): `python3 single_run_demo.py`
- Run the perturbation study (200 runs, writes `perturbation_results.csv` and `.png` to the current folder, overwriting any previous ones): `python3 perturbation_study.py`
- Switch which ODE system the study runs against: edit `MODEL_NAME` at the top of `perturbation_study.py`.

## Architecture

**`ode_solver/` is the engine and must stay model-agnostic.** It must never import `ode_catalog` or `models`, and `tests/test_catalog.py` enforces that. All system-specific detail lives in `ode_catalog/` or `models/`.

- **`spec.py`** — the contract. `ODESystemSpec` is what every system implements: `state_names`, `param_names`, `rhs(t, y, **params)`, optional `input_names`, `conserved_total`, and defaults. `SimulationRequest` is one run's input and `SimulationResult` its output.
- **`solve.py`** — `validate_request(spec, request)` returns a list of problems, including unknown parameters. `solve_ode(spec, request)` wraps `scipy.integrate.solve_ivp` and returns a `SimulationResult`; check its `success` flag. Inputs left out of a request fall back to `spec.default_inputs` and are listed in `metadata["defaulted_inputs"]`.
- **`sweep.py`** — `generate_dataset(spec, request, ...)` runs many parameter sets and returns one long-format `DataFrame` (`run`, `t`, `success`, state columns, param columns). Failed runs are kept with `success=False` and trigger a warning, so filter on `success` before using the data.
- **`analysis.py`** — `analyze_result(result, spec)` computes per-state stats, plus conservation drift only when the spec declares `conserved_total`.
- **`plotting.py`** — `plot_results(result, spec, save_path=None, show=False)` draws one subplot per state.
- **`pipeline.py`** — `Pipeline.run(spec, request, ...)` runs validate, solve, analyze and plot in order, logging each step and stopping at the first failure.

**`ode_catalog/` is the model collection.** `schema.py` defines `ModelEntry`, which wraps an `ODESystemSpec` with `Parameter`, `Reference`, `Scenario` and `ModelStatus`. Each model has its own folder with `model.py`, `README.md`, `reference_data/` and `test_validation.py`. Follow `ADDING_A_MODEL.md` to add one. Registered: `brain_heart`.

- `ModelStatus` uses the BioModels terms: executable, numerically verified, scientifically reproduced. Don't raise a status without the evidence that definition requires. Checking that equations were transcribed correctly isn't enough.
- `reference=None` means a model isn't from a publication. Status and source are recorded separately.
- Record values and units as the source prints them, and put conversions or apparent errors in `Parameter.note`. Don't invent citation details, and don't commit papers or their images.

**`models/` is temporary**, holding systems not yet migrated to the catalog: `bioprocess`, `metabolite_pathway`, `sird`. It's left out of the installed package.

**`tests/`** — `test_model_contract.py` runs every model in `MODELS` and `CATALOG` through the full pipeline. `test_catalog.py` checks entries and scenarios. `test_baseline.py` compares results against `tests/baseline/`, recorded before the restructuring; don't regenerate those files to make a failing test pass.

`perturbation_results.csv`/`.png` are generated artifacts, not source of truth — they reflect whichever `MODEL_NAME` was active the last time `perturbation_study.py` was run, and go stale relative to the code whenever the model or params change without a rerun.

This project directory lives under iCloud Desktop sync — rapid successive writes to the same filename (e.g. scripted `sed`+`mv`) can produce iCloud conflict copies (`file 2.py`) instead of overwriting in place.
