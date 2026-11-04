# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

No linter is configured. Install with `pip install -e ".[test]"`. Runtime dependencies are `numpy`, `scipy`, `pandas`, `matplotlib`, and the tests need `pytest`.

- Run the tests: `python3 -m pytest`. This covers `tests/` and each catalog model's `test_validation.py`.
- Run one scenario: `python3 experiments/single_run.py --model <name> [--scenario <name>]`
- Run a perturbation study: `python3 experiments/perturbation_study.py --model <name>`
- Regenerate the catalog index after changing an entry: `python3 -m ode_catalog.index`. A test fails if it's out of date.

Experiment scripts write to `outputs/<experiment>/<model>/` by default, which git ignores. Use `--output-dir` when testing so you don't overwrite someone's results.

## Architecture

**`ode_solver/` is the engine and must stay model-agnostic.** It must never import `ode_catalog`, and `tests/test_catalog.py` enforces that.

- **`spec.py`** — the contract. `ODESystemSpec` is what every system implements: `state_names`, `param_names`, `rhs(t, y, **params)`, optional `input_names`, `conserved_total`, and defaults. `SimulationRequest` is one run's input and `SimulationResult` its output.
- **`solve.py`** — `validate_request(spec, request)` returns a list of problems, including unknown parameters. `solve_ode(spec, request)` wraps `scipy.integrate.solve_ivp` and returns a `SimulationResult`; check its `success` flag. Inputs left out of a request fall back to `spec.default_inputs` and are listed in `metadata["defaulted_inputs"]`.
- **`sweep.py`** — `generate_dataset(spec, request, ...)` runs many parameter sets and returns one long-format `DataFrame` (`run`, `t`, `success`, state columns, param columns). Failed runs are kept with `success=False` and trigger a warning, so filter on `success` before using the data.
- **`analysis.py`** — `analyze_result(result, spec)` computes per-state stats, plus conservation drift only when the spec declares `conserved_total`.
- **`plotting.py`** — `plot_results` and `plot_sweep` draw one subplot per state. They take `time_unit` and `state_units`, which `ModelEntry.units` provides.
- **`pipeline.py`** — `Pipeline.run(spec, request, ...)` runs validate, solve, analyze and plot in order, logging each step and stopping at the first failure.

**`ode_catalog/` is the model collection.** `schema.py` defines `ModelEntry`, which wraps an `ODESystemSpec` with `Parameter`, `Reference`, `Scenario` and `ModelStatus`. Each model has its own folder with `model.py`, `README.md` and `test_validation.py`. Follow `ADDING_A_MODEL.md` to add one. `INDEX.md` lists the registered models.

- `ModelStatus` uses the BioModels terms: executable, numerically verified, scientifically reproduced. Don't raise a status without the evidence that definition requires. Conservation checks or checking that equations were transcribed correctly aren't enough.
- `reference=None` means a model isn't from a publication, or that no source is recorded. Status and source are recorded separately.
- Record values and units as the source prints them, and put conversions or apparent errors in `Parameter.note`. Use `None` for units that aren't recorded. Don't invent citation details or units, and don't commit papers or their images.
- `validation.py` holds `assert_agrees`, `rk4_fixed_step`, and the figure-comparison tools. Only a person can record `overlay_reviewed_by`; don't fill it in or fabricate digitized data.

**`tests/`** — `test_model_contract.py` runs every model in `CATALOG` through the full pipeline. `test_catalog.py` checks entries and scenarios. `test_baseline.py` compares results against `tests/baseline/`, recorded before the restructuring; don't regenerate those files to make a failing test pass. `test_validation_tools.py` uses synthetic data only. `test_experiments.py` runs the scripts for every model.

This project directory lives under iCloud Desktop sync — rapid successive writes to the same filename (e.g. scripted `sed`+`mv`) can produce iCloud conflict copies (`file 2.py`) instead of overwriting in place.
