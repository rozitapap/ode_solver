# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

No build step or linter is configured in this repo. No dependency manifest exists either (no `requirements.txt`/`pyproject.toml`) — runtime dependencies are `numpy`, `scipy`, `pandas`, `matplotlib`, and the tests need `pytest`.

- Run the tests: `python3 -m pytest`. `tests/test_model_contract.py` runs every model in `MODELS` through the full pipeline, so a new model is covered once registered.
- Run a single-simulation demo (bioprocess model, prints final state, shows a plot): `python3 ODE_solver.py`
- Run the perturbation study (200 runs, writes `perturbation_results.csv` and `.png`): `python3 perturbation_study.py`
- Switch which ODE system the study runs against: edit `MODEL_NAME` at the top of `perturbation_study.py` (currently `"bioprocess"`; `"sird"` is also registered).

## Architecture

This codebase is ODE-system-agnostic: `ODE_solver.py` and `data_generator.py` contain no knowledge of any specific system. All system-specific detail (state names, parameter names, the RHS equations, default values) lives in `models/`.

- **`models/base.py`** — the shared contract, byte-identical to `ode_agent/schema.py` in the PINNSLearning repo (change both together). `ODESystemSpec` is what every system implements: `state_names`, `param_names`, `rhs(t, y, **params)`, and defaults for params/initial conditions/time grid. `SimulationRequest` is one run's input and `SimulationResult` its output.
- **`models/__init__.py`** — `MODELS` dict registers each system by name (`"bioprocess"`, `"sird"`). This is the single extension point: adding a new ODE system means adding `models/<name>.py` with one `ODESystemSpec` instance and one registry entry — no changes needed to the solver or generator.
- **`ODE_solver.py`** — `validate_request(spec, request)` returns a list of problems, including unknown parameters. `solve_ode(spec, request)` is the only integration entry point, wrapping `scipy.integrate.solve_ivp` and returning a `SimulationResult`; check its `success` flag. `plot_results(result, spec)` draws one subplot per `spec.state_names`.
- **`data_generator.py`** — `generate_dataset(spec, request, ...)` calls `solve_ode()` across many parameter sets and returns one long-format `DataFrame` (`run`, `t`, `success`, state columns, param columns). Failed runs are kept with `success=False` and trigger a warning, so filter on `success` before using the data. Two sampling modes: `"random"` (jitter `request.params` by `perturb_frac`) and `"grid"` (full factorial over `param_ranges`). Optional relative Gaussian noise on state values simulates experimental measurement noise.
- **`perturbation_study.py`** — the driver script. Picks a model from `MODELS` via `MODEL_NAME`, runs `generate_dataset` in `"random"` mode, overlays all per-run trajectories per state, and writes `perturbation_results.csv`/`.png`.

Registered models:
- `bioprocess` — 5-state cell-culture kinetics (`VCD`, `Glc`, `Gln`, `Lac`, `Amm`), rates `mu, qGlc, qGln, qLac, qAmm`. No substrate-limiting feedback, so `VCD` grows as `exp(mu*t)` unconditionally.
- `sird` — 4-state epidemic compartmental model (`S`, `I`, `R`, `D`), rates `beta, gamma, mu`. Total population `N` is computed as `S+I+R+D` (constant, since the compartments sum to zero net change) rather than being a free parameter.

`perturbation_results.csv`/`.png` are generated artifacts, not source of truth — they reflect whichever `MODEL_NAME` was active the last time `perturbation_study.py` was run, and go stale relative to the code whenever the model or params change without a rerun.

This project directory lives under iCloud Desktop sync — rapid successive writes to the same filename (e.g. scripted `sed`+`mv`) can produce iCloud conflict copies (`file 2.py`) instead of overwriting in place.
