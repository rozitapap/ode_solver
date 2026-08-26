# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

No build step, linter, or test suite is configured in this repo. No dependency manifest exists either (no `requirements.txt`/`pyproject.toml`) — runtime dependencies are `numpy`, `scipy`, `pandas`, `matplotlib`.

- Run a single-simulation demo (bioprocess model, prints final state, shows a plot): `python3 ODE_solver.py`
- Run the perturbation study (200 runs, writes `perturbation_results.csv` and `.png`): `python3 perturbation_study.py`
- Switch which ODE system the study runs against: edit `MODEL_NAME` at the top of `perturbation_study.py` (currently `"bioprocess"`; `"sird"` is also registered).

## Architecture

This codebase is ODE-system-agnostic: `ODE_solver.py` and `data_generator.py` contain no knowledge of any specific system. All system-specific detail (state names, parameter names, the RHS equations, default values) lives in `models/`.

- **`models/base.py`** — the `ODEModel` dataclass is the contract every system implements: `state_keys`, `param_names`, `rhs(X, t, *params)`, and defaults for params/initial conditions/time grid.
- **`models/__init__.py`** — `MODELS` dict registers each system by name (`"bioprocess"`, `"sird"`). This is the single extension point: adding a new ODE system means adding `models/<name>.py` with one `ODEModel` instance and one registry entry — no changes needed to the solver or generator.
- **`ODE_solver.py`** — `solve_ode(model, params, U_init, t_pred)` is the only integration entry point, wrapping `scipy.integrate.odeint`. `params` is a dict; it's flattened to a positional tuple via `model.param_names` right at the `odeint` call site — this ordering is deliberate, to avoid param-mismatch bugs across models with different param counts/orders. `plot_results(results, model)` draws one subplot per `model.state_keys`.
- **`data_generator.py`** — `generate_dataset(model, ...)` calls `solve_ode()` across many parameter sets and returns one long-format `DataFrame` (`run`, `t`, state columns, param columns). Two sampling modes: `"random"` (jitter `base_params` by `perturb_frac`) and `"grid"` (full factorial over `param_ranges`). Optional relative Gaussian noise on state values simulates experimental measurement noise.
- **`perturbation_study.py`** — the driver script. Picks a model from `MODELS` via `MODEL_NAME`, runs `generate_dataset` in `"random"` mode, overlays all per-run trajectories per state, and writes `perturbation_results.csv`/`.png`.

Registered models:
- `bioprocess` — 5-state cell-culture kinetics (`VCD`, `Glc`, `Gln`, `Lac`, `Amm`), rates `mu, qGlc, qGln, qLac, qAmm`. No substrate-limiting feedback, so `VCD` grows as `exp(mu*t)` unconditionally.
- `sird` — 4-state epidemic compartmental model (`S`, `I`, `R`, `D`), rates `beta, gamma, mu`. Total population `N` is computed as `S+I+R+D` (constant, since the compartments sum to zero net change) rather than being a free parameter.

`perturbation_results.csv`/`.png` are generated artifacts, not source of truth — they reflect whichever `MODEL_NAME` was active the last time `perturbation_study.py` was run, and go stale relative to the code whenever the model or params change without a rerun.

This project directory lives under iCloud Desktop sync — rapid successive writes to the same filename (e.g. scripted `sed`+`mv`) can produce iCloud conflict copies (`file 2.py`) instead of overwriting in place.
