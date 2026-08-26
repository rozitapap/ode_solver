# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this directory.

## Purpose

Each file here defines one ODE system as an `ODEModel` instance (see `base.py`). Nothing outside this directory should hardcode a state name, parameter name, or RHS equation — that knowledge belongs here, not in `ODE_solver.py` or `data_generator.py`.

## Adding a new ODE system

1. Create `models/<name>.py` following `bioprocess.py`/`sird.py` as templates:
   - an `rhs(X, t, *params)` function returning a list of derivatives, one per state, in the same order as `state_keys`
   - a module-level `MODEL = ODEModel(...)` instance with `name`, `state_keys`, `param_names` (order must match the positional args of `rhs` after `X, t`), `rhs`, `default_params` (dict keyed by `param_names`), `default_U_init`, `default_t_pred`
2. Register it in `__init__.py`: add `from . import <name>` and `"<name>": <name>.MODEL` to `MODELS`.

No changes to `ODE_solver.py` or `data_generator.py` are needed — both are driven entirely by the `ODEModel` passed in.

## Conventions

- `param_names` order and `rhs`'s positional parameter order after `(X, t, ...)` must match exactly — this is the one place a mismatch would silently produce wrong results, since `ODE_solver.solve_ode` flattens the `params` dict to a positional tuple via `param_names` right before calling `odeint`.
- Keep `default_params`/`default_U_init` physically sensible for the system (e.g. non-negative concentrations/populations) since they're used as the base case for perturbation studies.
