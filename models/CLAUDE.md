# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this directory.

## Purpose

Each file here defines one ODE system as an `ODESystemSpec` instance (see `base.py`). `base.py` is kept byte-identical to `ode_agent/schema.py` in the PINNSLearning repo, so change both together. Nothing outside this directory should hardcode a state name, parameter name, or RHS equation — that knowledge belongs here, not in `ODE_solver.py` or `data_generator.py`.

## Adding a new ODE system

1. Create `models/<name>.py` following `bioprocess.py`/`sird.py` as templates:
   - an `rhs(t, y, **params)` function, with parameters named as in `param_names`, returning a list of derivatives, one per state, in the same order as `state_names`
   - a module-level `MODEL = ODESystemSpec(...)` instance with `name`, `state_names`, `param_names`, `rhs`, `default_params` (dict keyed by `param_names`), `default_initial_conditions` (dict keyed by `<state>0`), `default_t_span`, `default_n_points`
   - `conserved_total=True` only if the sum of all states really is constant
2. Register it in `__init__.py`: add `from . import <name>` and `"<name>": <name>.MODEL` to `MODELS`. `tests/test_model_contract.py` then covers it automatically.

No changes to `ODE_solver.py` or `data_generator.py` are needed — both are driven entirely by the `ODESystemSpec` passed in.

## Conventions

- Parameters are passed to `rhs` by name, so `param_names` order doesn't matter, but the names must match `rhs`'s argument names exactly.
- Keep `default_params`/`default_initial_conditions` physically sensible for the system (e.g. non-negative concentrations/populations) since they're used as the base case for perturbation studies.
