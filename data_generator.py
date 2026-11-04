"""
Dataset Generator
=============================

Runs solve_ode() many times to build a synthetic dataset of trajectories
for a range of rate parameters, for any ODESystemSpec (see models/) -- the
"digital twin" data source for downstream modeling or a future
simulation chatbot.

Two sampling modes:
    "random" - request.params each scaled by +/- perturb_frac (uniform),
               matching a perturbation/sensitivity study.
    "grid"   - full factorial grid over param_ranges
               {name: (low, high, n_levels)}.

Optional Gaussian measurement noise (relative to each state value) can be
added to the trajectories to mimic real experimental data.

Failed integrations are kept and flagged rather than dropped: filter on the
"success" column before using the data.
"""

import dataclasses
import itertools
import warnings

import numpy as np
import pandas as pd

from ODE_solver import solve_ode, validate_request


def generate_dataset(
    spec,
    request,
    param_ranges=None,
    n_runs=200,
    mode="random",
    perturb_frac=0.1,
    noise_std=0.0,
    seed=None,
):
    """
    Simulate many parameter sets for `spec` and collect all trajectories.

    Parameters
    ----------
    spec : ODESystemSpec
    request : SimulationRequest
        Base case. Its params are the base values for mode="random"; its
        initial conditions, time grid, inputs and solver options apply to
        every run.
    param_ranges : dict mapping spec.param_names to (low, high, n_levels) (mode="grid")

    Returns
    -------
    pandas.DataFrame in long format with columns:
        run, t, success, <spec.state_names>, <spec.param_names>
    A run that fails before producing any point gets one row with NaN t and states.
    """
    problems = validate_request(spec, request)
    if problems:
        raise ValueError(f"Invalid request for '{spec.name}': {problems}")

    rng = np.random.default_rng(seed)
    param_names = list(spec.param_names)

    if mode == "random":
        base_array = np.array([request.params[name] for name in param_names], dtype=float)
        param_sets = base_array * rng.uniform(
            1 - perturb_frac, 1 + perturb_frac, size=(n_runs, len(base_array))
        )
    elif mode == "grid":
        axes = [np.linspace(low, high, n) for low, high, n in
                (param_ranges[name] for name in param_names)]
        param_sets = np.array(list(itertools.product(*axes)))
    else:
        raise ValueError(f"unknown mode: {mode!r}")

    rows = []
    n_failed = 0
    for run, param_values in enumerate(param_sets):
        params = dict(zip(param_names, param_values))
        result = solve_ode(spec, dataclasses.replace(request, params=params))
        if not result.success:
            n_failed += 1
        if len(result.t) == 0:
            nan_states = [np.nan] * len(spec.state_names)
            rows.append([run, np.nan, result.success, *nan_states, *param_values])
        for i, t in enumerate(result.t):
            state_vals = [result.trajectories[key][i] for key in spec.state_names]
            if noise_std:
                state_vals = [
                    v + rng.normal(0, noise_std * abs(v)) for v in state_vals
                ]
            rows.append([run, t, result.success, *state_vals, *param_values])

    if n_failed:
        warnings.warn(
            f"{n_failed} of {len(param_sets)} runs failed to integrate; "
            f"see the 'success' column."
        )

    columns = ["run", "t", "success"] + list(spec.state_names) + param_names
    return pd.DataFrame(rows, columns=columns)
