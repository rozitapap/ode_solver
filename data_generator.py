"""
Dataset Generator
=============================

Runs solve_ode() many times to build a synthetic dataset of trajectories
for a range of rate parameters, for any ODEModel (see models/) -- the
"digital twin" data source for downstream modeling or a future
simulation chatbot.

Two sampling modes:
    "random" - base_params each scaled by +/- perturb_frac (uniform),
               matching a perturbation/sensitivity study.
    "grid"   - full factorial grid over param_ranges
               {name: (low, high, n_levels)}.

Optional Gaussian measurement noise (relative to each state value) can be
added to the trajectories to mimic real experimental data.
"""

import itertools

import numpy as np
import pandas as pd

from ODE_solver import solve_ode


def generate_dataset(
    model,
    U_init,
    t_pred,
    base_params=None,
    param_ranges=None,
    n_runs=200,
    mode="random",
    perturb_frac=0.1,
    noise_std=0.0,
    seed=None,
):
    """
    Simulate many parameter sets for `model` and collect all trajectories.

    Parameters
    ----------
    model : ODEModel
    base_params : dict mapping model.param_names to a value (mode="random")
    param_ranges : dict mapping model.param_names to (low, high, n_levels) (mode="grid")

    Returns
    -------
    pandas.DataFrame in long format with columns:
        run, t, <model.state_keys>, <model.param_names>
    """
    rng = np.random.default_rng(seed)
    param_names = model.param_names

    if mode == "random":
        base_array = np.array([base_params[name] for name in param_names], dtype=float)
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
    for run, param_values in enumerate(param_sets):
        params = dict(zip(param_names, param_values))
        results = solve_ode(model, params, U_init, t_pred)
        for i, t in enumerate(t_pred):
            state_vals = [results[key][i] for key in model.state_keys]
            if noise_std:
                state_vals = [
                    v + rng.normal(0, noise_std * abs(v)) for v in state_vals
                ]
            rows.append([run, t, *state_vals, *param_values])

    columns = ["run", "t"] + model.state_keys + list(param_names)
    return pd.DataFrame(rows, columns=columns)
