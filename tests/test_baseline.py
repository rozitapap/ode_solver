"""
Restructuring must not change results: compare current output against the
numbers recorded before it started. See tests/baseline/README.md.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from models import MODELS
from models.brain_heart import psp_constant, psp_profile
from ode_solver import SimulationRequest, generate_dataset, solve_ode

BASELINE = Path(__file__).parent / "baseline"


def brain_heart_request(**kwargs):
    spec = MODELS["brain_heart"]
    return SimulationRequest(
        system_name=spec.name,
        params=dict(spec.default_params),
        initial_conditions={"H0": 100 / 60, "P0": 90},
        t_span=(0, 60),
        n_points=600,
        **kwargs,
    )


CASES = {
    "bioprocess": ("bioprocess", MODELS["bioprocess"].default_request()),
    "sird_population_1000": (
        "sird",
        SimulationRequest(
            system_name="sird",
            params={"beta": 0.3, "gamma": 0.1, "mu": 0.01},
            initial_conditions={"S0": 990.0, "I0": 10.0, "R0": 0.0, "D0": 0.0},
            t_span=(0, 160),
            n_points=161,
        ),
    ),
    "sird_population_1000000": ("sird", MODELS["sird"].default_request()),
    "metabolite_pathway": ("metabolite_pathway", MODELS["metabolite_pathway"].default_request()),
    "brain_heart_resting": ("brain_heart", brain_heart_request(inputs={"Psp": psp_constant})),
    "brain_heart_setpoint_step": (
        "brain_heart",
        brain_heart_request(max_step=0.1, inputs={"Psp": psp_profile}),
    ),
}


@pytest.mark.parametrize("name", CASES)
def test_trajectories_match_baseline(name):
    model_name, request = CASES[name]
    result = solve_ode(MODELS[model_name], request)
    assert result.success, result.message

    expected = pd.read_csv(BASELINE / f"{name}.csv")
    np.testing.assert_allclose(result.t, expected["t"], rtol=1e-7, atol=1e-7)
    assert list(result.trajectories) == list(expected.columns[1:])
    for state, series in result.trajectories.items():
        np.testing.assert_allclose(series, expected[state], rtol=1e-7, atol=1e-7, err_msg=state)


def test_sweep_matches_baseline():
    spec = MODELS["bioprocess"]
    actual = generate_dataset(spec, spec.default_request(), n_runs=3, seed=0)
    expected = pd.read_csv(BASELINE / "sweep_bioprocess.csv")
    assert list(actual.columns) == list(expected.columns)
    assert actual["success"].all() and expected["success"].all()
    numeric = [column for column in expected.columns if column != "success"]
    np.testing.assert_allclose(actual[numeric], expected[numeric], rtol=1e-7, atol=1e-7)
