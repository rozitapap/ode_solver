"""
Restructuring must not change results: compare current output against the
numbers recorded before it started. See tests/baseline/README.md.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ode_catalog import CATALOG
from ode_solver import generate_dataset, solve_ode

BASELINE = Path(__file__).parent / "baseline"

# Baseline file -> (model, scenario)
CASES = {
    "bioprocess": ("bioprocess", "batch_culture"),
    "sird_population_1000": ("sird", "small_population"),
    "sird_population_1000000": ("sird", "outbreak"),
    "metabolite_pathway": ("metabolite_pathway", "from_m1"),
    "brain_heart_resting": ("brain_heart", "resting"),
    "brain_heart_setpoint_step": ("brain_heart", "setpoint_step"),
}


@pytest.mark.parametrize("name", CASES)
def test_trajectories_match_baseline(name):
    model_name, scenario = CASES[name]
    entry = CATALOG[model_name]
    result = solve_ode(entry.spec, entry.request(scenario))
    assert result.success, result.message

    expected = pd.read_csv(BASELINE / f"{name}.csv")
    np.testing.assert_allclose(result.t, expected["t"], rtol=1e-7, atol=1e-7)
    assert list(result.trajectories) == list(expected.columns[1:])
    for state, series in result.trajectories.items():
        np.testing.assert_allclose(series, expected[state], rtol=1e-7, atol=1e-7, err_msg=state)


def test_sweep_matches_baseline():
    entry = CATALOG["bioprocess"]
    actual = generate_dataset(entry.spec, entry.request(), n_runs=3, seed=0)
    expected = pd.read_csv(BASELINE / "sweep_bioprocess.csv")
    assert list(actual.columns) == list(expected.columns)
    assert actual["success"].all() and expected["success"].all()
    numeric = [column for column in expected.columns if column != "success"]
    np.testing.assert_allclose(actual[numeric], expected[numeric], rtol=1e-7, atol=1e-7)
