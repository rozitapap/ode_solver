"""
Evidence for numerically verified. SIRD has no closed form in time, but it has
exact relations between its states, and an independent integrator to compare with.
"""

import numpy as np
import pytest

from ode_catalog.validation import assert_agrees, rk4_fixed_step
from ode_solver import solve_ode

from .model import ENTRY

TOLERANCE = 1e-5  # maximum error, relative to the state's largest value
SCENARIOS = [scenario.name for scenario in ENTRY.scenarios]


@pytest.mark.parametrize("scenario", SCENARIOS)
def test_susceptibles_match_the_exact_relation(scenario):
    request = ENTRY.request(scenario)
    result = solve_ode(ENTRY.spec, request)
    p = request.params
    y0 = request.initial_conditions
    N = sum(y0.values())
    removed = result.trajectories["R"] + result.trajectories["D"] - y0["R0"] - y0["D0"]

    expected = y0["S0"] * np.exp(-p["beta"] / N / (p["gamma"] + p["mu"]) * removed)

    assert_agrees(result.trajectories["S"], expected, TOLERANCE, label="S")


@pytest.mark.parametrize("scenario", SCENARIOS)
def test_deaths_and_recoveries_keep_their_exact_ratio(scenario):
    request = ENTRY.request(scenario)
    result = solve_ode(ENTRY.spec, request)
    p = request.params
    expected = result.trajectories["R"] * p["mu"] / p["gamma"]
    assert_agrees(result.trajectories["D"], expected, TOLERANCE, label="D")


@pytest.mark.parametrize("scenario", SCENARIOS)
def test_total_population_is_conserved(scenario):
    request = ENTRY.request(scenario)
    result = solve_ode(ENTRY.spec, request)
    total = sum(result.trajectories.values())
    expected = np.full_like(total, sum(request.initial_conditions.values()), dtype=float)
    assert_agrees(total, expected, TOLERANCE, label="total")


@pytest.mark.parametrize("scenario", SCENARIOS)
def test_agrees_with_an_independent_integrator(scenario):
    request = ENTRY.request(scenario)
    result = solve_ode(ENTRY.spec, request)
    p = request.params

    def separately_typed_rhs(t, y):
        S, I, R, D = y
        new_infections = p["beta"] * S * I / (S + I + R + D)
        return [
            -new_infections,
            new_infections - (p["gamma"] + p["mu"]) * I,
            p["gamma"] * I,
            p["mu"] * I,
        ]

    y0 = [request.initial_conditions[name + "0"] for name in ENTRY.spec.state_names]
    other = rk4_fixed_step(separately_typed_rhs, y0, result.t, substeps=20)

    for i, name in enumerate(ENTRY.spec.state_names):
        assert_agrees(result.trajectories[name], other[:, i], TOLERANCE, label=name)
