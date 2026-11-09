"""
Evidence for numerically verified: in every damping regime the model has an
exact closed-form solution. Also checks the paper's printed solutions and
reported numbers against the model. None of this compares against the paper's
figures, so it doesn't support scientifically reproduced.
"""

import math

import numpy as np
import pytest

from ode_catalog.validation import assert_agrees
from ode_solver import solve_ode

from .model import ENTRY

TOLERANCE = 1e-5  # maximum error, relative to the state's largest value
SCENARIOS = [scenario.name for scenario in ENTRY.scenarios]


def exact_solution(t, m, c, k, y0, v0):
    """Displacement and velocity of m y'' + c y' + k y = 0, from the characteristic roots."""
    alpha = c / (2 * m)
    discriminant = c * c - 4 * m * k
    if discriminant < 0:
        w = math.sqrt(-discriminant) / (2 * m)
        A, B = y0, (v0 + alpha * y0) / w
        e = np.exp(-alpha * t)
        y = e * (A * np.cos(w * t) + B * np.sin(w * t))
        v = e * ((B * w - alpha * A) * np.cos(w * t) - (A * w + alpha * B) * np.sin(w * t))
    elif discriminant == 0:
        C1, C2 = y0, v0 + alpha * y0
        e = np.exp(-alpha * t)
        y = (C1 + C2 * t) * e
        v = (C2 - alpha * (C1 + C2 * t)) * e
    else:
        beta = math.sqrt(discriminant) / (2 * m)
        l1, l2 = -alpha + beta, -alpha - beta
        C1 = (v0 - l2 * y0) / (l1 - l2)
        C2 = y0 - C1
        y = C1 * np.exp(l1 * t) + C2 * np.exp(l2 * t)
        v = C1 * l1 * np.exp(l1 * t) + C2 * l2 * np.exp(l2 * t)
    return y, v


def simulate(scenario):
    request = ENTRY.request(scenario)
    result = solve_ode(ENTRY.spec, request)
    assert result.success, result.message
    return request, result


@pytest.mark.parametrize("scenario", SCENARIOS)
def test_matches_the_exact_solution(scenario):
    request, result = simulate(scenario)
    p, ic = request.params, request.initial_conditions
    y, v = exact_solution(result.t, p["m"], p["c"], p["k"], ic["y0"], ic["v0"])
    assert_agrees(result.trajectories["y"], y, TOLERANCE, label="y")
    assert_agrees(result.trajectories["v"], v, TOLERANCE, label="v")


@pytest.mark.parametrize(
    "scenario, regime",
    [("no_damping", "none"), ("underdamping", "under"), ("critical_damping", "critical"),
     ("overdamping", "over")],
)
def test_scenario_is_in_its_damping_regime(scenario, regime):
    p = ENTRY.request(scenario).params
    d = p["c"] ** 2 - 4 * p["m"] * p["k"]
    actual = "none" if p["c"] == 0 else "under" if d < 0 else "critical" if d == 0 else "over"
    assert actual == regime


# The paper's printed solutions, with its x read as t.
PRINTED = {
    "no_damping": lambda t: np.cos(4 * t),  # Eq. (4.2)
    "overdamping": lambda t: 4 / 3 * np.exp(-2 * t) - 1 / 3 * np.exp(-8 * t),  # Eq. (4.11)
}
CORRECTED = {
    # Eq. (4.5) prints e^{-t} cos(sqrt(15) t), i.e. B = 0; y'(0) = 0 needs B = 1/sqrt(15).
    "underdamping": lambda t: np.exp(-t) * (np.cos(math.sqrt(15) * t) + np.sin(math.sqrt(15) * t) / math.sqrt(15)),
    # Eq. (4.8) prints (1+4) e^{-4x}; C1 = 1, C2 = 4 in Eq. (2.2.7) give (1 + 4t) e^{-4t}.
    "critical_damping": lambda t: (1 + 4 * t) * np.exp(-4 * t),
}


@pytest.mark.parametrize("scenario", list(PRINTED) + list(CORRECTED))
def test_agrees_with_the_papers_solutions(scenario):
    _, result = simulate(scenario)
    expected = {**PRINTED, **CORRECTED}[scenario](result.t)
    assert_agrees(result.trajectories["y"], expected, TOLERANCE, label=scenario)


def test_printed_underdamped_solution_does_not_satisfy_the_initial_velocity():
    """Records why Eq. (4.5) is corrected: as printed, y'(0) = -1, not 0."""
    printed = lambda t: np.exp(-t) * np.cos(math.sqrt(15) * t)  # noqa: E731
    h = 1e-6
    assert (printed(h) - printed(-h)) / (2 * h) == pytest.approx(-1.0, abs=1e-6)


def test_reported_period_and_frequency_of_the_undamped_case():
    p = ENTRY.request("no_damping").params
    w0 = math.sqrt(p["k"] / p["m"])
    assert w0 == 4
    assert round(2 * math.pi / w0, 2) == 1.57  # T = 1.57
    assert round(w0 / (2 * math.pi), 2) == 0.64  # f = 0.64 Hz


@pytest.mark.parametrize("scenario", SCENARIOS)
def test_energy_is_conserved_only_without_damping(scenario):
    request, result = simulate(scenario)
    p = request.params
    y, v = result.trajectories["y"], result.trajectories["v"]
    energy = 0.5 * p["m"] * v**2 + 0.5 * p["k"] * y**2
    if p["c"] == 0:
        assert_agrees(energy, np.full_like(energy, energy[0]), TOLERANCE, label="energy")
    else:
        assert np.all(np.diff(energy) <= 1e-12 * energy[0])
        assert energy[-1] < 1e-6 * energy[0]
