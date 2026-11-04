"""
Checks the model against the reference's equations, and its results against a
separate implementation of those equations (numerically verified). None of this
compares against a published simulation result, so it doesn't support
scientifically reproduced.
"""

import itertools
import math

import numpy as np
import pytest

from ode_catalog.validation import assert_agrees, rk4_fixed_step
from ode_solver import solve_ode

from .model import ENTRY, rhs

PARAMS = ENTRY.spec.default_params
TOLERANCE = 1e-5  # maximum error, relative to the state's largest value


def equations_3_and_4(H, P, Psp):
    """Eqs. (3) and (4), typed from the paper separately from model.py."""
    Tp = 1 / (1 + math.exp(-PARAMS["asp"] * (P - Psp)))
    Ts = 1 - Tp
    dH = (
        PARAMS["betaH"] * Ts / (1 + PARAMS["gamma"] * Tp)
        - PARAMS["VH"] * Tp
        + PARAMS["deltaH"] * (PARAMS["H0"] - H)
    )
    dP = (
        -P / (PARAMS["Ra"] * (1 + PARAMS["alpha"] * Ts) * PARAMS["Ca"])
        + H * PARAMS["dV"] / PARAMS["Ca"]
    )
    return [dH, dP]


@pytest.mark.parametrize(
    "H, P, Psp",
    list(itertools.product([0.8, 100 / 60, 3.0], [50.0, 90.0, 140.0], [70.0, 90.0, 120.0])),
)
def test_rhs_matches_equations_3_and_4(H, P, Psp):
    expected_dH, expected_dP = equations_3_and_4(H, P, Psp)
    dH, dP = rhs(0.0, [H, P], **PARAMS, Psp=lambda t: Psp)
    assert dH == pytest.approx(expected_dH, rel=1e-12)
    assert dP == pytest.approx(expected_dP, rel=1e-12)


@pytest.mark.parametrize("scenario", ["resting", "setpoint_step"])
def test_trajectory_satisfies_equation_13(scenario):
    request = ENTRY.request(scenario)
    result = solve_ode(ENTRY.spec, request)
    psp = request.inputs["Psp"]

    for t, H, P in zip(result.t, result.trajectories["H"], result.trajectories["P"]):
        dH, _ = rhs(t, [H, P], **PARAMS, Psp=psp)
        s = 1 / (1 + np.exp(-PARAMS["asp"] * (P - psp(t))))
        Y = dH - PARAMS["deltaH"] * (PARAMS["H0"] - H)
        Q = (
            PARAMS["VH"] * PARAMS["gamma"] * s**2
            + s * (PARAMS["betaH"] + PARAMS["VH"] + PARAMS["gamma"] * Y)
            + (Y - PARAMS["betaH"])
        )
        assert Q == pytest.approx(0.0, abs=1e-12)


def test_resting_agrees_with_an_independent_integrator():
    result = solve_ode(ENTRY.spec, ENTRY.request("resting"))
    other = rk4_fixed_step(
        lambda t, y: equations_3_and_4(y[0], y[1], 90.0), [100 / 60, 90.0], result.t, substeps=20
    )
    assert_agrees(result.trajectories["H"], other[:, 0], TOLERANCE, label="H")
    assert_agrees(result.trajectories["P"], other[:, 1], TOLERANCE, label="P")


def test_setpoint_step_agrees_with_an_independent_integrator():
    result = solve_ode(ENTRY.spec, ENTRY.request("setpoint_step"))

    # Integrate each side of the jump in Psp separately, so no step straddles it.
    before_t = np.linspace(0, 30, 6001)
    after_t = np.linspace(30, 60, 6001)
    before = rk4_fixed_step(
        lambda t, y: equations_3_and_4(y[0], y[1], 90.0), [100 / 60, 90.0], before_t
    )
    after = rk4_fixed_step(lambda t, y: equations_3_and_4(y[0], y[1], 110.0), before[-1], after_t)
    t = np.r_[before_t, after_t[1:]]
    other = np.r_[before, after[1:]]

    assert_agrees(result.trajectories["H"], np.interp(result.t, t, other[:, 0]), TOLERANCE, label="H")
    assert_agrees(result.trajectories["P"], np.interp(result.t, t, other[:, 1]), TOLERANCE, label="P")


@pytest.mark.skip(reason="No forward-simulation result to compare against yet; see README.md")
def test_reproduces_published_result():
    """Replace this once reference data is in reference_data/."""
