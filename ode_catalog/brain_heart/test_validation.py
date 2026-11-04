"""
Checks the model against the reference's equations. These confirm the
equations were entered correctly; they aren't a comparison of simulation
results, so they don't raise the model's status.
"""

import itertools
import math

import numpy as np
import pytest

from ode_solver import solve_ode

from .model import ENTRY, rhs

PARAMS = ENTRY.spec.default_params


@pytest.mark.parametrize(
    "H, P, Psp",
    list(itertools.product([0.8, 100 / 60, 3.0], [50.0, 90.0, 140.0], [70.0, 90.0, 120.0])),
)
def test_rhs_matches_equations_3_and_4(H, P, Psp):
    Tp = 1 / (1 + math.exp(-PARAMS["asp"] * (P - Psp)))
    Ts = 1 - Tp
    expected_dH = (
        PARAMS["betaH"] * Ts / (1 + PARAMS["gamma"] * Tp)
        - PARAMS["VH"] * Tp
        + PARAMS["deltaH"] * (PARAMS["H0"] - H)
    )
    expected_dP = (
        -P / (PARAMS["Ra"] * (1 + PARAMS["alpha"] * Ts) * PARAMS["Ca"])
        + H * PARAMS["dV"] / PARAMS["Ca"]
    )

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


@pytest.mark.skip(reason="No forward-simulation result to compare against yet; see README.md")
def test_reproduces_published_result():
    """Replace this once reference data is in reference_data/."""
