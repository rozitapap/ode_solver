"""Evidence for numerically verified: with constant rates, the model has an exact solution."""

import numpy as np

from ode_catalog.validation import assert_agrees
from ode_solver import solve_ode

from .model import ENTRY

TOLERANCE = 1e-5  # maximum error, relative to the state's largest value


def test_matches_the_exact_solution():
    request = ENTRY.request()
    result = solve_ode(ENTRY.spec, request)
    p = request.params
    y0 = request.initial_conditions

    growth = np.exp(p["mu"] * result.t)
    produced_per_q = y0["VCD0"] * (growth - 1) / p["mu"]  # integral of VCD over time

    expected = {
        "VCD": y0["VCD0"] * growth,
        "Glc": y0["Glc0"] - p["qGlc"] * produced_per_q,
        "Gln": y0["Gln0"] - p["qGln"] * produced_per_q,
        "Lac": y0["Lac0"] + p["qLac"] * produced_per_q,
        "Amm": y0["Amm0"] + p["qAmm"] * produced_per_q,
    }
    for name, values in expected.items():
        assert_agrees(result.trajectories[name], values, TOLERANCE, label=name)
