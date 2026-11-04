"""Checks specific to this model. Replace these with your model's own evidence."""

import numpy as np

from ode_solver import solve_ode

from .model import ENTRY


def test_matches_the_exact_solution():
    result = solve_ode(ENTRY.spec, ENTRY.request())
    k = ENTRY.spec.default_params["k"]
    exact = np.exp(-k * result.t)
    np.testing.assert_allclose(result.trajectories["A"], exact, rtol=1e-6)
