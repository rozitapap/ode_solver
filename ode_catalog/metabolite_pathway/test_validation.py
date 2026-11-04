"""Evidence for numerically verified: the chain has an exact solution."""

import math

import numpy as np

from ode_catalog.validation import assert_agrees
from ode_solver import solve_ode

from .model import ENTRY

TOLERANCE = 1e-5  # maximum error, relative to the state's largest value


def exact_solution(t, A, k):
    """Bateman equations for material A starting in M1, with distinct rates k."""
    species = []
    for n in range(1, len(k) + 1):
        total = np.zeros_like(t)
        for i in range(n):
            denominator = math.prod(k[j] - k[i] for j in range(n) if j != i)
            total += np.exp(-k[i] * t) / denominator
        species.append(A * math.prod(k[: n - 1]) * total)
    species.append(A - sum(species))
    return species


def test_matches_the_exact_solution():
    request = ENTRY.request()
    result = solve_ode(ENTRY.spec, request)
    k = [request.params[name] for name in ("k1", "k2", "k3", "k4")]

    expected = exact_solution(result.t, request.initial_conditions["M10"], k)

    for name, values in zip(ENTRY.spec.state_names, expected):
        assert_agrees(result.trajectories[name], values, TOLERANCE, label=name)
