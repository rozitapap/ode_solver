"""
Tools for producing the evidence a ModelStatus requires. The engine doesn't
use these, and running them doesn't change any model's status.
"""

import numpy as np


def assert_agrees(actual, expected, tolerance, label=""):
    """Fail unless the maximum error, relative to the largest expected value, is within tolerance."""
    actual = np.asarray(actual, dtype=float)
    expected = np.asarray(expected, dtype=float)
    scale = max(np.abs(expected).max(), np.finfo(float).tiny)
    error = np.abs(actual - expected).max() / scale
    assert error <= tolerance, f"{label}: relative error {error:.2e} exceeds {tolerance:.2e}"


def rk4_fixed_step(f, y0, t_grid, substeps=1):
    """
    Integrate dy/dt = f(t, y) with the classical fixed-step Runge-Kutta method.

    This deliberately shares no code with solve_ivp, so it can serve as an
    independent implementation when checking a solver's results.

    Returns an array of shape (len(t_grid), len(y0)).
    """
    t_grid = np.asarray(t_grid, dtype=float)
    y = np.asarray(y0, dtype=float)
    out = np.empty((len(t_grid), len(y)))
    out[0] = y

    for i in range(len(t_grid) - 1):
        h = (t_grid[i + 1] - t_grid[i]) / substeps
        t = t_grid[i]
        for _ in range(substeps):
            k1 = np.asarray(f(t, y), dtype=float)
            k2 = np.asarray(f(t + h / 2, y + h / 2 * k1), dtype=float)
            k3 = np.asarray(f(t + h / 2, y + h / 2 * k2), dtype=float)
            k4 = np.asarray(f(t + h, y + h * k3), dtype=float)
            y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            t += h
        out[i + 1] = y

    return out
