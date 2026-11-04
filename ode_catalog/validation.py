"""
Tools for producing the evidence a ModelStatus requires. The engine doesn't
use these, and running them doesn't change any model's status.
"""

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd

REQUIRED_METADATA = (
    "paper",
    "figure",
    "scenario",
    "state",
    "x_unit",
    "y_unit",
    "x_scale",
    "y_scale",
    "tool",
    "digitized_by",
)


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


@dataclass(frozen=True)
class ReferenceCurve:
    """One published curve, digitized at least twice. x and y are in the model's units."""

    metadata: dict
    runs: Sequence[np.ndarray]  # each of shape (n_points, 2)

    @property
    def scenario(self):
        return self.metadata["scenario"]

    @property
    def state(self):
        return self.metadata["state"]

    @property
    def reviewed_by(self) -> Optional[str]:
        """Who confirmed the overlay shows the right figure and conditions, if anyone."""
        return self.metadata.get("overlay_reviewed_by")


@dataclass(frozen=True)
class CurveComparison:
    """How a simulated curve compares with a ReferenceCurve. Errors are relative."""

    curve_error: float
    peak_height_error: float
    peak_time_error: float
    final_value_error: float
    tolerance: float

    @property
    def passed(self):
        return max(
            self.curve_error,
            self.peak_height_error,
            self.peak_time_error,
            self.final_value_error,
        ) <= self.tolerance


def load_reference_curve(folder):
    """
    Load a digitized curve from a folder holding metadata.json and run*.csv files.

    Each CSV has columns x and y in the figure's units. x_scale and y_scale in
    the metadata convert them to the model's units.
    """
    folder = Path(folder)
    metadata = json.loads((folder / "metadata.json").read_text(encoding="utf-8"))
    missing = [key for key in REQUIRED_METADATA if key not in metadata]
    if missing:
        raise ValueError(f"{folder}: metadata is missing {missing}")

    runs = []
    for path in sorted(folder.glob("run*.csv")):
        data = pd.read_csv(path)
        x = data["x"].to_numpy(dtype=float) * metadata["x_scale"]
        y = data["y"].to_numpy(dtype=float) * metadata["y_scale"]
        order = np.argsort(x)
        runs.append(np.column_stack([x[order], y[order]]))
    if len(runs) < 2:
        raise ValueError(f"{folder}: needs at least two digitizing runs, found {len(runs)}")

    return ReferenceCurve(metadata=metadata, runs=runs)


def digitization_error(curve):
    """Largest RMS difference from the first digitizing run, relative to the curve's range."""
    first = curve.runs[0]
    scale = max(np.ptp(first[:, 1]), np.finfo(float).tiny)
    errors = []
    for run in curve.runs[1:]:
        difference = np.interp(first[:, 0], run[:, 0], run[:, 1]) - first[:, 1]
        errors.append(np.sqrt(np.mean(difference**2)) / scale)
    return float(max(errors))


def compare_to_curve(t, y, curve, safety_factor=2.0, minimum_tolerance=0.01):
    """
    Compare a simulated state with a ReferenceCurve.

    The tolerance is safety_factor times the measured digitization error, so a
    less precise figure gets a wider tolerance. minimum_tolerance stops two
    near-identical digitizing runs from setting an impossibly tight one.

    Values are compared relative to the curve's range, and the peak time
    relative to its duration. Passing is necessary for scientifically
    reproduced, but a person must also review the overlay from save_overlay.
    """
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    reference = curve.runs[0]
    x_ref, y_ref = reference[:, 0], reference[:, 1]
    value_scale = max(np.ptp(y_ref), np.finfo(float).tiny)
    time_scale = max(np.ptp(x_ref), np.finfo(float).tiny)

    y_at_ref = np.interp(x_ref, t, y)
    in_window = (t >= x_ref[0]) & (t <= x_ref[-1])
    t_window, y_window = t[in_window], y[in_window]

    return CurveComparison(
        curve_error=float(np.sqrt(np.mean((y_at_ref - y_ref) ** 2)) / value_scale),
        peak_height_error=float(abs(y_window.max() - y_ref.max()) / value_scale),
        peak_time_error=float(abs(t_window[np.argmax(y_window)] - x_ref[np.argmax(y_ref)]) / time_scale),
        final_value_error=float(abs(y_at_ref[-1] - y_ref[-1]) / value_scale),
        tolerance=max(safety_factor * digitization_error(curve), minimum_tolerance),
    )


def save_overlay(t, y, curve, path):
    """Save the simulation over the digitized points, for a person to review."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(t, y, label="simulation")
    for i, run in enumerate(curve.runs, start=1):
        ax.plot(run[:, 0], run[:, 1], "x", markersize=4, label=f"digitized run {i}")
    ax.set_xlabel("Time")
    ax.set_ylabel(curve.state)
    ax.set_title(f"{curve.metadata['figure']}: {curve.scenario}")
    ax.legend(frameon=False)
    fig.tight_layout()
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path
