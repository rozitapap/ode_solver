"""Summary statistics for a SimulationResult."""

import numpy as np


def analyze_result(result, spec):
    """
    Compute summary statistics (peaks, finals) for every state.

    A "_conservation" entry is added only for systems whose spec declares
    conserved_total.
    """
    stats = {}
    for name, series in result.trajectories.items():
        peak_idx = int(np.argmax(series))
        trough_idx = int(np.argmin(series))
        stats[name] = {
            "initial": float(series[0]),
            "final": float(series[-1]),
            "peak_value": float(series[peak_idx]),
            "peak_time": float(result.t[peak_idx]),
            "min_value": float(series[trough_idx]),
            "min_time": float(result.t[trough_idx]),
        }

    if spec.conserved_total:
        total_initial = sum(float(series[0]) for series in result.trajectories.values())
        total_final = sum(float(series[-1]) for series in result.trajectories.values())
        stats["_conservation"] = {
            "total_initial": total_initial,
            "total_final": total_final,
            "abs_drift": abs(total_final - total_initial),
        }
    return stats
