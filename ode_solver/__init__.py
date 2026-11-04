"""Model-agnostic ODE engine: define a system once, then solve, sweep, analyze and plot it."""

from .analysis import analyze_result
from .pipeline import Pipeline, PipelineRun
from .plotting import plot_results, plot_sweep
from .solve import solve_ode, validate_request
from .spec import ODESystemSpec, SimulationRequest, SimulationResult
from .sweep import generate_dataset

__all__ = [
    "ODESystemSpec",
    "Pipeline",
    "PipelineRun",
    "SimulationRequest",
    "SimulationResult",
    "analyze_result",
    "generate_dataset",
    "plot_results",
    "plot_sweep",
    "solve_ode",
    "validate_request",
]
