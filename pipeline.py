"""
Fixed simulation pipeline:

    validate_request -> solve_ode -> analyze_result -> plot_results

Each step is logged, and the run stops at the first failure.
"""

from dataclasses import dataclass, field
from typing import Optional

from models.base import SimulationRequest, SimulationResult
from ODE_solver import analyze_result, plot_results, solve_ode, validate_request


@dataclass
class PipelineRun:
    """Everything produced by one pipeline run."""

    request: SimulationRequest
    result: Optional[SimulationResult] = None
    analysis: Optional[dict] = None
    plot_path: Optional[str] = None
    logs: list = field(default_factory=list)
    ok: bool = False


def _summarize(state_names, analysis):
    """Turn raw per-state stats into a short, human-readable summary."""
    lines = ["[pipeline] summary:"]
    for name in state_names:
        st = analysis[name]
        lines.append(
            f"  - {name}: {st['initial']:.2f} -> {st['final']:.2f}, "
            f"peak {st['peak_value']:.2f} at t={st['peak_time']:.1f}"
        )
    return "\n".join(lines)


class Pipeline:
    """Turns an ODESystemSpec and a SimulationRequest into a full report."""

    def __init__(self, verbose=True):
        self.verbose = verbose

    def _log(self, run, message):
        run.logs.append(message)
        if self.verbose:
            print(message)

    def run(self, spec, request, plot_path=None, show_plot=False):
        run = PipelineRun(request=request)

        problems = validate_request(spec, request)
        if problems:
            self._log(run, "[pipeline] validation failed:")
            for problem in problems:
                self._log(run, f"  - {problem}")
            return run
        self._log(run, f"[pipeline] validation passed for '{spec.name}' "
                       f"({len(spec.state_names)} states, {len(spec.param_names)} params)")

        result = solve_ode(spec, request)
        run.result = result
        for name in result.metadata["defaulted_inputs"]:
            self._log(run, f"[pipeline] input '{name}' not provided, using default")
        if not result.success:
            self._log(run, f"[pipeline] integration failed: {result.message}")
            return run
        self._log(run, f"[pipeline] integration succeeded "
                       f"({result.metadata['n_func_evals']} function evaluations)")

        analysis = analyze_result(result, spec)
        run.analysis = analysis
        if "_conservation" in analysis:
            drift = analysis["_conservation"]["abs_drift"]
            rel_drift = drift / max(analysis["_conservation"]["total_initial"], 1e-12)
            self._log(run, f"[pipeline] conservation drift: {drift:.3e} ({rel_drift:.2%})")
        self._log(run, _summarize(spec.state_names, analysis))

        run.plot_path = plot_results(result, spec, save_path=plot_path, show=show_plot)
        if run.plot_path:
            self._log(run, f"[pipeline] plot saved to {run.plot_path}")

        run.ok = True
        self._log(run, "[pipeline] run complete.")
        return run
