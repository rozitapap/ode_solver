"""
Generic ODE Solver
======================

Integrates any ODE system described by an ODESystemSpec (see models/base.py)
using scipy.integrate.solve_ivp. Contains no knowledge of any specific
system -- see models/ for the bioprocess and SIRD models.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

from models.base import SimulationResult


def validate_request(spec, request):
    """Check a SimulationRequest against spec; returns a list of problems (empty if OK)."""
    problems = []

    missing_params = [name for name in spec.param_names if name not in request.params]
    if missing_params:
        problems.append(f"Missing parameters: {sorted(missing_params)}")

    unknown_params = [name for name in request.params if name not in spec.param_names]
    if unknown_params:
        problems.append(f"Unknown parameters: {sorted(unknown_params)}")

    missing_ics = [name + "0" for name in spec.state_names
                   if name + "0" not in request.initial_conditions]
    if missing_ics:
        problems.append(f"Missing initial conditions: {sorted(missing_ics)}")

    missing_inputs = [name for name in spec.input_names if name not in request.inputs]
    if missing_inputs:
        problems.append(f"Missing time-varying inputs: {sorted(missing_inputs)}")

    if request.t_span[1] <= request.t_span[0]:
        problems.append(f"t_span must be increasing, got {request.t_span}")

    if request.n_points < 2:
        problems.append(f"n_points must be >= 2, got {request.n_points}")

    if not problems and spec.invariant_check is not None:
        problems.extend(spec.invariant_check(request.initial_conditions, request.params))

    return problems


def solve_ode(spec, request):
    """
    Integrate an ODE system over a request's time grid.

    Parameters
    ----------
    spec : ODESystemSpec
    request : SimulationRequest

    Returns
    -------
    SimulationResult. If integration fails, success is False and the
    trajectories hold only the points reached before the failure.
    """
    y0 = [request.initial_conditions[name + "0"] for name in spec.state_names]
    inputs = {name: request.inputs[name] for name in spec.input_names}

    # solve_ivp never returns if the derivative at the first point is inf/nan.
    f0 = spec.rhs(request.t_span[0], y0, **request.params, **inputs)
    if not (np.isfinite(y0).all() and np.isfinite(f0).all()):
        return SimulationResult(
            system_name=spec.name,
            t=np.array([]),
            trajectories={name: np.array([]) for name in spec.state_names},
            params=dict(request.params),
            initial_conditions=dict(request.initial_conditions),
            success=False,
            message="Non-finite initial state or derivative.",
            metadata={"n_func_evals": 1, "method": request.method},
        )

    sol = solve_ivp(
        lambda t, y: spec.rhs(t, y, **request.params, **inputs),
        request.t_span,
        y0,
        method=request.method,
        t_eval=request.t_eval(),
        rtol=request.rtol,
        atol=request.atol,
        max_step=request.max_step,
    )

    # solve_ivp can report success on a trajectory that has gone inf/nan.
    finite = bool(np.isfinite(sol.y).all())

    return SimulationResult(
        system_name=spec.name,
        t=sol.t,
        trajectories={name: sol.y[i] for i, name in enumerate(spec.state_names)},
        params=dict(request.params),
        initial_conditions=dict(request.initial_conditions),
        success=bool(sol.success) and finite,
        message=sol.message if finite else "Trajectory contains non-finite values.",
        metadata={"n_func_evals": sol.nfev, "method": request.method},
    )


def plot_results(result, spec, save_path=None):
    """Plot each state trajectory; show interactively or save to save_path."""
    n = len(spec.state_names)
    fig, axes = plt.subplots(n, 1, figsize=(6, 2.4 * n), sharex=True)
    for ax, key in zip(axes, spec.state_names):
        ax.plot(result.t, result.trajectories[key])
        ax.set_ylabel(key)
    axes[-1].set_xlabel("Time (h)")
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150)
    else:
        plt.show()
    return fig


if __name__ == "__main__":
    from models import MODELS

    spec = MODELS["bioprocess"]
    result = solve_ode(spec, spec.default_request())

    print(f"Model: {spec.name}")
    print("Final state at t =", result.t[-1])
    for key in spec.state_names:
        print(f"  {key}: {result.trajectories[key][-1]:.4f}")

    plot_results(result, spec)
