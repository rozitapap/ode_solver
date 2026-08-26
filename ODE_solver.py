"""
Generic ODE Solver
======================

Integrates any ODE system described by an ODEModel spec (see
models/base.py) using scipy.integrate.odeint. Contains no knowledge of
any specific system -- see models/ for the bioprocess and SIRD models.
"""

import matplotlib.pyplot as plt
from scipy.integrate import odeint


def solve_ode(model, params, U_init, t_pred):
    """
    Integrate an ODE model over t_pred.

    Parameters
    ----------
    model : ODEModel
    params : dict mapping each name in model.param_names to a value
    U_init : array-like, shape (len(model.state_keys),)
        Initial state.
    t_pred : 1D array
        Time points at which to evaluate the solution.

    Returns
    -------
    dict with keys "t" and each of model.state_keys
    """
    param_tuple = tuple(params[name] for name in model.param_names)
    sol = odeint(model.rhs, U_init, t_pred, args=param_tuple)

    results = {"t": t_pred}
    for i, key in enumerate(model.state_keys):
        results[key] = sol[:, i]
    return results


def plot_results(results, model, save_path=None):
    """Plot each state trajectory; show interactively or save to save_path."""
    n = len(model.state_keys)
    fig, axes = plt.subplots(n, 1, figsize=(6, 2.4 * n), sharex=True)
    for ax, key in zip(axes, model.state_keys):
        ax.plot(results["t"], results[key])
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

    model = MODELS["bioprocess"]
    results = solve_ode(model, model.default_params, model.default_U_init, model.default_t_pred)

    print(f"Model: {model.name}")
    print("Final state at t =", model.default_t_pred[-1])
    for key in model.state_keys:
        print(f"  {key}: {results[key][-1]:.4f}")

    plot_results(results, model)
