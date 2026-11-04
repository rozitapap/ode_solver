"""Plots for a SimulationResult."""

import os


def plot_results(result, spec, save_path=None, show=False):
    """Plot each state trajectory on its own axis; optionally save and/or show it."""
    import matplotlib.pyplot as plt

    n = len(spec.state_names)
    fig, axes = plt.subplots(n, 1, figsize=(6, 2.4 * n), sharex=True)
    for ax, key in zip(axes, spec.state_names):
        ax.plot(result.t, result.trajectories[key])
        ax.set_ylabel(key)
    axes[-1].set_xlabel("Time")
    fig.tight_layout()
    if save_path:
        os.makedirs(os.path.dirname(save_path) or ".", exist_ok=True)
        fig.savefig(save_path, dpi=150)
    if show:
        plt.show()
    else:
        plt.close(fig)
    return save_path
