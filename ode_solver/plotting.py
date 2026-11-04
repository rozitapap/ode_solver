"""Plots for a SimulationResult or a sweep."""

import os


def _label(name, unit):
    return f"{name} ({unit})" if unit else name


def _state_axes(spec, time_unit, state_units):
    """Create one labelled axis per state, sharing a time axis."""
    import matplotlib.pyplot as plt

    state_units = state_units or {}
    n = len(spec.state_names)
    fig, axes = plt.subplots(n, 1, figsize=(6, 2.4 * n), sharex=True, squeeze=False)
    axes = axes[:, 0]
    for ax, key in zip(axes, spec.state_names):
        ax.set_ylabel(_label(key, state_units.get(key)))
    axes[-1].set_xlabel(_label("Time", time_unit))
    return fig, axes


def _finish(fig, save_path, show):
    import matplotlib.pyplot as plt

    fig.tight_layout()
    if save_path:
        os.makedirs(os.path.dirname(save_path) or ".", exist_ok=True)
        fig.savefig(save_path, dpi=150)
    if show:
        plt.show()
    else:
        plt.close(fig)
    return save_path


def plot_results(result, spec, save_path=None, show=False, time_unit=None, state_units=None):
    """Plot each state trajectory on its own axis; optionally save and/or show it."""
    fig, axes = _state_axes(spec, time_unit, state_units)
    for ax, key in zip(axes, spec.state_names):
        ax.plot(result.t, result.trajectories[key])
    return _finish(fig, save_path, show)


def plot_sweep(df, spec, save_path=None, show=False, time_unit=None, state_units=None, title=None):
    """Overlay every successful run of a generate_dataset() result, one axis per state."""
    fig, axes = _state_axes(spec, time_unit, state_units)
    for _, run_df in df[df["success"]].groupby("run"):
        for ax, key in zip(axes, spec.state_names):
            ax.plot(run_df["t"], run_df[key], alpha=0.2, color="C0")
    if title:
        axes[0].set_title(title)
    return _finish(fig, save_path, show)
