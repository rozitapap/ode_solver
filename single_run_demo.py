"""Run one simulation of the bioprocess model, print its final state, and show a plot."""

from models import MODELS
from ode_solver import plot_results, solve_ode

spec = MODELS["bioprocess"]
result = solve_ode(spec, spec.default_request())

print(f"Model: {spec.name}")
print("Final state at t =", result.t[-1])
for key in spec.state_names:
    print(f"  {key}: {result.trajectories[key][-1]:.4f}")

plot_results(result, spec, show=True)
