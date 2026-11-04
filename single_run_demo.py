"""Run one simulation of the bioprocess model, print its final state, and show a plot."""

from ode_catalog import CATALOG
from ode_solver import plot_results, solve_ode

entry = CATALOG["bioprocess"]
spec = entry.spec
result = solve_ode(spec, entry.request())

print(f"Model: {spec.name}")
print("Final state at t =", result.t[-1])
for key in spec.state_names:
    print(f"  {key}: {result.trajectories[key][-1]:.4f}")

plot_results(result, spec, show=True)
