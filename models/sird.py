"""
SIRD Epidemic Model
======================

4-state compartmental epidemic model:
    S - Susceptible
    I - Infected
    R - Recovered
    D - Deceased

Rates:
    beta  - transmission rate
    gamma - recovery rate
    mu    - mortality rate (of infected individuals)

N (total population) is not a free parameter -- it's taken as
S + I + R + D, which is constant since the compartments sum to zero
change (dS+dI+dR+dD = 0). If you need N to vary independently of the
initial condition, promote it to a parameter instead.
"""

from .base import ODESystemSpec


def rhs(t, X, beta, gamma, mu):
    """RHS of the SIRD compartmental ODE system."""
    S, I, R, D = X
    N = S + I + R + D

    dSdt = -beta * S * I / N
    dIdt = beta * S * I / N - gamma * I - mu * I
    dRdt = gamma * I
    dDdt = mu * I

    return [dSdt, dIdt, dRdt, dDdt]


def invariants(y0, params):
    """Initial conditions and rates must be non-negative, and I must be able to leave."""
    problems = []
    for key in ("S0", "I0", "R0", "D0"):
        if y0.get(key, 0) < 0:
            problems.append(f"Initial condition {key} must be non-negative.")
    for key in ("beta", "gamma", "mu"):
        if params.get(key, 0) < 0:
            problems.append(f"Parameter {key} must be non-negative.")
    if (params.get("gamma", 0) + params.get("mu", 0)) == 0:
        problems.append("gamma + mu is zero: infected population never leaves the I compartment.")
    return problems


MODEL = ODESystemSpec(
    name="sird",
    state_names=["S", "I", "R", "D"],
    param_names=["beta", "gamma", "mu"],
    rhs=rhs,
    conserved_total=True,
    invariant_check=invariants,
    default_params={
        "beta": 0.40,
        "gamma": 0.10,
        "mu": 0.01,
    },
    default_initial_conditions={"S0": 999990, "I0": 20, "R0": 0, "D0": 0},
    default_t_span=(0, 160),
    default_n_points=400,
)
