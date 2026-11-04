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


MODEL = ODESystemSpec(
    name="sird",
    state_names=["S", "I", "R", "D"],
    param_names=["beta", "gamma", "mu"],
    rhs=rhs,
    conserved_total=True,
    default_params={
        "beta": 0.3,
        "gamma": 0.1,
        "mu": 0.01,
    },
    default_initial_conditions={"S0": 990.0, "I0": 10.0, "R0": 0.0, "D0": 0.0},
    default_t_span=(0, 160),
    default_n_points=161,
)
