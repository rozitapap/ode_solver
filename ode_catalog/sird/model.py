"""
SIRD epidemic model: S (Susceptible), I (Infected), R (Recovered), D (Deceased).

N (total population) is not a free parameter -- it's taken as
S + I + R + D, which is constant since the compartments sum to zero
change (dS+dI+dR+dD = 0). If you need N to vary independently of the
initial condition, promote it to a parameter instead.

See README.md in this folder for context.
"""

from ode_catalog.schema import ModelEntry, ModelStatus, Parameter, Scenario
from ode_solver.spec import ODESystemSpec


def rhs(t, y, beta, gamma, mu):
    """RHS of the SIRD compartmental ODE system."""
    S, I, R, D = y
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


PARAMETERS = [
    Parameter("beta", 0.40, unit="day^-1", description="transmission rate"),
    Parameter("gamma", 0.10, unit="day^-1", description="recovery rate"),
    Parameter("mu", 0.01, unit="day^-1", description="mortality rate of infected individuals"),
]

OUTBREAK = Scenario(
    name="outbreak",
    description="About 1,000,000 people with 20 initial infections.",
    initial_conditions={"S0": 999990, "I0": 20, "R0": 0, "D0": 0},
    t_span=(0, 160),
    n_points=400,
)

SMALL_POPULATION = Scenario(
    name="small_population",
    description="1,000 people with 10 initial infections and a lower transmission rate.",
    initial_conditions={"S0": 990.0, "I0": 10.0, "R0": 0.0, "D0": 0.0},
    t_span=(0, 160),
    n_points=161,
    param_overrides={"beta": 0.3},
)

ENTRY = ModelEntry(
    spec=ODESystemSpec(
        name="sird",
        state_names=["S", "I", "R", "D"],
        param_names=[parameter.name for parameter in PARAMETERS],
        rhs=rhs,
        description="Compartmental epidemic model with constant total population.",
        invariant_check=invariants,
        conserved_total=True,
        default_params={parameter.name: parameter.value for parameter in PARAMETERS},
        default_initial_conditions=OUTBREAK.initial_conditions,
        default_t_span=OUTBREAK.t_span,
        default_n_points=OUTBREAK.n_points,
    ),
    title="SIRD epidemic model",
    parameters=PARAMETERS,
    scenarios=[OUTBREAK, SMALL_POPULATION],
    default_scenario="outbreak",
    time_unit="day",
    state_units={"S": "people", "I": "people", "R": "people", "D": "people"},
    status=ModelStatus.NUMERICALLY_VERIFIED,
    evidence=(
        "Matches the exact S-versus-removed and D-versus-R relations, and a separate RK4 implementation, within 1e-5 relative error."
    ),
    reference=None,
)
