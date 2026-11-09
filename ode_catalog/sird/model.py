"""
SIRD epidemic model: S (Susceptible), I (Infected), R (Recovered), D (Deceased).

N (total population) is not a free parameter -- it's taken as
S + I + R + D, which is constant since the compartments sum to zero
change (dS+dI+dR+dD = 0). If you need N to vary independently of the
initial condition, promote it to a parameter instead.

See README.md in this folder for context.
"""

from ode_catalog.schema import (
    Auxiliary,
    ConservedQuantity,
    ModelEntry,
    ModelStatus,
    Parameter,
    Scale,
    Scenario,
    StateVariable,
)


def rhs(t, y, beta, gamma, mu):
    """RHS of the SIRD compartmental ODE system."""
    S, I, R, D = y
    N = S + I + R + D

    dSdt = -beta * S * I / N
    dIdt = beta * S * I / N - gamma * I - mu * I
    dRdt = gamma * I
    dDdt = mu * I

    return [dSdt, dIdt, dRdt, dDdt]


def extra_checks(y0, params):
    """Infected people must be able to leave the I compartment."""
    if params["gamma"] + params["mu"] == 0:
        return ["gamma + mu is zero: infected population never leaves the I compartment."]
    return []


POPULATION_RANGE = (0, 2_000_000)
POPULATION_BASIS = "heuristic: 0 to about twice the default population"
RATE_BASIS = "heuristic: 0.1x to 10x the default"


def _compartment(key, symbol, description, initial, equation):
    return StateVariable(
        key=key, symbol=symbol, description=description, unit="people", initial=initial,
        equation=equation, ui_range=POPULATION_RANGE, ui_range_basis=POPULATION_BASIS, lower=0,
    )


STATES = [
    _compartment("S", "S", "susceptible people", 999990, r"-\frac{\beta S I}{N}"),
    _compartment("I", "I", "infected people", 20, r"\frac{\beta S I}{N} - \gamma I - \mu I"),
    _compartment("R", "R", "recovered people", 0, r"\gamma I"),
    _compartment("D", "D", "deceased people", 0, r"\mu I"),
]


def _rate(key, symbol, default, description):
    return Parameter(
        key=key, default=default, symbol=symbol, description=description, unit="day^-1",
        ui_range=(default / 10, default * 10), ui_range_basis=RATE_BASIS, scale=Scale.LOG, lower=0,
    )


PARAMETERS = [
    _rate("beta", r"\beta", 0.40, "transmission rate"),
    _rate("gamma", r"\gamma", 0.10, "recovery rate"),
    _rate("mu", r"\mu", 0.01, "mortality rate of infected individuals"),
]

OUTBREAK = Scenario(
    name="outbreak",
    description="About 1,000,000 people with 20 initial infections.",
    t_span=(0, 160),
    n_points=400,
)

SMALL_POPULATION = Scenario(
    name="small_population",
    description="1,000 people with 10 initial infections and a lower transmission rate.",
    t_span=(0, 160),
    n_points=161,
    initial_overrides={"S": 990.0, "I": 10.0, "R": 0.0, "D": 0.0},
    param_overrides={"beta": 0.3},
)

ENTRY = ModelEntry(
    id="sird",
    name="SIRD epidemic model",
    version="1.0.0",
    domain="epidemiology",
    tags=("compartmental", "infectious disease"),
    summary="Compartmental epidemic model with constant total population.",
    rhs=rhs,
    time_unit="day",
    states=STATES,
    parameters=PARAMETERS,
    auxiliaries=[Auxiliary("N", "S + I + R + D", "total population")],
    conserved_quantities=[ConservedQuantity("total population", ("S", "I", "R", "D"))],
    scenarios=[OUTBREAK, SMALL_POPULATION],
    default_scenario="outbreak",
    status=ModelStatus.NUMERICALLY_VERIFIED,
    evidence=(
        "Matches the exact S-versus-removed and D-versus-R relations, and a separate RK4 implementation, within 1e-5 relative error."
    ),
    extra_checks=extra_checks,
    todo=("Slider ranges (ui_range) are heuristic, not from a source.",),
)
