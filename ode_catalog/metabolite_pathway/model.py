"""
Metabolite pathway: a linear 5-species mass-action chain M1 -> M2 -> M3 -> M4 -> M5,
each step first-order with rate k_i.

See README.md in this folder for context.
"""

from ode_catalog.schema import ModelEntry, ModelStatus, Parameter, Scenario
from ode_solver.spec import ODESystemSpec


def rhs(t, y, k1, k2, k3, k4):
    """RHS of the linear first-order chain."""
    M1, M2, M3, M4, M5 = y

    dM1 = -k1 * M1
    dM2 = k1 * M1 - k2 * M2
    dM3 = k2 * M2 - k3 * M3
    dM4 = k3 * M3 - k4 * M4
    dM5 = k4 * M4

    return [dM1, dM2, dM3, dM4, dM5]


def invariants(y0, params):
    """Initial conditions and rate constants must be non-negative."""
    problems = []
    for key in ("M10", "M20", "M30", "M40", "M50"):
        if y0.get(key, 0) < 0:
            problems.append(f"Initial condition {key} must be non-negative.")
    for key in ("k1", "k2", "k3", "k4"):
        if params.get(key, 0) < 0:
            problems.append(f"Rate constant {key} must be non-negative.")
    return problems


PARAMETERS = [
    Parameter("k1", 0.50, description="rate constant of M1 -> M2"),
    Parameter("k2", 0.30, description="rate constant of M2 -> M3"),
    Parameter("k3", 0.20, description="rate constant of M3 -> M4"),
    Parameter("k4", 0.10, description="rate constant of M4 -> M5"),
]

FROM_M1 = Scenario(
    name="from_m1",
    description="All of the material starts as M1.",
    initial_conditions={"M10": 100, "M20": 0, "M30": 0, "M40": 0, "M50": 0},
    t_span=(0, 50),
    n_points=300,
)

ENTRY = ModelEntry(
    spec=ODESystemSpec(
        name="metabolite_pathway",
        state_names=["M1", "M2", "M3", "M4", "M5"],
        param_names=[parameter.name for parameter in PARAMETERS],
        rhs=rhs,
        description="Linear first-order chain M1 -> M2 -> M3 -> M4 -> M5.",
        invariant_check=invariants,
        conserved_total=True,
        default_params={parameter.name: parameter.value for parameter in PARAMETERS},
        default_initial_conditions=FROM_M1.initial_conditions,
        default_t_span=FROM_M1.t_span,
        default_n_points=FROM_M1.n_points,
    ),
    title="Linear metabolite pathway",
    parameters=PARAMETERS,
    scenarios=[FROM_M1],
    default_scenario="from_m1",
    time_unit=None,
    state_units={"M1": None, "M2": None, "M3": None, "M4": None, "M5": None},
    status=ModelStatus.NUMERICALLY_VERIFIED,
    evidence=(
        "Matches the exact Bateman solution within 1e-5 relative error."
    ),
    reference=None,
)
