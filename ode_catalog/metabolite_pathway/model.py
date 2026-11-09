"""
Metabolite pathway: a linear 5-species mass-action chain M1 -> M2 -> M3 -> M4 -> M5,
each step first-order with rate k_i.

See README.md in this folder for context.
"""

from ode_catalog.schema import (
    ConservedQuantity,
    ModelEntry,
    ModelStatus,
    Parameter,
    Scale,
    Scenario,
    StateVariable,
)


def rhs(t, y, k1, k2, k3, k4):
    """RHS of the linear first-order chain."""
    M1, M2, M3, M4, M5 = y

    dM1 = -k1 * M1
    dM2 = k1 * M1 - k2 * M2
    dM3 = k2 * M2 - k3 * M3
    dM4 = k3 * M3 - k4 * M4
    dM5 = k4 * M4

    return [dM1, dM2, dM3, dM4, dM5]


AMOUNT_RANGE = (0, 1000)
AMOUNT_BASIS = "heuristic: 0 to 10x the default total amount"
RATE_BASIS = "heuristic: 0.1x to 10x the default"

EQUATIONS = [
    r"-k_1 M_1",
    r"k_1 M_1 - k_2 M_2",
    r"k_2 M_2 - k_3 M_3",
    r"k_3 M_3 - k_4 M_4",
    r"k_4 M_4",
]

# TODO(units): no units are recorded for the amounts or the rate constants.
STATES = [
    StateVariable(
        key=f"M{i}", symbol=f"M_{i}", description=f"amount of metabolite {i}", unit=None,
        initial=100 if i == 1 else 0, equation=EQUATIONS[i - 1],
        ui_range=AMOUNT_RANGE, ui_range_basis=AMOUNT_BASIS, lower=0,
    )
    for i in range(1, 6)
]

PARAMETERS = [
    Parameter(
        key=f"k{i}", default=default, symbol=f"k_{i}",
        description=f"rate constant of M{i} -> M{i + 1}", unit=None,
        ui_range=(default / 10, default * 10), ui_range_basis=RATE_BASIS, scale=Scale.LOG, lower=0,
    )
    for i, default in enumerate([0.50, 0.30, 0.20, 0.10], start=1)
]

FROM_M1 = Scenario(
    name="from_m1",
    description="All of the material starts as M1.",
    t_span=(0, 50),
    n_points=300,
)

ENTRY = ModelEntry(
    id="metabolite_pathway",
    name="Linear metabolite pathway",
    version="1.0.0",
    domain="biochemistry",
    tags=("mass action", "reaction chain"),
    summary="Linear first-order chain M1 -> M2 -> M3 -> M4 -> M5.",
    rhs=rhs,
    time_unit=None,  # TODO(units): not recorded
    states=STATES,
    parameters=PARAMETERS,
    conserved_quantities=[ConservedQuantity("total material", ("M1", "M2", "M3", "M4", "M5"))],
    scenarios=[FROM_M1],
    default_scenario="from_m1",
    status=ModelStatus.NUMERICALLY_VERIFIED,
    evidence=(
        "Matches the exact Bateman solution within 1e-5 relative error."
    ),
    todo=(
        "Time unit isn't recorded.",
        "Amount units (M1-M5) aren't recorded.",
        "Rate-constant units (k1-k4) aren't recorded.",
        "Slider ranges (ui_range) are heuristic, not from a source.",
    ),
)
