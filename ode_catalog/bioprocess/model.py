"""
Bioprocess model: 5-state cell-culture kinetics.

    VCD - Viable Cell Density
    Glc - Glucose concentration
    Gln - Glutamine concentration
    Lac - Lactate concentration
    Amm - Ammonia concentration

Rates (mu, qGlc, qGln, qLac, qAmm) are constant parameters rather than
being predicted by a model. No substrate-limiting feedback, so VCD grows
as exp(mu*t) unconditionally -- keep mu small and qGlc/qGln low enough
that consumption doesn't outrun the initial Glc/Gln pools over t_span.

See README.md in this folder for context.
"""

from ode_catalog.schema import ModelEntry, ModelStatus, Parameter, Scenario
from ode_solver.spec import ODESystemSpec


def rhs(t, y, mu, qGlc, qGln, qLac, qAmm):
    """RHS of the batch mass-balance ODE system."""
    VCD, Glc, Gln, Lac, Amm = y

    dVCDdt = mu * VCD
    dGlcdt = -qGlc * VCD
    dGlndt = -qGln * VCD
    dLacdt = qLac * VCD
    dAmmdt = qAmm * VCD

    return [dVCDdt, dGlcdt, dGlndt, dLacdt, dAmmdt]


PARAMETERS = [
    Parameter("mu", 0.02, description="specific growth rate"),
    Parameter("qGlc", 0.015, description="specific glucose consumption rate"),
    Parameter("qGln", 0.002, description="specific glutamine consumption rate"),
    Parameter("qLac", 0.01, description="specific lactate production rate"),
    Parameter("qAmm", 0.001, description="specific ammonia production rate"),
]

BATCH_CULTURE = Scenario(
    name="batch_culture",
    description="Batch culture over 200 time units.",
    initial_conditions={"VCD0": 0.5, "Glc0": 30.0, "Gln0": 4.0, "Lac0": 0.1, "Amm0": 0.1},
    t_span=(0, 200),
    n_points=201,
)

ENTRY = ModelEntry(
    spec=ODESystemSpec(
        name="bioprocess",
        state_names=["VCD", "Glc", "Gln", "Lac", "Amm"],
        param_names=[parameter.name for parameter in PARAMETERS],
        rhs=rhs,
        description="Batch cell-culture mass balance with constant specific rates.",
        default_params={parameter.name: parameter.value for parameter in PARAMETERS},
        default_initial_conditions=BATCH_CULTURE.initial_conditions,
        default_t_span=BATCH_CULTURE.t_span,
        default_n_points=BATCH_CULTURE.n_points,
    ),
    title="Batch cell-culture kinetics",
    parameters=PARAMETERS,
    scenarios=[BATCH_CULTURE],
    default_scenario="batch_culture",
    time_unit="h",
    state_units={"VCD": None, "Glc": None, "Gln": None, "Lac": None, "Amm": None},
    status=ModelStatus.NUMERICALLY_VERIFIED,
    evidence=(
        "Matches the exact closed-form solution within 1e-5 relative error."
    ),
    reference=None,
)
