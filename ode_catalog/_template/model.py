"""
Template catalog entry: a two-species conversion A -> B at rate k.

Copy this folder and replace the contents. See ode_catalog/ADDING_A_MODEL.md.
"""

from ode_catalog.schema import ModelEntry, ModelStatus, Parameter, Scenario
from ode_solver.spec import ODESystemSpec


def rhs(t, y, k):
    """Derivatives in the same order as state_names, with parameters named as in PARAMETERS."""
    A, B = y
    return [-k * A, k * A]


PARAMETERS = [
    Parameter("k", 0.5, unit="s^-1", description="conversion rate", source=None),
]

BASE_CASE = Scenario(
    name="base_case",
    description="All of the material starts as A.",
    initial_conditions={"A0": 1.0, "B0": 0.0},
    t_span=(0, 10),
    n_points=101,
    source=None,
)

ENTRY = ModelEntry(
    spec=ODESystemSpec(
        name="template",
        state_names=["A", "B"],
        param_names=[parameter.name for parameter in PARAMETERS],
        rhs=rhs,
        description="Two-species conversion A -> B. dA/dt=-k*A, dB/dt=k*A.",
        conserved_total=True,
        default_params={parameter.name: parameter.value for parameter in PARAMETERS},
        default_initial_conditions=BASE_CASE.initial_conditions,
        default_t_span=BASE_CASE.t_span,
        default_n_points=BASE_CASE.n_points,
    ),
    title="Template model",
    parameters=PARAMETERS,
    scenarios=[BASE_CASE],
    default_scenario="base_case",
    time_unit="s",
    state_units={"A": "mol", "B": "mol"},
    status=ModelStatus.EXECUTABLE,
    reference=None,  # None means the model isn't from a publication
)
