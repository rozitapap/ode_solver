"""
Template catalog entry: a two-species conversion A -> B at rate k.

Copy this folder and replace every value marked REPLACE. The result must pass
validate_entry(), which runs when ENTRY is created. See ode_catalog/ADDING_A_MODEL.md.
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


# REPLACE: derivatives in the same order as STATES; arguments named as PARAMETERS
# (and INPUTS, if any), after t and y.
def rhs(t, y, k):
    A, B = y
    return [-k * A, k * A]


# REPLACE: one StateVariable per state, in rhs order. `equation` is the LaTeX
# right-hand side of d(key)/dt and may only use declared symbols and t.
STATES = [
    StateVariable(
        key="A", symbol="A", description="amount of A", unit="mol", initial=1.0,
        equation=r"-k A", ui_range=(0, 10), ui_range_basis="heuristic: 0 to 10x the default",
        lower=0,
    ),
    StateVariable(
        key="B", symbol="B", description="amount of B", unit="mol", initial=0.0,
        equation=r"k A", ui_range=(0, 10), ui_range_basis="heuristic: 0 to 10x A's default",
        lower=0,
    ),
]

# REPLACE: one Parameter per rhs argument. Give a source and reference if the
# value comes from a publication; put conversions in `note`.
PARAMETERS = [
    Parameter(
        key="k", default=0.5, symbol="k", description="conversion rate", unit="s^-1",
        ui_range=(0.05, 5), ui_range_basis="heuristic: 0.1x to 10x the default",
        scale=Scale.LOG, lower=0, source=None, reference=None,
    ),
]

# REPLACE: at least one scenario. Anything not overridden uses the defaults above.
BASE_CASE = Scenario(
    name="base_case",
    description="All of the material starts as A.",
    t_span=(0, 10),
    n_points=101,
    source=None,
)

ENTRY = ModelEntry(
    id="template",  # REPLACE: the folder name; never change it once published
    name="Template model",  # REPLACE
    version="1.0.0",
    domain="example",  # REPLACE: e.g. epidemiology, bioprocess, physiology, biochemistry
    tags=(),
    summary="Two-species conversion A -> B.",  # REPLACE: one line
    rhs=rhs,
    time_unit="s",  # REPLACE: None if it isn't recorded
    states=STATES,
    parameters=PARAMETERS,
    conserved_quantities=[ConservedQuantity("total amount", ("A", "B"))],  # REPLACE or remove
    scenarios=[BASE_CASE],
    default_scenario="base_case",
    references=(),  # REPLACE: Reference(...) entries, if the model is from a publication
    equations_reference=None,  # REPLACE: the id of the reference the equations come from
    status=ModelStatus.EXECUTABLE,  # start here; see "Choosing a status"
    todo=(),  # REPLACE: anything you couldn't verify
)
