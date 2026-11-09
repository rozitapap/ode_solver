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

from ode_catalog.schema import ModelEntry, ModelStatus, Parameter, Scale, Scenario, StateVariable


def rhs(t, y, mu, qGlc, qGln, qLac, qAmm):
    """RHS of the batch mass-balance ODE system."""
    VCD, Glc, Gln, Lac, Amm = y

    dVCDdt = mu * VCD
    dGlcdt = -qGlc * VCD
    dGlndt = -qGln * VCD
    dLacdt = qLac * VCD
    dAmmdt = qAmm * VCD

    return [dVCDdt, dGlcdt, dGlndt, dLacdt, dAmmdt]


STATE_BASIS = "heuristic: 0 to 10x the default"
RATE_BASIS = "heuristic: 0.1x to 10x the default"


def _state(key, description, initial, equation):
    # TODO(units): no state units are recorded for this model.
    return StateVariable(
        key=key, symbol=rf"\mathrm{{{key}}}", description=description, unit=None, initial=initial,
        equation=equation, ui_range=(0, 10 * initial), ui_range_basis=STATE_BASIS, lower=0,
    )


STATES = [
    _state("VCD", "viable cell density", 0.5, r"\mu \mathrm{VCD}"),
    _state("Glc", "glucose concentration", 30.0, r"-q_{\mathrm{Glc}} \mathrm{VCD}"),
    _state("Gln", "glutamine concentration", 4.0, r"-q_{\mathrm{Gln}} \mathrm{VCD}"),
    _state("Lac", "lactate concentration", 0.1, r"q_{\mathrm{Lac}} \mathrm{VCD}"),
    _state("Amm", "ammonia concentration", 0.1, r"q_{\mathrm{Amm}} \mathrm{VCD}"),
]


def _rate(key, symbol, default, description, lower):
    # TODO(units): no parameter units are recorded for this model.
    return Parameter(
        key=key, default=default, symbol=symbol, description=description, unit=None,
        ui_range=(default / 10, default * 10), ui_range_basis=RATE_BASIS, scale=Scale.LOG,
        lower=lower,
    )


PARAMETERS = [
    # mu has no lower bound: a negative specific growth rate (net death) is meaningful.
    _rate("mu", r"\mu", 0.02, "specific growth rate", lower=None),
    _rate("qGlc", r"q_{\mathrm{Glc}}", 0.015, "specific glucose consumption rate", lower=0),
    _rate("qGln", r"q_{\mathrm{Gln}}", 0.002, "specific glutamine consumption rate", lower=0),
    _rate("qLac", r"q_{\mathrm{Lac}}", 0.01, "specific lactate production rate", lower=0),
    _rate("qAmm", r"q_{\mathrm{Amm}}", 0.001, "specific ammonia production rate", lower=0),
]

BATCH_CULTURE = Scenario(
    name="batch_culture",
    description="Batch culture over 200 time units.",
    t_span=(0, 200),
    n_points=201,
)

ENTRY = ModelEntry(
    id="bioprocess",
    name="Batch cell-culture kinetics",
    version="1.0.0",
    domain="bioprocess",
    tags=("cell culture", "mass balance"),
    summary="Batch cell-culture mass balance with constant specific rates.",
    rhs=rhs,
    time_unit="h",  # TODO(units): only the time label on the model's earlier plots supports hours
    states=STATES,
    parameters=PARAMETERS,
    scenarios=[BATCH_CULTURE],
    default_scenario="batch_culture",
    status=ModelStatus.NUMERICALLY_VERIFIED,
    evidence=(
        "Matches the exact closed-form solution within 1e-5 relative error."
    ),
    todo=(
        "State units (VCD, Glc, Gln, Lac, Amm) aren't recorded.",
        "Parameter units (mu, qGlc, qGln, qLac, qAmm) aren't recorded.",
        "Time unit 'h' comes only from the time label on earlier plots.",
        "No source is recorded for the parameter values.",
        "Slider ranges (ui_range) are heuristic, not from a source.",
    ),
)
