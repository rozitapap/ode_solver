"""
Mass-spring-damper: an unforced, damped linear oscillator m y'' + c y' + k y = 0,
written as the first-order system of Eq. (3.5) of the reference:

    y'  = v
    v'  = (-c v - k y) / m

y is the displacement from static equilibrium (positive downwards) and v = y'.
The damping coefficient c sets the regime: c = 0 undamped, c^2 < 4mk
underdamped, c^2 = 4mk critically damped, c^2 > 4mk overdamped.

See README.md in this folder for context and sources.
"""

from ode_catalog.schema import (
    ModelEntry,
    ModelStatus,
    Parameter,
    Reference,
    Scale,
    Scenario,
    StateVariable,
)


def rhs(t, y, m, c, k):
    """RHS of Eq. (3.5): y1' = y2, y2' = (-c*y2 - k*y1)/m."""
    displacement, velocity = y
    return [velocity, (-c * velocity - k * displacement) / m]


KRCHEVA_2024 = Reference(
    id="krcheva2024",
    citation=(
        'V. Krcheva, "Modelling a mass-spring system using a second-order homogeneous '
        'linear ordinary differential equation with constant coefficients," Balkan Journal '
        "of Applied Mathematics and Informatics (BJAMI), vol. 7, no. 1, pp. 7-18, 2024."
    ),
    authors=("Violeta Krcheva",),
    title=(
        "Modelling a mass-spring system using a second-order homogeneous linear ordinary "
        "differential equation with constant coefficients"
    ),
    year=2024,
    venue="Balkan Journal of Applied Mathematics and Informatics (BJAMI), vol. 7, no. 1, pp. 7-18",
    url="https://js.ugd.edu.mk/index.php/bjami/en/article/view/5932",
    note=(
        "ISSN 2545-4803. No article DOI is printed. The article is dated September 08, 2023 "
        "and appears in the 2024 issue."
    ),
)

STATES = [
    StateVariable(
        key="y", symbol="y", description="displacement from static equilibrium, positive downwards",
        unit=None,  # TODO(units): not printed; the SI parameter units imply m
        initial=1.0, equation="y'",
        ui_range=(-5, 5), ui_range_basis="heuristic: 5x the paper's initial displacement either way",
    ),
    StateVariable(
        key="v", symbol="y'", description="velocity",
        unit=None,  # TODO(units): not printed; the SI parameter units imply m/s
        initial=0.0, equation=r"\frac{-c y' - k y}{m}",
        ui_range=(-20, 20), ui_range_basis="heuristic",
    ),
]

PARAMETERS = [
    Parameter(
        key="m", default=1, symbol="m", description="mass", unit="kg",
        ui_range=(0.1, 10), ui_range_basis="heuristic: 0.1x to 10x the default", scale=Scale.LOG,
        lower=0, lower_exclusive=True, source="Section 3", reference="krcheva2024",
    ),
    Parameter(
        key="c", default=2, symbol="c", description="damping coefficient", unit="Ns/m",
        ui_range=(0, 20), ui_range_basis="heuristic: 0 to twice the paper's overdamped value",
        lower=0, source="Section 3 (underdamping case)", reference="krcheva2024",
        note=(
            "The paper uses c = 0, 2, 8 and 10 for its four cases; 2 (underdamping) is the "
            "default here. Section 2.2 requires c > 0 for a damped system; c = 0 is the "
            "undamped model of Eq. (2.1.1)."
        ),
    ),
    Parameter(
        key="k", default=16, symbol="k", description="spring constant", unit="N/m",
        ui_range=(1.6, 160), ui_range_basis="heuristic: 0.1x to 10x the default", scale=Scale.LOG,
        lower=0, lower_exclusive=True, source="Section 3", reference="krcheva2024",
    ),
]


def _case(name, c, description, figure):
    return Scenario(
        name=name,
        description=f"{description} m = 1, c = {c}, k = 16, y(0) = 1, y'(0) = 0.",
        t_span=(0, 20),
        n_points=401,
        param_overrides={} if c == 2 else {"c": c},
        source=f"Section 3 and {figure}",
    )


SCENARIOS = [
    _case("underdamping", 2, "Oscillation with decaying amplitude, Eq. (4.3).", "Fig. 5"),
    _case("no_damping", 0, "Harmonic oscillation at the natural frequency, Eq. (4.1).", "Fig. 4"),
    _case("critical_damping", 8, "Fastest return to equilibrium without oscillating, Eq. (4.6).",
          "Fig. 6"),
    _case("overdamping", 10, "Slow return to equilibrium without oscillating, Eq. (4.9).", "Fig. 7"),
]

ENTRY = ModelEntry(
    id="mass_spring_damper",
    name="Mass-spring-damper system",
    version="1.0.0",
    domain="mechanics",
    tags=("oscillator", "damping", "linear ODE"),
    summary="Unforced, damped linear oscillator m y'' + c y' + k y = 0, in four damping regimes.",
    rhs=rhs,
    time_unit="s",  # TODO(units): implied by the paper's frequencies in Hz, not printed for t
    states=STATES,
    parameters=PARAMETERS,
    scenarios=SCENARIOS,
    default_scenario="underdamping",
    references=[KRCHEVA_2024],
    equations_reference="krcheva2024",
    status=ModelStatus.NUMERICALLY_VERIFIED,
    evidence=(
        "Matches the exact closed-form solution in all four damping regimes within 1e-5 relative error."
    ),
    todo=(
        "Displacement and velocity units aren't printed; the SI parameter units imply m and m/s.",
        "Time unit 's' is implied by the paper's frequencies in Hz, not printed for t.",
        "The paper's printed solutions (4.5) and (4.8) don't match its own equations; see README.md.",
        "Slider ranges (ui_range) are heuristic, not from a source.",
    ),
)
