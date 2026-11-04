"""
Brain-heart model: baroreflex control of heart rate and arterial pressure.

State H = heart rate, P = mean arterial blood pressure.
Input Psp(t) = baroreflex set-point; output y = H.
Time is in seconds, so H is in beats/s and H0 = 100/60.
Note the parameter H0 (intrinsic heart rate) is separate from the
initial-condition key "H0" (H at t=0); they live in different dicts.

See README.md in this folder for context and sources.
"""

import numpy as np

from ode_catalog.schema import ModelEntry, ModelStatus, Parameter, Reference, Scenario
from ode_solver.spec import ODESystemSpec


def rhs(t, y, Ra, Ca, dV, betaH, gamma, VH, deltaH, H0, alpha, asp, Psp):
    """RHS of the baroreflex ODE system, Eqs. (3) and (4) of the reference."""
    H, P = y  # H(t), P(t) at the current time t
    dH = (
        betaH * (1.0 - 1.0 / (1.0 + np.exp(-asp * (P - Psp(t)))))
        / (1.0 + gamma / (1.0 + np.exp(-asp * (P - Psp(t)))))
        - VH / (1.0 + np.exp(-asp * (P - Psp(t))))
        + deltaH * (H0 - H)
    )
    dP = (
        -P / (Ra * (1.0 + alpha * (1.0 - 1.0 / (1.0 + np.exp(-asp * (P - Psp(t)))))) * Ca)
        + H * dV / Ca
    )
    return [dH, dP]


def invariants(y0, params):
    """Initial conditions and parameters must be physically sensible."""
    problems = []
    for key in ("H0", "P0"):
        if y0.get(key, 0) < 0:
            problems.append(f"Initial condition {key} must be non-negative.")
    for key in ("Ra", "Ca", "asp", "H0"):
        if params.get(key, 1) <= 0:
            problems.append(f"Parameter {key} must be positive.")
    for key in ("dV", "betaH", "gamma", "VH", "deltaH", "alpha"):
        if params.get(key, 0) < 0:
            problems.append(f"Parameter {key} must be non-negative.")
    return problems


def psp_constant(t):
    """Resting baroreflex set-point [mmHg]: the default when no Psp is given."""
    return 90.0


def psp_profile(t):
    """Baroreflex set-point [mmHg]: 90 mmHg, stepping up to 110 mmHg at t=30 s."""
    return 90.0 if t < 30.0 else 110.0


PARAMETERS = [
    Parameter(
        "Ra", 0.6, unit="mmHg·g·ml^-1", description="minimum arterial resistance", source="Table I",
        note="Unit as printed. Eq. (4) is only dimensionally consistent with mmHg·s·ml^-1.",
    ),
    Parameter(
        "Ca", 1.55, unit="mmHg·g^-1", description="arterial compliance", source="Table I",
        note="Unit as printed. Eq. (4) is only dimensionally consistent with ml·mmHg^-1.",
    ),
    Parameter("dV", 50, unit="ml", description="stroke volume", source="Table I"),
    Parameter("betaH", 0.84, unit="s^-2", description="sympathetic control of HR", source="Table I"),
    Parameter("gamma", 0.2, description="parasympathetic damping of betaH", source="Table I"),
    Parameter("VH", 1.17, unit="s^-1", description="parasympathetic control of HR", source="Table I"),
    Parameter(
        "deltaH", 1.7, unit="s^-1", description="relaxation time", source="Table I",
        note="Table I calls this a time, but its unit and its role in Eq. (3) are those of a rate.",
    ),
    Parameter(
        "H0", 100 / 60, unit="s^-1", description="intrinsic HR", source="Table I",
        note="Table I gives 100 min^-1; converted because this model's time is in seconds.",
    ),
    Parameter("alpha", 1.3, description="sympathetic effect on Ra", source="Table I"),
    Parameter(
        "asp", 0.07, unit="mmHg^-1",
        description="slope of sigmoid determining baroreflex response", source="Table I",
    ),
]

# The reference treats Psp as unknown and estimates it from data, so it doesn't
# publish these inputs or initial conditions. Their source is therefore None.
INITIAL_CONDITIONS = {"H0": 100 / 60, "P0": 90}

RESTING = Scenario(
    name="resting",
    description="Constant 90 mmHg set-point. P settles from 90 to about 85.7 mmHg within 10 s.",
    initial_conditions=INITIAL_CONDITIONS,
    t_span=(0, 60),
    n_points=600,
    inputs={"Psp": psp_constant},
    source=None,
)

SETPOINT_STEP = Scenario(
    name="setpoint_step",
    description="Set-point follows psp_profile: 90 mmHg, stepping up to 110 mmHg at t=30 s.",
    initial_conditions=INITIAL_CONDITIONS,
    t_span=(0, 60),
    n_points=600,
    inputs={"Psp": psp_profile},
    max_step=0.1,  # keep the solver from stepping over the jump in Psp
    source=None,
)

ENTRY = ModelEntry(
    spec=ODESystemSpec(
        name="brain_heart",
        state_names=["H", "P"],
        param_names=[parameter.name for parameter in PARAMETERS],
        input_names=["Psp"],
        rhs=rhs,
        description=(
            "Delay-free baroreflex model of heart rate H and mean arterial pressure P "
            "driven by the set-point Psp(t)."
        ),
        invariant_check=invariants,
        default_params={parameter.name: parameter.value for parameter in PARAMETERS},
        default_initial_conditions=RESTING.initial_conditions,
        default_t_span=RESTING.t_span,
        default_n_points=RESTING.n_points,
        default_inputs=RESTING.inputs,
    ),
    title="Autonomic cardiac regulation (delay-free baroreflex model)",
    parameters=PARAMETERS,
    scenarios=[RESTING, SETPOINT_STEP],
    default_scenario="resting",
    time_unit="s",
    state_units={"H": "s^-1", "P": "mmHg"},
    status=ModelStatus.EXECUTABLE,
    reference=Reference(
        citation=(
            "S. N. Sadoun, G. A. D'Inverno, A. Boutin, F. Cottin, and T.-M. Laleg-Kirati, "
            '"Physics-Informed Neural Estimation of State and Unknown Input in Autonomic '
            'Cardiac Dynamics with Left-Invertibility Constraints," 2026 European Control '
            "Conference (ECC), Reykjavík, Iceland, July 7-10, 2026, pp. 1795-1800."
        ),
    ),
    related_references=[
        Reference(
            citation=(
                'J. T. Ottesen, "Modelling of the baroreflex-feedback mechanism with '
                'time-delay," Journal of Mathematical Biology, vol. 36, no. 1, 1997.'
            ),
            note="Origin of the model, with a delay. Cited as reference [7]; not yet checked directly.",
        ),
        Reference(
            citation=(
                'P. Ataee, G. A. Dumont, and W. T. Boyce, "Baroreflex modeling in the genesis '
                'of stress reactivity using sigmoidal characteristic," Proceedings of the 2010 '
                "American Control Conference, IEEE, 2010."
            ),
            note="Added the operating set-point. Cited as reference [8]; not yet checked directly.",
        ),
    ],
)
