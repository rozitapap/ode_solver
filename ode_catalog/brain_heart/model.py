"""
Brain-heart model: baroreflex control of heart rate and arterial pressure.

State H = heart rate, P = mean arterial blood pressure.
Input Psp(t) = baroreflex set-point; output y = H.
Time is in seconds, so H is in beats/s and H0 = 100/60.
Note the parameter H0 (intrinsic heart rate) is separate from the
engine's initial-condition key "H0" (H at t=0). In this entry the initial
value lives on the state "H", so the two never meet.

See README.md in this folder for context and sources.
"""

import numpy as np

from ode_catalog.schema import (
    Auxiliary,
    InputVariable,
    ModelEntry,
    ModelStatus,
    Parameter,
    Reference,
    Scale,
    Scenario,
    Signal,
    StateVariable,
)


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


SADOUN_2026 = Reference(
    id="sadoun2026",
    citation=(
        "S. N. Sadoun, G. A. D'Inverno, A. Boutin, F. Cottin, and T.-M. Laleg-Kirati, "
        '"Physics-Informed Neural Estimation of State and Unknown Input in Autonomic '
        'Cardiac Dynamics with Left-Invertibility Constraints," 2026 European Control '
        "Conference (ECC), Reykjavík, Iceland, July 7-10, 2026, pp. 1795-1800."
    ),
    authors=("S. N. Sadoun", "G. A. D'Inverno", "A. Boutin", "F. Cottin", "T.-M. Laleg-Kirati"),
    title=(
        "Physics-Informed Neural Estimation of State and Unknown Input in Autonomic "
        "Cardiac Dynamics with Left-Invertibility Constraints"
    ),
    year=2026,
    venue="2026 European Control Conference (ECC), Reykjavík, Iceland, pp. 1795-1800",
)

OTTESEN_1997 = Reference(
    id="ottesen1997",
    citation=(
        'J. T. Ottesen, "Modelling of the baroreflex-feedback mechanism with '
        'time-delay," Journal of Mathematical Biology, vol. 36, pp. 41-63, 1997.'
    ),
    authors=("J. T. Ottesen",),
    title="Modelling of the baroreflex-feedback mechanism with time-delay",
    year=1997,
    venue="Journal of Mathematical Biology, vol. 36, pp. 41-63",
    note=(
        "Cited as the model's origin, reference [7]. Checked directly: it's a related "
        "three-state delay model, not this model with the delay removed. See README.md."
    ),
)

ATAEE_2010 = Reference(
    id="ataee2010",
    citation=(
        'P. Ataee, G. A. Dumont, and W. T. Boyce, "Baroreflex modeling in the genesis '
        'of stress reactivity using sigmoidal characteristic," Proceedings of the 2010 '
        "American Control Conference, IEEE, 2010."
    ),
    authors=("P. Ataee", "G. A. Dumont", "W. T. Boyce"),
    title="Baroreflex modeling in the genesis of stress reactivity using sigmoidal characteristic",
    year=2010,
    venue="Proceedings of the 2010 American Control Conference, IEEE",
    note="Added the operating set-point. Cited as reference [8]; not yet checked directly.",
)

TP = r"\frac{1}{1 + \exp\left(-a_{sp} (P - P_{sp}(t))\right)}"

AUXILIARIES = [
    Auxiliary("T_p", TP, "parasympathetic activity (logistic tone)"),
    Auxiliary("T_s", r"1 - T_p", "sympathetic activity"),
]

# The reference treats Psp as unknown and estimates it from data, so it doesn't
# publish initial conditions or inputs. Those, and every slider range, are ours.
STATES = [
    StateVariable(
        key="H", symbol="H", description="heart rate", unit="s^-1", initial=100 / 60,
        equation=r"\frac{\beta_H T_s}{1 + \gamma T_p} - V_H T_p + \delta_H (H_0 - H)",
        ui_range=(0.5, 4.0), ui_range_basis="heuristic: 30 to 240 beats/min", lower=0,
    ),
    StateVariable(
        key="P", symbol="P", description="mean arterial blood pressure", unit="mmHg", initial=90,
        equation=r"-\frac{P}{R_a (1 + \alpha T_s) C_a} + \frac{H \Delta V}{C_a}",
        ui_range=(40, 200), ui_range_basis="heuristic: physiological pressure range", lower=0,
    ),
]

RANGE_BASIS = "heuristic: 0.1x to 10x the default"


def _parameter(key, symbol, default, unit, description, positive=False, **extra):
    return Parameter(
        key=key, default=default, symbol=symbol, description=description, unit=unit,
        ui_range=(default / 10, default * 10), ui_range_basis=RANGE_BASIS, scale=Scale.LOG,
        lower=0, lower_exclusive=positive, source="Table I", reference="sadoun2026", **extra,
    )


PARAMETERS = [
    _parameter(
        "Ra", "R_a", 0.6, "mmHg·g·ml^-1", "minimum arterial resistance", positive=True,
        note=(
            "Unit as printed. Eq. (4) needs mmHg·s·ml^-1, the unit Ottesen (1997) gives for "
            "resistance. Ottesen's peripheral resistance is 1.05, not 0.6."
        ),
    ),
    _parameter(
        "Ca", "C_a", 1.55, "mmHg·g^-1", "arterial compliance", positive=True,
        note="Unit as printed. Eq. (4) needs ml·mmHg^-1, and Ottesen (1997) gives 1.55 ml/mmHg.",
    ),
    _parameter(
        "dV", r"\Delta V", 50, "ml", "stroke volume",
        note="Ottesen (1997) uses 67.9 ml at rest.",
    ),
    _parameter("betaH", r"\beta_H", 0.84, "s^-2", "sympathetic control of HR"),
    _parameter(
        "gamma", r"\gamma", 0.2, None, "parasympathetic damping of betaH",
        note=(
            "Table I prints no unit; Eq. (3) implies it's dimensionless. "
            "Ottesen (1997) uses 0 for the equivalent parameter."
        ),
    ),
    _parameter(
        "VH", "V_H", 1.17, "s^-1", "parasympathetic control of HR",
        note="Unit as printed. Eq. (3) needs s^-2, and Ottesen (1997) gives 1.17 sec^-2.",
    ),
    _parameter(
        "deltaH", r"\delta_H", 1.7, "s^-1", "relaxation time",
        note="Table I calls this a time, but its unit and its role in Eq. (3) are those of a rate.",
    ),
    _parameter(
        "H0", "H_0", 100 / 60, "s^-1", "intrinsic HR", positive=True,
        note="Table I gives 100 min^-1; converted because this model's time is in seconds.",
    ),
    _parameter(
        "alpha", r"\alpha", 1.3, None, "sympathetic effect on Ra",
        note="Table I prints no unit; Eq. (4) implies it's dimensionless.",
    ),
    _parameter(
        "asp", "a_{sp}", 0.07, "mmHg^-1", "slope of sigmoid determining baroreflex response",
        positive=True,
    ),
]

PSP_RESTING = Signal.constant(90.0)
PSP_STEP = Signal.step(before=90.0, after=110.0, at=30.0)

INPUTS = [
    InputVariable(
        key="Psp", symbol="P_{sp}", description="baroreflex set-point", unit="mmHg",
        default=PSP_RESTING, ui_range=(60, 160), ui_range_basis="heuristic",
    ),
]

# Callables kept for code that passes inputs directly to the engine.
psp_constant = PSP_RESTING.to_callable()  # 90 mmHg
psp_profile = PSP_STEP.to_callable()  # 90 mmHg, stepping up to 110 mmHg at t=30 s

RESTING = Scenario(
    name="resting",
    description="Constant 90 mmHg set-point. P settles from 90 to about 85.7 mmHg within 10 s.",
    t_span=(0, 60),
    n_points=600,
    inputs={"Psp": PSP_RESTING},
    source=None,
)

SETPOINT_STEP = Scenario(
    name="setpoint_step",
    description="Set-point of 90 mmHg, stepping up to 110 mmHg at t=30 s.",
    t_span=(0, 60),
    n_points=600,
    inputs={"Psp": PSP_STEP},
    solver_overrides={"max_step": 0.1},  # keep the solver from stepping over the jump in Psp
    source=None,
)

ENTRY = ModelEntry(
    id="brain_heart",
    name="Autonomic cardiac regulation (delay-free baroreflex model)",
    version="1.0.0",
    domain="physiology",
    tags=("baroreflex", "cardiovascular", "autonomic nervous system"),
    summary=(
        "Delay-free baroreflex model of heart rate H and mean arterial pressure P "
        "driven by the set-point Psp(t)."
    ),
    rhs=rhs,
    time_unit="s",
    states=STATES,
    parameters=PARAMETERS,
    inputs=INPUTS,
    auxiliaries=AUXILIARIES,
    scenarios=[RESTING, SETPOINT_STEP],
    default_scenario="resting",
    references=[SADOUN_2026, OTTESEN_1997, ATAEE_2010],
    equations_reference="sadoun2026",
    status=ModelStatus.NUMERICALLY_VERIFIED,
    evidence=(
        "Matches a separate RK4 implementation of Eqs. (3)-(4), typed separately from model.py, within 1e-5 relative error. Both were written in this project, so this doesn't replace an external implementation."
    ),
    todo=(
        "Units of gamma and alpha aren't printed in Table I (dimensionless by Eqs. 3-4).",
        "Table I units for Ra, Ca and VH don't fit the equations; recorded as printed, see notes.",
        "Initial conditions, Psp signals and scenarios aren't from the paper.",
        "Slider ranges (ui_range) are heuristic, not from a source.",
    ),
)
