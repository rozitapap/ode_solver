"""
Brain-Heart Model
======================

Baroreflex control of heart rate and arterial pressure.
State H = heart rate, P = mean arterial blood pressure.
Input Psp(t) = baroreflex set-point; output y = H.
Time is in seconds, so H is in beats/s and H0 = 100/60.
Note the parameter H0 (intrinsic heart rate) is separate from the
initial-condition key "H0" (H at t=0); they live in different dicts.
"""

import numpy as np

from .base import ODESystemSpec


def rhs(t, y, Ra, Ca, dV, betaH, gamma, VH, deltaH, H0, alpha, asp, Psp):
    """RHS of the baroreflex ODE system."""
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


# Table I values, time in seconds.
PARAMS = {
    "Ra": 0.6,        # minimum arterial resistance
    "Ca": 1.55,       # arterial compliance
    "dV": 50,         # stroke volume [ml]
    "betaH": 0.84,    # sympathetic control of HR [s^-2]
    "gamma": 0.2,     # parasympathetic damping of betaH
    "VH": 1.17,       # parasympathetic control of HR [s^-1]
    "deltaH": 1.7,    # relaxation rate [s^-1]
    "H0": 100 / 60,   # intrinsic HR: 100 min^-1 in beats/s
    "alpha": 1.3,     # sympathetic effect on Ra
    "asp": 0.07,      # slope of baroreflex sigmoid [mmHg^-1]
}


def psp_constant(t):
    """Resting baroreflex set-point [mmHg]: the default when no Psp is given."""
    return 90.0


def psp_profile(t):
    """Baroreflex set-point [mmHg]: 90 mmHg, stepping up to 110 mmHg at t=30 s."""
    return 90.0 if t < 30.0 else 110.0


MODEL = ODESystemSpec(
    name="brain_heart",
    state_names=["H", "P"],
    param_names=list(PARAMS),
    input_names=["Psp"],
    rhs=rhs,
    description=(
        "Baroreflex model of heart rate H and mean arterial pressure P driven "
        "by the set-point Psp(t)."
    ),
    invariant_check=invariants,
    default_params=PARAMS,
    default_initial_conditions={"H0": 100 / 60, "P0": 90},
    default_t_span=(0, 60),
    default_n_points=600,
    default_inputs={"Psp": psp_constant},
)
