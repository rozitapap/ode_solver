"""
Bioprocess ODE Model
======================

5-state cell-culture kinetic model:
    VCD - Viable Cell Density
    Glc - Glucose concentration
    Gln - Glutamine concentration
    Lac - Lactate concentration
    Amm - Ammonia concentration

Rates (mu, qGlc, qGln, qLac, qAmm) are constant parameters rather than
being predicted by a model. No substrate-limiting feedback, so VCD grows
as exp(mu*t) unconditionally -- keep mu small and qGlc/qGln low enough
that consumption doesn't outrun the initial Glc/Gln pools over t_span.
"""

from .base import ODESystemSpec


def rhs(t, X, mu, qGlc, qGln, qLac, qAmm):
    """RHS of the batch mass-balance ODE system."""
    VCD, Glc, Gln, Lac, Amm = X

    dVCDdt = mu * VCD
    dGlcdt = -qGlc * VCD
    dGlndt = -qGln * VCD
    dLacdt = qLac * VCD
    dAmmdt = qAmm * VCD

    return [dVCDdt, dGlcdt, dGlndt, dLacdt, dAmmdt]


MODEL = ODESystemSpec(
    name="bioprocess",
    state_names=["VCD", "Glc", "Gln", "Lac", "Amm"],
    param_names=["mu", "qGlc", "qGln", "qLac", "qAmm"],
    rhs=rhs,
    default_params={
        "mu": 0.02,
        "qGlc": 0.015,
        "qGln": 0.002,
        "qLac": 0.01,
        "qAmm": 0.001,
    },
    default_initial_conditions={"VCD0": 0.5, "Glc0": 30.0, "Gln0": 4.0, "Lac0": 0.1, "Amm0": 0.1},
    default_t_span=(0, 200),
    default_n_points=201,
)
