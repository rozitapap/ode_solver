"""
Metabolite Pathway Model
======================

Linear 5-species mass-action chain M1 -> M2 -> M3 -> M4 -> M5, each step
first-order with rate k_i.
"""

from ode_solver.spec import ODESystemSpec


def rhs(t, y, k1, k2, k3, k4):
    """RHS of the linear first-order chain."""
    M1, M2, M3, M4, M5 = y

    dM1 = -k1 * M1
    dM2 = k1 * M1 - k2 * M2
    dM3 = k2 * M2 - k3 * M3
    dM4 = k3 * M3 - k4 * M4
    dM5 = k4 * M4

    return [dM1, dM2, dM3, dM4, dM5]


def invariants(y0, params):
    """Initial conditions and rate constants must be non-negative."""
    problems = []
    for key in ("M10", "M20", "M30", "M40", "M50"):
        if y0.get(key, 0) < 0:
            problems.append(f"Initial condition {key} must be non-negative.")
    for key in ("k1", "k2", "k3", "k4"):
        if params.get(key, 0) < 0:
            problems.append(f"Rate constant {key} must be non-negative.")
    return problems


MODEL = ODESystemSpec(
    name="metabolite_pathway",
    state_names=["M1", "M2", "M3", "M4", "M5"],
    param_names=["k1", "k2", "k3", "k4"],
    rhs=rhs,
    description=(
        "Linear 5-metabolite mass-action chain M1->M2->M3->M4->M5. "
        "dM1/dt=-k1*M1, dM2/dt=k1*M1-k2*M2, dM3/dt=k2*M2-k3*M3, "
        "dM4/dt=k3*M3-k4*M4, dM5/dt=k4*M4."
    ),
    invariant_check=invariants,
    conserved_total=True,
    default_params={"k1": 0.50, "k2": 0.30, "k3": 0.20, "k4": 0.10},
    default_initial_conditions={"M10": 100, "M20": 0, "M30": 0, "M40": 0, "M50": 0},
    default_t_span=(0, 50),
    default_n_points=300,
)
