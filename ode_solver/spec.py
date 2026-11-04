"""
Shared ODE contract: how a system is defined, requested, and reported.

Any new ODE system only needs to satisfy this contract to plug into the
solver and everything built on it.
"""

from dataclasses import dataclass, field
from typing import Callable, Optional, Sequence

import numpy as np


@dataclass(frozen=True)
class ODESystemSpec:
    """Definition of a system of ODEs, independent of any particular run.

    `rhs` must have the signature `rhs(t, y, **params) -> Sequence[float]`,
    matching the state ordering in `state_names` and returning derivatives
    in that same order. This is the shape `scipy.integrate.solve_ivp` expects.
    """

    name: str
    state_names: Sequence[str]
    param_names: Sequence[str]
    rhs: Callable
    description: str = ""
    # Optional model-specific invariant check, e.g. conservation of population.
    # Signature: check(y0: dict, params: dict) -> list[str] of problems (empty if OK).
    invariant_check: Optional[Callable] = None
    # Names of time-varying inputs, each passed to `rhs` as a keyword
    # argument holding a function of time (e.g. `Psp(t)`).
    input_names: Sequence[str] = ()
    # True only if the sum of all states is constant for this system. Summing
    # states is meaningless otherwise, so drift is reported only when declared.
    conserved_total: bool = False
    # A runnable base case, used by default_request().
    default_params: dict = field(default_factory=dict)
    default_initial_conditions: dict = field(default_factory=dict)
    default_t_span: Optional[tuple] = None
    default_n_points: int = 300
    default_inputs: dict = field(default_factory=dict)

    def default_request(self) -> "SimulationRequest":
        """Build a SimulationRequest from this spec's default values."""
        if self.default_t_span is None:
            raise ValueError(f"System '{self.name}' has no default_t_span.")
        return SimulationRequest(
            system_name=self.name,
            params=dict(self.default_params),
            initial_conditions=dict(self.default_initial_conditions),
            t_span=self.default_t_span,
            n_points=self.default_n_points,
            inputs=dict(self.default_inputs),
        )


@dataclass
class SimulationRequest:
    """User-facing input: which system, what parameters, what conditions."""

    system_name: str
    params: dict
    initial_conditions: dict
    t_span: tuple
    n_points: int = 300
    method: str = "RK45"
    rtol: float = 1e-8
    atol: float = 1e-10
    max_step: float = np.inf
    # Time-varying inputs keyed by the spec's `input_names`, each a callable f(t).
    inputs: dict = field(default_factory=dict)

    def t_eval(self) -> np.ndarray:
        return np.linspace(self.t_span[0], self.t_span[1], self.n_points)


@dataclass
class SimulationResult:
    """Structured output of a completed (or failed) integration run."""

    system_name: str
    t: np.ndarray
    trajectories: dict
    params: dict
    initial_conditions: dict
    success: bool
    message: str
    metadata: dict = field(default_factory=dict)

    def summary_table(self, n_rows: int = 5) -> str:
        """Return a small preformatted text table of the first/last rows."""
        column_names = ["t"]
        for state_name in self.trajectories:
            column_names.append(state_name)

        header_cells = []
        for column_name in column_names:
            header_cells.append(format(column_name, ">12"))
        lines = ["  ".join(header_cells)]

        total_rows = len(self.t)
        first_count = min(n_rows, total_rows)

        for i in range(first_count):
            row_cells = [f"{self.t[i]:12.3f}"]
            for state_name in self.trajectories:
                row_cells.append(f"{self.trajectories[state_name][i]:12.3f}")
            lines.append("  ".join(row_cells))

        if total_rows > 2 * n_rows:
            dots_cells = []
            for _ in column_names:
                dots_cells.append(format("...", ">12"))
            lines.append("  ".join(dots_cells))

        start = max(total_rows - n_rows, first_count)
        for i in range(start, total_rows):
            row_cells = [f"{self.t[i]:12.3f}"]
            for state_name in self.trajectories:
                row_cells.append(f"{self.trajectories[state_name][i]:12.3f}")
            lines.append("  ".join(row_cells))

        return "\n".join(lines)
