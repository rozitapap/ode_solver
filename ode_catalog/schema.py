"""
What a catalog entry records about a model, on top of the ODESystemSpec the
engine needs: its source, parameter provenance, runnable scenarios, and how
far it has been validated.

The engine never imports this package, so it stays model-agnostic.
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Mapping, Optional, Sequence

import numpy as np

from ode_solver.spec import ODESystemSpec, SimulationRequest


class ModelStatus(Enum):
    """Empirical quality control finding, using the BioModels terms."""

    # Simulation runs without crashing.
    EXECUTABLE = "executable"
    # Results agree across independent implementations or against a trusted
    # reference within defined tolerances.
    NUMERICALLY_VERIFIED = "numerically verified"
    # Simulation reproduces a documented, published result.
    SCIENTIFICALLY_REPRODUCED = "scientifically reproduced"


@dataclass(frozen=True)
class Reference:
    """A publication. Leave out anything you can't confirm from the paper itself."""

    citation: str
    doi: Optional[str] = None
    url: Optional[str] = None
    note: str = ""


@dataclass(frozen=True)
class Parameter:
    """One parameter value and where it came from."""

    name: str
    value: float
    unit: Optional[str] = None
    description: str = ""
    source: Optional[str] = None  # e.g. "Table I"; None if it isn't from the reference
    note: str = ""  # conversions, or anything that differs from the source


@dataclass(frozen=True)
class Scenario:
    """Everything needed to repeat one simulation."""

    name: str
    description: str
    initial_conditions: Mapping[str, float]
    t_span: tuple
    n_points: int = 300
    inputs: Mapping = field(default_factory=dict)
    param_overrides: Mapping[str, float] = field(default_factory=dict)
    method: str = "RK45"
    rtol: float = 1e-8
    atol: float = 1e-10
    max_step: float = np.inf
    source: Optional[str] = None  # e.g. "Fig. 3"; None if it isn't from the reference


@dataclass(frozen=True)
class ModelEntry:
    """One model in the catalog."""

    spec: ODESystemSpec
    title: str
    parameters: Sequence[Parameter]
    scenarios: Sequence[Scenario]
    default_scenario: str
    time_unit: Optional[str]  # None if it isn't recorded
    state_units: Mapping[str, Optional[str]]  # one entry per state; None if it isn't recorded
    status: ModelStatus = ModelStatus.EXECUTABLE
    reference: Optional[Reference] = None  # None means the model isn't from a publication
    related_references: Sequence[Reference] = ()
    evidence: str = ""  # what supports a status above executable, and its tolerance

    def __post_init__(self):
        problems = []

        names = [parameter.name for parameter in self.parameters]
        if sorted(names) != sorted(self.spec.param_names):
            problems.append(f"parameters {names} don't match spec.param_names")
        values = {parameter.name: parameter.value for parameter in self.parameters}
        if values != dict(self.spec.default_params):
            problems.append("parameter values don't match spec.default_params")

        if sorted(self.state_units) != sorted(self.spec.state_names):
            problems.append("state_units must have one entry per state")

        scenario_names = [scenario.name for scenario in self.scenarios]
        if len(set(scenario_names)) != len(scenario_names):
            problems.append(f"duplicate scenario names: {scenario_names}")
        if self.default_scenario not in scenario_names:
            problems.append(f"default_scenario '{self.default_scenario}' isn't a scenario")

        if self.status is not ModelStatus.EXECUTABLE and not self.evidence:
            problems.append(f"status '{self.status.value}' needs evidence")

        if self.status is ModelStatus.SCIENTIFICALLY_REPRODUCED and self.reference is None:
            problems.append("a scientifically reproduced model needs a reference")

        if problems:
            raise ValueError(f"Invalid catalog entry '{self.spec.name}': {problems}")

    @property
    def units(self) -> dict:
        """Keyword arguments giving the plotting functions this model's units."""
        return {"time_unit": self.time_unit, "state_units": dict(self.state_units)}

    def scenario(self, name: Optional[str] = None) -> Scenario:
        """Look up a scenario by name, or the default one."""
        name = name or self.default_scenario
        for scenario in self.scenarios:
            if scenario.name == name:
                return scenario
        available = ", ".join(s.name for s in self.scenarios)
        raise KeyError(f"Unknown scenario '{name}' for '{self.spec.name}'. Available: {available}")

    def request(self, name: Optional[str] = None) -> SimulationRequest:
        """Build the SimulationRequest for a scenario, or the default one."""
        scenario = self.scenario(name)
        return SimulationRequest(
            system_name=self.spec.name,
            params={**self.spec.default_params, **scenario.param_overrides},
            initial_conditions=dict(scenario.initial_conditions),
            t_span=scenario.t_span,
            n_points=scenario.n_points,
            method=scenario.method,
            rtol=scenario.rtol,
            atol=scenario.atol,
            max_step=scenario.max_step,
            inputs=dict(scenario.inputs),
        )
