"""
The universal structure of a catalog entry.

A ModelEntry holds everything a person or a website needs to show, configure
and run a model: identity, description, citations, declared states and
parameters with units and bounds, LaTeX equations, scenarios, solver hints and
reference data. The engine's ODESystemSpec is built from it, so nothing is
written twice.

The engine never imports this package, so it stays model-agnostic.
"""

import csv
import inspect
import math
import re
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Callable, Mapping, Optional, Sequence

from ode_solver.spec import ODESystemSpec, SimulationRequest

from .latex import symbols_used

ID_PATTERN = re.compile(r"^[a-z][a-z0-9_]*$")
VERSION_PATTERN = re.compile(r"^\d+\.\d+\.\d+$")
SOLVER_METHODS = {"RK45", "RK23", "DOP853", "Radau", "BDF", "LSODA"}
DEFAULT_ALIAS = "default"


class ModelStatus(Enum):
    """Empirical quality control finding, using the BioModels terms."""

    # Simulation runs without crashing.
    EXECUTABLE = "executable"
    # Results agree across independent implementations or against a trusted
    # reference within defined tolerances.
    NUMERICALLY_VERIFIED = "numerically verified"
    # Simulation reproduces a documented, published result.
    SCIENTIFICALLY_REPRODUCED = "scientifically reproduced"


class Scale(Enum):
    """How a UI slider maps position to value."""

    LINEAR = "linear"
    LOG = "log"


class DataPurpose(Enum):
    """What a reference dataset is for."""

    VALIDATION = "validation"  # compared against in test_validation.py
    DISPLAY = "display"  # shown alongside simulations


@dataclass(frozen=True)
class Reference:
    """
    A publication. `citation` is the full citation as recorded; the structured
    fields are optional, and only for details confirmed from the paper itself.
    """

    id: str  # local key, e.g. "sadoun2026"; referenced by equations_reference and Parameter.reference
    citation: str
    authors: tuple = ()
    title: Optional[str] = None
    year: Optional[int] = None
    venue: Optional[str] = None
    doi: Optional[str] = None
    url: Optional[str] = None
    note: str = ""


@dataclass(frozen=True)
class StateVariable:
    """One state of the ODE system, in rhs order."""

    key: str  # "S"; engine initial-condition key is key + "0"
    symbol: str  # LaTeX
    description: str
    unit: Optional[str]  # None if it isn't recorded
    initial: float  # default initial value
    equation: str  # LaTeX right-hand side of d(key)/dt
    ui_range: tuple  # (min, max) for the initial-value slider
    ui_range_basis: str  # why that range, e.g. "heuristic: 0 to 10x the default"
    lower: Optional[float] = None  # hard bound on the initial value, enforced on every request
    upper: Optional[float] = None
    scale: Scale = Scale.LINEAR


@dataclass(frozen=True)
class Parameter:
    """One parameter, its default value, bounds and where it came from."""

    key: str
    default: float
    symbol: str  # LaTeX
    description: str
    unit: Optional[str]  # as the source prints it; None if it isn't recorded
    ui_range: tuple  # (min, max) for the slider
    ui_range_basis: str  # why that range, e.g. "heuristic: 0.1x to 10x the default"
    scale: Scale = Scale.LINEAR
    lower: Optional[float] = None  # hard validity bound, enforced on every request
    upper: Optional[float] = None
    lower_exclusive: bool = False  # True means the value must be strictly greater than lower
    source: Optional[str] = None  # e.g. "Table I"; None if it isn't from a reference
    reference: Optional[str] = None  # Reference.id the source refers to
    note: str = ""  # conversions, or anything that differs from the source


@dataclass(frozen=True)
class Auxiliary:
    """A named sub-expression used in the displayed equations, e.g. N = S + I + R + D."""

    symbol: str  # LaTeX
    equation: str  # LaTeX
    description: str = ""


@dataclass(frozen=True)
class Signal:
    """
    A declarative time-varying input, so a website can show and edit it.

    constant:            values=(v,)
    step:                values=(before, after), times=(t_step,)
    piecewise_constant:  values=(v0, v1, ..., vn), times=(t1, ..., tn), increasing
    """

    kind: str
    values: tuple
    times: tuple = ()

    @classmethod
    def constant(cls, value):
        return cls("constant", (value,))

    @classmethod
    def step(cls, before, after, at):
        return cls("step", (before, after), (at,))

    def problems(self):
        if self.kind == "constant":
            ok = len(self.values) == 1 and not self.times
        elif self.kind == "step":
            ok = len(self.values) == 2 and len(self.times) == 1
        elif self.kind == "piecewise_constant":
            ok = len(self.values) == len(self.times) + 1 and list(self.times) == sorted(self.times)
        else:
            return [f"unknown signal kind {self.kind!r}"]
        return [] if ok else [f"malformed {self.kind} signal: values={self.values}, times={self.times}"]

    def to_callable(self) -> Callable:
        """Return f(t). Each value applies from its switching time onwards."""
        values, times = self.values, self.times
        if self.kind == "constant":
            value = values[0]
            return lambda t: value
        if self.kind == "step":
            before, after = values
            at = times[0]
            return lambda t: before if t < at else after

        def piecewise(t):
            for switch, value in zip(times, values):
                if t < switch:
                    return value
            return values[-1]

        return piecewise

    def to_dict(self):
        return {"kind": self.kind, "values": list(self.values), "times": list(self.times)}


@dataclass(frozen=True)
class InputVariable:
    """A time-varying input passed to rhs as a function of t."""

    key: str
    symbol: str  # LaTeX
    description: str
    unit: Optional[str]
    default: Signal
    ui_range: tuple  # (min, max) for the signal's values
    ui_range_basis: str


@dataclass(frozen=True)
class ConservedQuantity:
    """A sum of states that the dynamics keep constant."""

    name: str
    states: tuple


@dataclass(frozen=True)
class SolverHints:
    """Default solver settings, and limits for a web backend."""

    method: str = "RK45"
    rtol: float = 1e-8
    atol: float = 1e-10
    max_step: float = math.inf
    stiff: bool = False
    max_points: int = 5000  # largest n_points a backend should accept


@dataclass(frozen=True)
class ReferenceDataset:
    """
    Data to compare simulations with, stored under the model's reference_data/.

    `path` is either a CSV file, or a digitized-curve folder in the format of
    ode_catalog.validation.load_reference_curve (metadata.json and run*.csv,
    with time_column "x" and columns {"y": state}).
    """

    id: str
    path: str  # relative to reference_data/
    description: str
    purpose: frozenset  # of DataPurpose
    time_column: str
    time_unit: Optional[str]
    columns: Mapping[str, str]  # data column -> state key
    column_units: Mapping[str, Optional[str]]  # data column -> unit
    scenario: Optional[str] = None  # the scenario this data should be compared with
    reference: Optional[str] = None  # Reference.id the data came from


@dataclass(frozen=True)
class Scenario:
    """A named preset. Anything not overridden comes from the entry's defaults."""

    name: str
    description: str
    t_span: tuple
    n_points: int = 300
    initial_overrides: Mapping[str, float] = field(default_factory=dict)  # by state key
    param_overrides: Mapping[str, float] = field(default_factory=dict)
    inputs: Mapping[str, Signal] = field(default_factory=dict)
    solver_overrides: Mapping[str, Any] = field(default_factory=dict)  # SolverHints fields
    source: Optional[str] = None  # e.g. "Fig. 3"; None if it isn't from a reference


@dataclass(frozen=True)
class ModelEntry:
    """One model in the catalog. See ADDING_A_MODEL.md for what each field means."""

    id: str
    name: str
    version: str
    domain: str
    summary: str
    rhs: Callable
    time_unit: Optional[str]
    states: Sequence[StateVariable]
    parameters: Sequence[Parameter]
    scenarios: Sequence[Scenario]
    default_scenario: str
    tags: tuple = ()
    inputs: Sequence[InputVariable] = ()
    auxiliaries: Sequence[Auxiliary] = ()
    conserved_quantities: Sequence[ConservedQuantity] = ()
    solver: SolverHints = SolverHints()
    references: Sequence[Reference] = ()
    equations_reference: Optional[str] = None  # Reference.id; None means not from a publication
    reference_data: Sequence[ReferenceDataset] = ()
    status: ModelStatus = ModelStatus.EXECUTABLE
    evidence: str = ""  # what supports a status above executable, and its tolerance
    todo: tuple = ()  # values that couldn't be verified; exported so a UI can flag them
    extra_checks: Optional[Callable] = None  # check(y0: dict, params: dict) -> list[str]
    folder: Optional[Path] = None  # defaults to the folder of the module defining rhs
    spec: ODESystemSpec = field(init=False, repr=False, compare=False)

    def __post_init__(self):
        problems = validate_entry(self)
        if problems:
            raise ValueError(f"Invalid catalog entry '{self.id}': {problems}")
        object.__setattr__(self, "spec", self._build_spec())

    # -- lookups -----------------------------------------------------------

    @property
    def model_folder(self) -> Path:
        return Path(self.folder) if self.folder else Path(inspect.getfile(self.rhs)).parent

    @property
    def reference(self) -> Optional[Reference]:
        """The reference the equations come from, or None."""
        return _by_id(self.references).get(self.equations_reference)

    @property
    def related_references(self):
        return [r for r in self.references if r.id != self.equations_reference]

    @property
    def state_units(self) -> dict:
        return {state.key: state.unit for state in self.states}

    @property
    def units(self) -> dict:
        """Keyword arguments giving the plotting functions this model's units."""
        return {"time_unit": self.time_unit, "state_units": self.state_units}

    @property
    def long_description(self) -> str:
        """The model's README.md, which is the single source for the long description."""
        readme = self.model_folder / "README.md"
        return readme.read_text(encoding="utf-8") if readme.exists() else ""

    def scenario(self, name: Optional[str] = None) -> Scenario:
        """Look up a scenario by name; None or "default" gives the default scenario."""
        by_name = {scenario.name: scenario for scenario in self.scenarios}
        if name in by_name:
            return by_name[name]
        if name is None or name == DEFAULT_ALIAS:
            return by_name[self.default_scenario]
        available = ", ".join(by_name)
        raise KeyError(f"Unknown scenario '{name}' for '{self.id}'. Available: {available}")

    def request(self, name: Optional[str] = None) -> SimulationRequest:
        """Build the SimulationRequest for a scenario, or the default one."""
        scenario = self.scenario(name)
        initial = {state.key: state.initial for state in self.states}
        initial.update(scenario.initial_overrides)
        params = {parameter.key: parameter.default for parameter in self.parameters}
        params.update(scenario.param_overrides)
        signals = {variable.key: variable.default for variable in self.inputs}
        signals.update(scenario.inputs)
        solver = {**_solver_fields(self.solver), **scenario.solver_overrides}
        return SimulationRequest(
            system_name=self.id,
            params=params,
            initial_conditions={key + "0": value for key, value in initial.items()},
            t_span=tuple(scenario.t_span),
            n_points=scenario.n_points,
            method=solver["method"],
            rtol=solver["rtol"],
            atol=solver["atol"],
            max_step=solver["max_step"],
            inputs={key: signal.to_callable() for key, signal in signals.items()},
        )

    # -- engine spec -------------------------------------------------------

    def _build_spec(self) -> ODESystemSpec:
        default = self.request()
        all_states = {state.key for state in self.states}
        return ODESystemSpec(
            name=self.id,
            state_names=[state.key for state in self.states],
            param_names=[parameter.key for parameter in self.parameters],
            rhs=self.rhs,
            description=self.summary,
            invariant_check=self._check_request_values,
            input_names=[variable.key for variable in self.inputs],
            conserved_total=any(set(q.states) == all_states for q in self.conserved_quantities),
            default_params=default.params,
            default_initial_conditions=default.initial_conditions,
            default_t_span=default.t_span,
            default_n_points=default.n_points,
            default_inputs={variable.key: variable.default.to_callable() for variable in self.inputs},
        )

    def _check_request_values(self, y0, params):
        """The engine's invariant check: declared hard bounds, then extra_checks."""
        problems = []
        for state in self.states:
            key = state.key + "0"
            if key in y0:
                problems += _bound_problems(f"Initial condition {key}", y0[key], state.lower, state.upper)
        for parameter in self.parameters:
            if parameter.key in params:
                problems += _bound_problems(
                    f"Parameter {parameter.key}", params[parameter.key],
                    parameter.lower, parameter.upper, parameter.lower_exclusive,
                )
        if not problems and self.extra_checks is not None:
            problems += self.extra_checks(y0, params)
        return problems

    # -- export ------------------------------------------------------------

    def summary_dict(self) -> dict:
        """The short record for a "list models" response."""
        return {
            "id": self.id,
            "name": self.name,
            "version": self.version,
            "domain": self.domain,
            "tags": list(self.tags),
            "status": self.status.value,
            "summary": self.summary,
            "n_states": len(self.states),
            "n_parameters": len(self.parameters),
            "has_reference": self.equations_reference is not None,
            "has_reference_data": bool(self.reference_data),
        }

    def to_dict(self) -> dict:
        """Everything except callables, as JSON-compatible values."""
        return {
            **self.summary_dict(),
            "long_description": self.long_description,
            "time_unit": self.time_unit,
            "states": [_plain(state) for state in self.states],
            "auxiliaries": [_plain(aux) for aux in self.auxiliaries],
            "parameters": [_plain(parameter) for parameter in self.parameters],
            "inputs": [_plain(variable) for variable in self.inputs],
            "conserved_quantities": [_plain(q) for q in self.conserved_quantities],
            "scenarios": [_plain(scenario) for scenario in self.scenarios],
            "default_scenario": self.default_scenario,
            "solver": _plain(self.solver),
            "references": [_plain(reference) for reference in self.references],
            "equations_reference": self.equations_reference,
            "reference_data": [_plain(dataset) for dataset in self.reference_data],
            "evidence": self.evidence,
            "todo": list(self.todo),
        }


# -- validation --------------------------------------------------------------


def validate_entry(entry: ModelEntry) -> list:
    """Return every problem with an entry; empty if it's valid."""
    problems = []
    add = problems.append

    # Identity
    if not ID_PATTERN.match(entry.id):
        add(f"id '{entry.id}' must match {ID_PATTERN.pattern}")
    if not VERSION_PATTERN.match(entry.version):
        add(f"version '{entry.version}' must look like 1.0.0")
    for name in ("name", "domain", "summary"):
        if not getattr(entry, name):
            add(f"{name} is empty")

    states = {state.key: state for state in entry.states}
    params = {parameter.key: parameter for parameter in entry.parameters}
    inputs = {variable.key: variable for variable in entry.inputs}
    refs = _by_id(entry.references)

    # Unique keys and symbols
    keys = list(states) + list(params) + list(inputs)
    if len(entry.states) + len(entry.parameters) + len(entry.inputs) != len(set(keys)):
        add(f"state, parameter and input keys must be unique: {_duplicates(_all_keys(entry))}")
    symbols = [x.symbol for x in (*entry.states, *entry.parameters, *entry.inputs, *entry.auxiliaries)]
    if len(symbols) != len(set(symbols)):
        add(f"symbols must be unique: {_duplicates(symbols)}")
    if len(refs) != len(entry.references):
        add("reference ids must be unique")

    # rhs signature matches the declared parameters and inputs
    try:
        arguments = list(inspect.signature(entry.rhs).parameters)[2:]
        if sorted(arguments) != sorted(list(params) + list(inputs)):
            add(f"rhs arguments {arguments} don't match parameters and inputs {list(params) + list(inputs)}")
    except (TypeError, ValueError):
        add("rhs must be a Python function rhs(t, y, **params)")

    # States
    for state in entry.states:
        if not state.description:
            add(f"state {state.key} has no description")
        if not state.equation.strip():
            add(f"state {state.key} has no equation")
        problems += _range_problems(f"state {state.key}", state.initial, state.ui_range,
                                    state.ui_range_basis, state.scale, state.lower, state.upper)

    # Parameters
    for parameter in entry.parameters:
        label = f"parameter {parameter.key}"
        if not parameter.description:
            add(f"{label} has no description")
        problems += _range_problems(label, parameter.default, parameter.ui_range,
                                    parameter.ui_range_basis, parameter.scale,
                                    parameter.lower, parameter.upper, parameter.lower_exclusive)
        if parameter.source and parameter.reference is None:
            add(f"{label} has a source but no reference")
        if parameter.reference is not None and parameter.reference not in refs:
            add(f"{label} cites unknown reference '{parameter.reference}'")

    # Inputs
    for variable in entry.inputs:
        problems += [f"input {variable.key}: {p}" for p in _signal_problems(variable.default, variable)]
        if not variable.ui_range_basis:
            add(f"input {variable.key} has no ui_range_basis")

    # Equations reference declared symbols, and every parameter and input is used
    declared = symbols
    used = set()
    equations = [(f"d{s.key}/dt", s.equation) for s in entry.states]
    equations += [(a.symbol, a.equation) for a in entry.auxiliaries]
    for label, equation in equations:
        found, latex_problems = symbols_used(equation, declared)
        used |= found
        problems += [f"equation {label}: {p}" for p in latex_problems]
    for x in (*entry.parameters, *entry.inputs, *entry.auxiliaries):
        if x.symbol not in used:
            add(f"symbol {x.symbol!r} ({getattr(x, 'key', 'auxiliary')}) appears in no equation")

    # Conserved quantities
    for quantity in entry.conserved_quantities:
        unknown = [key for key in quantity.states if key not in states]
        if unknown or not quantity.states:
            add(f"conserved quantity '{quantity.name}' has unknown states {unknown}")

    # Scenarios
    names = [scenario.name for scenario in entry.scenarios]
    if len(set(names)) != len(names):
        add(f"duplicate scenario names: {names}")
    if entry.default_scenario not in names:
        add(f"default_scenario '{entry.default_scenario}' isn't a scenario")
    solver_fields = set(_solver_fields(entry.solver))
    for scenario in entry.scenarios:
        label = f"scenario {scenario.name}"
        if not scenario.description:
            add(f"{label} has no description")
        if len(scenario.t_span) != 2 or not scenario.t_span[1] > scenario.t_span[0]:
            add(f"{label}: t_span must be increasing, got {scenario.t_span}")
        if not 2 <= scenario.n_points <= entry.solver.max_points:
            add(f"{label}: n_points {scenario.n_points} outside [2, {entry.solver.max_points}]")
        for key, value in scenario.initial_overrides.items():
            if key not in states:
                add(f"{label} sets unknown state '{key}'")
            else:
                s = states[key]
                problems += _value_problems(f"{label}: {key}", value, s.ui_range, s.lower, s.upper)
        for key, value in scenario.param_overrides.items():
            if key not in params:
                add(f"{label} sets unknown parameter '{key}'")
            else:
                p = params[key]
                problems += _value_problems(f"{label}: {key}", value, p.ui_range, p.lower, p.upper,
                                            p.lower_exclusive)
        for key, signal in scenario.inputs.items():
            if key not in inputs:
                add(f"{label} sets unknown input '{key}'")
            else:
                problems += [f"{label}: input {key}: {p}" for p in _signal_problems(signal, inputs[key])]
        unknown = set(scenario.solver_overrides) - solver_fields
        if unknown:
            add(f"{label} sets unknown solver fields {sorted(unknown)}")

    # Solver
    methods = {entry.solver.method} | {s.solver_overrides.get("method", entry.solver.method) for s in entry.scenarios}
    if methods - SOLVER_METHODS:
        add(f"unknown solver methods {sorted(methods - SOLVER_METHODS)}")

    # References and status
    if entry.equations_reference is not None and entry.equations_reference not in refs:
        add(f"equations_reference '{entry.equations_reference}' isn't in references")
    if entry.references and entry.equations_reference is None:
        add("equations_reference must name one of the references")
    if entry.status is not ModelStatus.EXECUTABLE and not entry.evidence:
        add(f"status '{entry.status.value}' needs evidence")
    if entry.status is ModelStatus.SCIENTIFICALLY_REPRODUCED and entry.equations_reference is None:
        add("a scientifically reproduced model needs a reference")

    # Reference data
    if entry.reference_data:
        problems += _reference_data_problems(entry, states, refs, set(names))

    return problems


def _reference_data_problems(entry, states, refs, scenario_names):
    problems = []
    root = entry.model_folder / "reference_data"
    for dataset in entry.reference_data:
        label = f"reference data '{dataset.id}'"
        if not dataset.purpose or not all(isinstance(p, DataPurpose) for p in dataset.purpose):
            problems.append(f"{label}: purpose must be a non-empty set of DataPurpose")
        for column, state in dataset.columns.items():
            if state not in states:
                problems.append(f"{label}: column '{column}' maps to unknown state '{state}'")
            if column not in dataset.column_units:
                problems.append(f"{label}: column '{column}' has no unit entry")
        if dataset.scenario is not None and dataset.scenario not in scenario_names:
            problems.append(f"{label}: unknown scenario '{dataset.scenario}'")
        if dataset.reference is not None and dataset.reference not in refs:
            problems.append(f"{label}: unknown reference '{dataset.reference}'")

        path = root / dataset.path
        if path.is_dir():
            if not (path / "metadata.json").exists() or not list(path.glob("run*.csv")):
                problems.append(f"{label}: folder needs metadata.json and run*.csv files")
            elif dataset.time_column != "x" or set(dataset.columns) != {"y"}:
                problems.append(f"{label}: a digitized-curve folder uses time_column 'x' and column 'y'")
        elif path.is_file():
            with path.open(newline="", encoding="utf-8") as handle:
                header = next(csv.reader(handle), [])
            missing = [c for c in [dataset.time_column, *dataset.columns] if c not in header]
            if missing:
                problems.append(f"{label}: columns {missing} aren't in {path.name}")
        else:
            problems.append(f"{label}: {path} doesn't exist")
    return problems


# -- helpers -----------------------------------------------------------------


def _by_id(items):
    return {item.id: item for item in items}


def _all_keys(entry):
    return [x.key for x in (*entry.states, *entry.parameters, *entry.inputs)]


def _duplicates(values):
    return sorted({v for v in values if list(values).count(v) > 1})


def _solver_fields(solver):
    return {"method": solver.method, "rtol": solver.rtol, "atol": solver.atol, "max_step": solver.max_step}


def _bound_problems(label, value, lower, upper, lower_exclusive=False):
    if lower is not None:
        if lower_exclusive and not value > lower:
            return [f"{label} must be greater than {lower:g}, got {value:g}."]
        if not lower_exclusive and not value >= lower:
            return [f"{label} must be at least {lower:g}, got {value:g}."]
    if upper is not None and not value <= upper:
        return [f"{label} must be at most {upper:g}, got {value:g}."]
    return []


def _value_problems(label, value, ui_range, lower, upper, lower_exclusive=False):
    problems = _bound_problems(label, value, lower, upper, lower_exclusive)
    if ui_range and len(ui_range) == 2 and not ui_range[0] <= value <= ui_range[1]:
        problems.append(f"{label} = {value:g} is outside ui_range {tuple(ui_range)}")
    return problems


def _range_problems(label, default, ui_range, basis, scale, lower, upper, lower_exclusive=False):
    problems = []
    if len(ui_range) != 2 or not ui_range[0] < ui_range[1]:
        return [f"{label}: ui_range must be (min, max) with min < max, got {ui_range}"]
    if not basis:
        problems.append(f"{label} has no ui_range_basis")
    if scale is Scale.LOG and ui_range[0] <= 0:
        problems.append(f"{label}: a log-scale ui_range must be positive")
    for end in ui_range:
        problems += _bound_problems(f"{label}: ui_range end", end, lower, upper, lower_exclusive)
    problems += _value_problems(f"{label}: default", default, ui_range, lower, upper, lower_exclusive)
    return problems


def _signal_problems(signal, variable):
    problems = signal.problems()
    low, high = variable.ui_range
    outside = [v for v in signal.values if not low <= v <= high]
    if outside:
        problems.append(f"values {outside} outside ui_range {tuple(variable.ui_range)}")
    return problems


def _plain(value):
    """Convert dataclasses, enums, tuples and inf into JSON-compatible values."""
    if isinstance(value, Signal):
        return value.to_dict()
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if hasattr(value, "__dataclass_fields__"):
        return {name: _plain(getattr(value, name)) for name in value.__dataclass_fields__}
    if isinstance(value, Mapping):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (frozenset, set)):
        return sorted(_plain(v) for v in value)
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    return value
