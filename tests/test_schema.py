"""
Every catalog entry has the universal structure, passes validate_entry(), and
round-trips through JSON; and validate_entry() rejects each kind of mistake.
"""

import dataclasses
import json

import numpy as np
import pytest

from ode_catalog import (
    CATALOG,
    Auxiliary,
    DataPurpose,
    InputVariable,
    ModelEntry,
    Parameter,
    Reference,
    ReferenceDataset,
    Scale,
    Scenario,
    Signal,
    StateVariable,
    validate_entry,
)
from ode_catalog._template import model as template
from ode_catalog.export import catalog_summary, catalog_to_dict, main as export_main
from ode_catalog.latex import symbols_used
from ode_solver import solve_ode, validate_request

ENTRIES = {**CATALOG, "template": template.ENTRY}
IDS = list(ENTRIES)


# -- every entry ---------------------------------------------------------------


@pytest.mark.parametrize("entry", ENTRIES.values(), ids=IDS)
def test_entry_passes_validation(entry):
    assert validate_entry(entry) == []


@pytest.mark.parametrize("entry", ENTRIES.values(), ids=IDS)
def test_entry_round_trips_through_json(entry):
    data = entry.to_dict()
    text = json.dumps(data, allow_nan=False)
    assert json.loads(text) == data


@pytest.mark.parametrize("entry", ENTRIES.values(), ids=IDS)
def test_export_has_every_declared_item(entry):
    data = entry.to_dict()
    assert [s["key"] for s in data["states"]] == list(entry.spec.state_names)
    assert [p["key"] for p in data["parameters"]] == list(entry.spec.param_names)
    assert [i["key"] for i in data["inputs"]] == list(entry.spec.input_names)
    assert data["default_scenario"] in [s["name"] for s in data["scenarios"]]
    assert all(s["equation"] and s["symbol"] for s in data["states"])


@pytest.mark.parametrize("entry", CATALOG.values(), ids=list(CATALOG))
def test_long_description_is_the_readme(entry):
    assert entry.long_description.startswith(f"# {entry.name}")


@pytest.mark.parametrize("entry", ENTRIES.values(), ids=IDS)
def test_declared_conserved_quantities_hold(entry):
    result = solve_ode(entry.spec, entry.request())
    for quantity in entry.conserved_quantities:
        total = sum(result.trajectories[key] for key in quantity.states)
        assert total[-1] == pytest.approx(total[0], rel=1e-6), quantity.name


@pytest.mark.parametrize("entry", ENTRIES.values(), ids=IDS)
def test_engine_enforces_declared_bounds(entry):
    bounded = [s for s in entry.states if s.lower is not None]
    if not bounded:
        pytest.skip("no bounded state")
    request = entry.request()
    request.initial_conditions[bounded[0].key + "0"] = bounded[0].lower - 1
    assert any("must be at least" in p for p in validate_request(entry.spec, request))


def test_catalog_export_is_valid_json(tmp_path):
    path = tmp_path / "catalog.json"
    export_main(["-o", str(path)])
    data = json.loads(path.read_text(encoding="utf-8"))
    assert data == catalog_to_dict()
    assert [m["id"] for m in data["models"]] == sorted(CATALOG)


def test_summary_lists_every_model_without_details():
    summary = catalog_summary()
    assert [m["id"] for m in summary] == sorted(CATALOG)
    assert all("states" not in m and "long_description" not in m for m in summary)


def test_catalog_keys_are_entry_ids():
    assert all(key == entry.id for key, entry in CATALOG.items())


# -- template ------------------------------------------------------------------


def test_filled_template_is_a_valid_entry():
    """Filling every placeholder in the template with another model gives a valid entry."""

    def rhs(t, y, r, K):
        (N,) = y
        return [r * N * (1 - N / K)]

    entry = dataclasses.replace(
        template.ENTRY,
        id="logistic",
        name="Logistic growth",
        domain="ecology",
        summary="Logistic population growth.",
        rhs=rhs,
        time_unit="day",
        states=[
            StateVariable(
                key="N", symbol="N", description="population", unit="individuals", initial=10,
                equation=r"r N \left(1 - \frac{N}{K}\right)", ui_range=(0, 1000),
                ui_range_basis="heuristic", lower=0,
            ),
        ],
        parameters=[
            Parameter(key="r", default=0.5, symbol="r", description="growth rate", unit="day^-1",
                      ui_range=(0.05, 5), ui_range_basis="heuristic", scale=Scale.LOG, lower=0),
            Parameter(key="K", default=500, symbol="K", description="carrying capacity",
                      unit="individuals", ui_range=(1, 1000), ui_range_basis="heuristic",
                      lower=0, lower_exclusive=True),
        ],
        conserved_quantities=[],
        scenarios=[Scenario("growth", "Growth from 10 individuals.", (0, 30), n_points=31)],
        default_scenario="growth",
    )
    assert validate_entry(entry) == []
    result = solve_ode(entry.spec, entry.request())
    assert result.success
    assert json.loads(json.dumps(entry.to_dict()))["id"] == "logistic"


# -- signals and LaTeX ---------------------------------------------------------


def test_step_signal_switches_at_its_time():
    f = Signal.step(before=90.0, after=110.0, at=30.0).to_callable()
    assert [f(0), f(29.999), f(30.0), f(60)] == [90.0, 90.0, 110.0, 110.0]


def test_piecewise_signal():
    f = Signal("piecewise_constant", (1.0, 2.0, 3.0), (10.0, 20.0)).to_callable()
    assert [f(0), f(10), f(15), f(20), f(99)] == [1.0, 2.0, 2.0, 3.0, 3.0]


@pytest.mark.parametrize(
    "signal",
    [Signal("ramp", (1.0,)), Signal("step", (1.0,), (2.0,)), Signal("piecewise_constant", (1.0, 2.0), (3.0, 1.0))],
)
def test_malformed_signals_are_reported(signal):
    assert signal.problems()


def test_latex_prefers_the_longest_symbol():
    used, problems = symbols_used(r"P - P_{sp}(t) + H_0", ["P", "P_{sp}", "H", "H_0"])
    assert problems == []
    assert used == {"P", "P_{sp}", "H_0"}


@pytest.mark.parametrize(
    "equation, message",
    [(r"-k A x", "undeclared symbol 'x'"), (r"\sin(A)", r"\\sin"), (r"k A @", "unexpected character")],
)
def test_latex_reports_undeclared_tokens(equation, message):
    _, problems = symbols_used(equation, ["A", "k"])
    assert any(message.replace("\\\\", "\\") in p for p in problems), problems


# -- validate_entry rejects each mistake -----------------------------------------

A, B = template.STATES
K = template.PARAMETERS[0]
BASE = template.BASE_CASE


def make(**overrides):
    return dataclasses.replace(template.ENTRY, **overrides)


@pytest.mark.parametrize(
    "overrides, message",
    [
        ({"id": "Bad-Id"}, "id 'Bad-Id'"),
        ({"version": "1"}, "version"),
        ({"summary": ""}, "summary is empty"),
        ({"states": [A, dataclasses.replace(B, key="A")]}, "keys must be unique"),
        ({"states": [A, dataclasses.replace(B, symbol="A")]}, "symbols must be unique"),
        ({"states": [A, dataclasses.replace(B, equation=" ")]}, "has no equation"),
        ({"states": [A, dataclasses.replace(B, equation="k A + z")]}, "undeclared symbol 'z'"),
        ({"states": [dataclasses.replace(A, equation="-A"), dataclasses.replace(B, equation="A")]},
         "appears in no equation"),
        ({"states": [dataclasses.replace(A, initial=50.0), B]}, "outside ui_range"),
        ({"states": [dataclasses.replace(A, ui_range=(-1, 10)), B]}, "must be at least 0"),
        ({"states": [dataclasses.replace(A, ui_range_basis=""), B]}, "no ui_range_basis"),
        ({"parameters": [dataclasses.replace(K, ui_range=(0, 5))]}, "log-scale ui_range must be positive"),
        ({"parameters": [dataclasses.replace(K, default=99.0)]}, "outside ui_range"),
        ({"parameters": [dataclasses.replace(K, source="Table 1")]}, "source but no reference"),
        ({"parameters": [dataclasses.replace(K, reference="nobody")]}, "unknown reference"),
        ({"scenarios": [dataclasses.replace(BASE, initial_overrides={"C": 1.0})]}, "unknown state 'C'"),
        ({"scenarios": [dataclasses.replace(BASE, param_overrides={"q": 1.0})]}, "unknown parameter 'q'"),
        ({"scenarios": [dataclasses.replace(BASE, param_overrides={"k": -1.0})]}, "must be at least 0"),
        ({"scenarios": [dataclasses.replace(BASE, inputs={"u": Signal.constant(1.0)})]}, "unknown input 'u'"),
        ({"scenarios": [dataclasses.replace(BASE, solver_overrides={"tol": 1})]}, "unknown solver fields"),
        ({"scenarios": [dataclasses.replace(BASE, t_span=(5, 1))]}, "t_span must be increasing"),
        ({"scenarios": [dataclasses.replace(BASE, n_points=10**6)]}, "n_points"),
        ({"references": [Reference("r1", "A paper.")]}, "equations_reference must name"),
        ({"equations_reference": "r1"}, "isn't in references"),
        ({"conserved_quantities": [template.ConservedQuantity("x", ("A", "Z"))]}, "unknown states"),
        ({"auxiliaries": [Auxiliary("N", "A + B")]}, "appears in no equation"),
    ],
)
def test_invalid_entries_are_rejected(overrides, message):
    with pytest.raises(ValueError, match=message):
        make(**overrides)


def test_rhs_must_match_declared_inputs():
    psp = InputVariable("u", "u", "forcing", None, Signal.constant(1.0), (0, 2), "heuristic")
    with pytest.raises(ValueError, match="rhs arguments"):
        make(inputs=[psp])


# -- reference data --------------------------------------------------------------


def dataset(**overrides):
    fields = dict(
        id="measured", path="measured.csv", description="Synthetic data.",
        purpose=frozenset({DataPurpose.DISPLAY}), time_column="time", time_unit="s",
        columns={"a": "A"}, column_units={"a": "mol"}, scenario="base_case",
    )
    return ReferenceDataset(**{**fields, **overrides})


@pytest.fixture
def model_folder(tmp_path):
    (tmp_path / "reference_data").mkdir()
    (tmp_path / "reference_data" / "measured.csv").write_text("time,a\n0,1\n1,0.6\n")
    return tmp_path


def test_valid_reference_data_is_accepted_and_exported(model_folder):
    entry = make(folder=model_folder, reference_data=[dataset()])
    exported = json.loads(json.dumps(entry.to_dict()))["reference_data"][0]
    assert exported["columns"] == {"a": "A"}
    assert exported["purpose"] == ["display"]
    assert entry.summary_dict()["has_reference_data"]


@pytest.mark.parametrize(
    "overrides, message",
    [
        ({"columns": {"a": "Z"}, "column_units": {"a": None}}, "unknown state 'Z'"),
        ({"columns": {"b": "A"}, "column_units": {"b": None}}, r"columns \['b'\] aren't in"),
        ({"column_units": {}}, "has no unit entry"),
        ({"path": "missing.csv"}, "doesn't exist"),
        ({"scenario": "typo"}, "unknown scenario"),
        ({"purpose": frozenset()}, "purpose"),
    ],
)
def test_invalid_reference_data_is_rejected(model_folder, overrides, message):
    with pytest.raises(ValueError, match=message):
        make(folder=model_folder, reference_data=[dataset(**overrides)])


def test_digitized_curve_folder_is_accepted(model_folder):
    curve = model_folder / "reference_data" / "fig2"
    curve.mkdir()
    (curve / "metadata.json").write_text("{}")
    (curve / "run1.csv").write_text("x,y\n0,1\n")
    entry = make(
        folder=model_folder,
        reference_data=[dataset(path="fig2", time_column="x", columns={"y": "A"}, column_units={"y": "mol"})],
    )
    assert np.isfinite(entry.request().t_eval()).all()
