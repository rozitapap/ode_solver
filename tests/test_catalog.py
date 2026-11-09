"""
Every catalog entry must be complete and every scenario runnable, and the
engine must stay independent of the catalog.
"""

import ast
import dataclasses
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import pytest

import ode_catalog
import ode_solver
from ode_catalog import CATALOG, ModelEntry, ModelStatus, Parameter, Scenario
from ode_catalog._template.model import ENTRY as TEMPLATE
from ode_catalog.index import INDEX_PATH, render_index
from ode_catalog.validation import assert_agrees, rk4_fixed_step
from ode_solver import ODESystemSpec, solve_ode, validate_request

ENTRIES = {**CATALOG, "template": TEMPLATE}
IDS = list(ENTRIES)
SCENARIOS = [(entry, scenario.name) for entry in ENTRIES.values() for scenario in entry.scenarios]
SCENARIO_IDS = [f"{entry.spec.name}-{name}" for entry, name in SCENARIOS]


def test_engine_never_imports_the_catalog():
    engine = Path(ode_solver.__file__).parent
    for path in engine.glob("*.py"):
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if isinstance(node, ast.Import):
                modules = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom) and node.level == 0:
                modules = [node.module]
            else:
                continue
            for module in modules:
                assert module.split(".")[0] != "ode_catalog", f"{path.name}: {module}"


@pytest.mark.parametrize("key", list(CATALOG))
def test_registry_key_matches_spec_name(key):
    assert CATALOG[key].spec.name == key


@pytest.mark.parametrize("entry", ENTRIES.values(), ids=IDS)
def test_parameters_are_documented(entry):
    for parameter in entry.parameters:
        assert parameter.description, parameter.name
        if entry.reference is not None:
            assert parameter.source, parameter.name


@pytest.mark.parametrize("entry", ENTRIES.values(), ids=IDS)
def test_default_scenario_matches_the_spec_defaults(entry):
    expected = entry.spec.default_request()
    actual = entry.request()
    assert actual.params == expected.params
    assert actual.initial_conditions == expected.initial_conditions
    assert actual.t_span == expected.t_span
    assert actual.n_points == expected.n_points
    assert set(actual.inputs) == set(entry.spec.default_inputs)


@pytest.mark.parametrize("entry, name", SCENARIOS, ids=SCENARIO_IDS)
def test_scenarios_run(entry, name):
    request = entry.request(name)
    assert validate_request(entry.spec, request) == []
    result = solve_ode(entry.spec, request)
    assert result.success, result.message
    assert len(result.t) == request.n_points


@pytest.mark.parametrize("method", ["DOP853", "Radau", "LSODA"])
@pytest.mark.parametrize("entry, name", SCENARIOS, ids=SCENARIO_IDS)
def test_solver_methods_agree(entry, name, method):
    request = entry.request(name)
    expected = solve_ode(entry.spec, request)
    actual = solve_ode(entry.spec, dataclasses.replace(request, method=method))
    assert actual.success, actual.message
    for state in entry.spec.state_names:
        assert_agrees(actual.trajectories[state], expected.trajectories[state], 1e-5, label=state)


def test_rk4_has_fourth_order_accuracy():
    t = [0.0, 1.0]
    exact = 2.718281828459045
    coarse = abs(rk4_fixed_step(lambda t, y: [y[0]], [1.0], t, substeps=10)[-1, 0] - exact)
    fine = abs(rk4_fixed_step(lambda t, y: [y[0]], [1.0], t, substeps=20)[-1, 0] - exact)
    assert coarse / fine == pytest.approx(16, rel=0.05)


def test_assert_agrees_uses_the_relative_error():
    assert_agrees([100.0, 200.0], [100.0, 200.001], 1e-5)
    with pytest.raises(AssertionError, match="relative error"):
        assert_agrees([100.0, 200.0], [100.0, 201.0], 1e-5, label="x")


def test_index_is_up_to_date():
    assert INDEX_PATH.read_text(encoding="utf-8") == render_index(), (
        "Run python3 -m ode_catalog.index"
    )


@pytest.mark.parametrize("entry", CATALOG.values(), ids=list(CATALOG))
def test_readme_status_matches_the_entry(entry):
    readme = Path(ode_catalog.__file__).parent / entry.spec.name / "README.md"
    assert f"**Status:** {entry.status.value}" in readme.read_text(encoding="utf-8")


def test_unknown_scenario_lists_the_available_ones():
    with pytest.raises(KeyError, match="base_case"):
        TEMPLATE.request("typo")


def make_entry(**overrides):
    fields = {
        "id": "test",
        "name": "Test",
        "version": "1.0.0",
        "domain": "example",
        "summary": "Test entry.",
        "rhs": TEMPLATE.rhs,
        "time_unit": "s",
        "states": TEMPLATE.states,
        "parameters": TEMPLATE.parameters,
        "scenarios": TEMPLATE.scenarios,
        "default_scenario": "base_case",
    }
    return ModelEntry(**{**fields, **overrides})


def test_valid_entry_is_accepted():
    assert make_entry().status is ModelStatus.EXECUTABLE


def test_default_alias_gives_the_default_scenario():
    assert TEMPLATE.scenario("default").name == TEMPLATE.default_scenario


def test_status_names_match_the_biomodels_terms():
    assert [status.value for status in ModelStatus] == [
        "executable",
        "numerically verified",
        "scientifically reproduced",
    ]


K = TEMPLATE.parameters[0]


@pytest.mark.parametrize(
    "overrides, message",
    [
        ({"parameters": []}, "rhs arguments"),
        ({"parameters": [dataclasses.replace(K, key="kk")]}, "rhs arguments"),
        ({"default_scenario": "typo"}, "isn't a scenario"),
        ({"scenarios": [TEMPLATE.scenarios[0]] * 2}, "duplicate scenario"),
        ({"status": ModelStatus.NUMERICALLY_VERIFIED}, "needs evidence"),
        ({"status": ModelStatus.SCIENTIFICALLY_REPRODUCED, "evidence": "x"}, "needs a reference"),
    ],
)
def test_invalid_entries_are_rejected(overrides, message):
    with pytest.raises(ValueError, match=message):
        make_entry(**overrides)


def test_scenario_overrides_apply():
    scenario = Scenario(
        "fast", "Faster rate.", (0, 1), initial_overrides={"A": 2.0}, param_overrides={"k": 5.0},
    )
    entry = make_entry(scenarios=[scenario], default_scenario="fast")
    request = entry.request()
    assert request.params == {"k": 5.0}
    assert request.initial_conditions == {"A0": 2.0, "B0": 0.0}
    assert entry.spec.default_params == {"k": 5.0}
    assert {p.key: p.default for p in entry.parameters} == {"k": 0.5}
