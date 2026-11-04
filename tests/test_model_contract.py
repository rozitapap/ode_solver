"""
Every registered model must run through the whole pipeline, and the pipeline
must handle systems unlike the registered ones. A new model is covered
automatically once it's in CATALOG.
"""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from ode_catalog import CATALOG
from ode_catalog.brain_heart.model import psp_constant, psp_profile
from ode_solver import (
    ODESystemSpec,
    Pipeline,
    SimulationRequest,
    generate_dataset,
    plot_results,
    plot_sweep,
    solve_ode,
    validate_request,
)

ALL_SPECS = {name: entry.spec for name, entry in CATALOG.items()}
SPECS = list(ALL_SPECS.values())
IDS = list(ALL_SPECS)


def decay_rhs(t, y, k):
    return [-k * y[0]]


# One state, which also has an unconventional name.
DECAY = ODESystemSpec(
    name="decay",
    state_names=["x"],
    param_names=["k"],
    rhs=decay_rhs,
    default_params={"k": 0.1},
    default_initial_conditions={"x0": 1.0},
    default_t_span=(0, 10),
    default_n_points=11,
)


def blowup_rhs(t, y, k):
    return [k * y[0] ** 2]


# Diverges at t = 1, so the integration can't reach the end of the time span.
BLOWUP = ODESystemSpec(
    name="blowup",
    state_names=["x"],
    param_names=["k"],
    rhs=blowup_rhs,
    default_params={"k": 1.0},
    default_initial_conditions={"x0": 1.0},
    default_t_span=(0, 2),
    default_n_points=21,
)


def glucose_insulin_rhs(t, y, k):
    S, I, _ = y
    return [-k * S * I, 0.1 * S - 0.2 * I, 0.0]


# Its state names collide with epidemic compartments, but it isn't an epidemic.
LOOKS_EPIDEMIC = ODESystemSpec(
    name="looks_epidemic",
    state_names=["S", "I", "D"],
    param_names=["k"],
    rhs=glucose_insulin_rhs,
    default_params={"k": 0.01},
    default_initial_conditions={"S0": 100.0, "I0": 10.0, "D0": 1.0},
    default_t_span=(0, 50),
)


@pytest.mark.parametrize("key, spec", zip(IDS, SPECS), ids=IDS)
def test_registry_key_matches_spec_name(key, spec):
    assert key == spec.name


@pytest.mark.parametrize("spec", SPECS, ids=IDS)
def test_default_request_is_valid(spec):
    assert validate_request(spec, spec.default_request()) == []


@pytest.mark.parametrize("spec", SPECS, ids=IDS)
def test_rhs_returns_one_derivative_per_state(spec):
    request = spec.default_request()
    y0 = [request.initial_conditions[name + "0"] for name in spec.state_names]
    derivatives = spec.rhs(request.t_span[0], y0, **request.params, **request.inputs)
    assert len(derivatives) == len(spec.state_names)


@pytest.mark.parametrize("spec", SPECS, ids=IDS)
def test_solves_over_the_full_time_grid(spec):
    request = spec.default_request()
    result = solve_ode(spec, request)
    assert result.success, result.message
    assert len(result.t) == request.n_points
    assert list(result.trajectories) == list(spec.state_names)
    for series in result.trajectories.values():
        assert np.isfinite(series).all()


@pytest.mark.parametrize("spec", SPECS, ids=IDS)
def test_conserved_total_holds_when_declared(spec):
    if not spec.conserved_total:
        pytest.skip("system does not declare a conserved total")
    result = solve_ode(spec, spec.default_request())
    totals = sum(result.trajectories.values())
    assert totals[-1] == pytest.approx(totals[0], rel=1e-6)


@pytest.mark.parametrize("spec", SPECS, ids=IDS)
def test_sweep_produces_flagged_rows_for_every_run(spec):
    request = spec.default_request()
    df = generate_dataset(spec, request, n_runs=2, seed=0)
    assert len(df) == 2 * request.n_points
    assert df["success"].all()
    assert set(spec.state_names) | set(spec.param_names) <= set(df.columns)


@pytest.mark.parametrize("spec", SPECS, ids=IDS)
def test_plots_and_saves(spec, tmp_path):
    path = tmp_path / "plot.png"
    plot_results(solve_ode(spec, spec.default_request()), spec, save_path=path)
    assert path.exists()


@pytest.mark.parametrize("spec", SPECS, ids=IDS)
def test_pipeline_runs_end_to_end(spec, tmp_path):
    plot_path = tmp_path / "plot.png"
    run = Pipeline(verbose=False).run(spec, spec.default_request(), plot_path=str(plot_path))
    assert run.ok, run.logs
    assert plot_path.exists()
    for state in spec.state_names:
        assert state in run.analysis
        assert any(f"- {state}:" in line for log in run.logs for line in log.splitlines())


@pytest.mark.parametrize("spec", SPECS, ids=IDS)
def test_drift_is_reported_only_when_declared(spec):
    run = Pipeline(verbose=False).run(spec, spec.default_request())
    assert ("_conservation" in run.analysis) == spec.conserved_total
    assert any("drift" in line for line in run.logs) == spec.conserved_total


def test_state_names_do_not_trigger_domain_language():
    run = Pipeline(verbose=False).run(LOOKS_EPIDEMIC, LOOKS_EPIDEMIC.default_request())
    logs = "\n".join(run.logs).lower()
    for word in ("infection", "susceptible", "deaths", "population"):
        assert word not in logs


def test_single_state_system_solves_and_sweeps():
    result = solve_ode(DECAY, DECAY.default_request())
    assert result.trajectories["x"][-1] == pytest.approx(np.exp(-1.0), rel=1e-6)
    assert generate_dataset(DECAY, DECAY.default_request(), n_runs=2, seed=0)["success"].all()


def test_single_state_system_plots(tmp_path):
    path = tmp_path / "plot.png"
    plot_results(solve_ode(DECAY, DECAY.default_request()), DECAY, save_path=path)
    assert path.exists()


def test_single_state_sweep_plots(tmp_path):
    path = tmp_path / "plot.png"
    df = generate_dataset(DECAY, DECAY.default_request(), n_runs=2, seed=0)
    plot_sweep(df, DECAY, save_path=path)
    assert path.exists()


@pytest.fixture
def open_figure(monkeypatch):
    """Keep the figure open after plotting, so its axes can be inspected."""
    monkeypatch.setattr(plt, "close", lambda *args: None)
    return plt.gcf


def test_axes_use_the_given_units_and_nothing_else(open_figure):
    entry = CATALOG["sird"]
    plot_results(solve_ode(entry.spec, entry.request()), entry.spec, **entry.units)
    axes = open_figure().axes
    assert [ax.get_ylabel() for ax in axes] == [f"{name} (people)" for name in "SIRD"]
    assert axes[-1].get_xlabel() == "Time (day)"


def test_axes_leave_out_units_that_arent_recorded(open_figure):
    plot_results(solve_ode(DECAY, DECAY.default_request()), DECAY)
    ax = open_figure().axes[0]
    assert ax.get_ylabel() == "x"
    assert ax.get_xlabel() == "Time"


def test_negative_values_are_not_clipped(open_figure):
    spec = ODESystemSpec(
        name="oscillator",
        state_names=["x", "v"],
        param_names=["w"],
        rhs=lambda t, y, w: [y[1], -w * w * y[0]],
        default_params={"w": 1.0},
        default_initial_conditions={"x0": 1.0, "v0": 0.0},
        default_t_span=(0, 10),
    )
    plot_results(solve_ode(spec, spec.default_request()), spec)
    assert open_figure().axes[0].get_ylim()[0] < -0.9


def test_failed_integration_is_reported():
    assert not solve_ode(BLOWUP, BLOWUP.default_request()).success
    run = Pipeline(verbose=False).run(BLOWUP, BLOWUP.default_request())
    assert not run.ok
    assert run.analysis is None


def test_failed_runs_are_flagged_in_the_dataset():
    with pytest.warns(UserWarning, match="2 of 2 runs failed"):
        df = generate_dataset(BLOWUP, BLOWUP.default_request(), n_runs=2, seed=0)
    assert not df["success"].any()
    assert df["run"].nunique() == 2


def test_run_with_no_points_still_appears_in_the_dataset():
    def nan_rhs(t, y, k):
        return [np.nan]

    spec = ODESystemSpec("nan", ["x"], ["k"], nan_rhs)
    with pytest.warns(UserWarning):
        df = generate_dataset(spec, DECAY.default_request(), n_runs=2, seed=0)
    assert df["run"].nunique() == 2
    assert not df["success"].any()


def test_non_finite_trajectory_is_a_failure():
    def inf_rhs(t, y, k):
        return [np.inf]

    spec = ODESystemSpec("inf", ["x"], ["k"], inf_rhs)
    assert not solve_ode(spec, DECAY.default_request()).success


def test_derivative_going_nan_mid_run_is_a_failure():
    def late_nan_rhs(t, y, k):
        return [np.nan if t > 5 else -k * y[0]]

    spec = ODESystemSpec("late_nan", ["x"], ["k"], late_nan_rhs)
    result = solve_ode(spec, DECAY.default_request())
    assert not result.success
    assert 0 < len(result.t) < DECAY.default_n_points


def test_unknown_parameter_is_rejected():
    request = DECAY.default_request()
    request.params["typo"] = 1.0
    assert any("typo" in problem for problem in validate_request(DECAY, request))
    with pytest.raises(ValueError, match="typo"):
        generate_dataset(DECAY, request, n_runs=1)
    assert not Pipeline(verbose=False).run(DECAY, request).ok


def test_param_names_order_does_not_change_the_result():
    sird = ALL_SPECS["sird"]
    reordered = ODESystemSpec("reordered", sird.state_names, ["mu", "gamma", "beta"], sird.rhs)
    expected = solve_ode(sird, sird.default_request())
    actual = solve_ode(reordered, sird.default_request())
    np.testing.assert_array_equal(actual.trajectories["I"], expected.trajectories["I"])


def brain_heart_request(**kwargs):
    spec = ALL_SPECS["brain_heart"]
    return SimulationRequest(
        system_name=spec.name,
        params=dict(spec.default_params),
        initial_conditions=dict(spec.default_initial_conditions),
        t_span=spec.default_t_span,
        **kwargs,
    )


def test_missing_input_falls_back_to_the_default_and_says_so():
    spec = ALL_SPECS["brain_heart"]
    run = Pipeline(verbose=False).run(spec, brain_heart_request())
    assert run.ok, run.logs
    assert run.result.metadata["defaulted_inputs"] == ["Psp"]
    assert any("'Psp' not provided, using default" in line for line in run.logs)

    explicit = solve_ode(spec, brain_heart_request(inputs={"Psp": psp_constant}))
    np.testing.assert_array_equal(run.result.trajectories["P"], explicit.trajectories["P"])


def test_provided_input_overrides_the_default():
    spec = ALL_SPECS["brain_heart"]
    run = Pipeline(verbose=False).run(spec, brain_heart_request(inputs={"Psp": psp_profile}))
    default = solve_ode(spec, brain_heart_request())
    assert run.result.metadata["defaulted_inputs"] == []
    assert not any("using default" in line for line in run.logs)
    assert run.result.trajectories["P"][-1] > default.trajectories["P"][-1]


def test_missing_input_without_a_default_fails_validation():
    spec = ODESystemSpec(
        name="needs_input",
        state_names=["x"],
        param_names=["k"],
        rhs=lambda t, y, k, u: [-k * y[0] + u(t)],
        input_names=["u"],
        default_params={"k": 1.0},
        default_initial_conditions={"x0": 1.0},
        default_t_span=(0, 1),
    )
    run = Pipeline(verbose=False).run(spec, spec.default_request())
    assert not run.ok
    assert any("Missing time-varying inputs" in line for line in run.logs)
