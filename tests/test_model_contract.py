"""
Every registered model must run through the whole pipeline, and the pipeline
must handle systems unlike the registered ones. A new model is covered
automatically once it's in MODELS.
"""

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest

from data_generator import generate_dataset
from models import MODELS
from models.base import ODESystemSpec
from ODE_solver import plot_results, solve_ode, validate_request

SPECS = list(MODELS.values())
IDS = list(MODELS)


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
    assert len(spec.rhs(request.t_span[0], y0, **request.params)) == len(spec.state_names)


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


def test_single_state_system_solves_and_sweeps():
    result = solve_ode(DECAY, DECAY.default_request())
    assert result.trajectories["x"][-1] == pytest.approx(np.exp(-1.0), rel=1e-6)
    assert generate_dataset(DECAY, DECAY.default_request(), n_runs=2, seed=0)["success"].all()


@pytest.mark.xfail(
    raises=TypeError,
    strict=True,
    reason="plot_results assumes at least 2 states (plt.subplots returns a bare Axes)",
)
def test_single_state_system_plots(tmp_path):
    plot_results(solve_ode(DECAY, DECAY.default_request()), DECAY, save_path=tmp_path / "plot.png")


def test_failed_integration_is_reported():
    assert not solve_ode(BLOWUP, BLOWUP.default_request()).success


def test_failed_runs_are_flagged_in_the_dataset():
    with pytest.warns(UserWarning, match="2 of 2 runs failed"):
        df = generate_dataset(BLOWUP, BLOWUP.default_request(), n_runs=2, seed=0)
    assert not df["success"].any()
    assert df["run"].nunique() == 2


def test_run_with_no_points_still_appears_in_the_dataset():
    def nan_rhs(t, y, k):
        return [np.nan]

    spec = ODESystemSpec("nan", ["x"], ["k"], nan_rhs)
    request = DECAY.default_request()
    with pytest.warns(UserWarning):
        df = generate_dataset(spec, request, n_runs=2, seed=0)
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


def test_param_names_order_does_not_change_the_result():
    sird = MODELS["sird"]
    reordered = ODESystemSpec("reordered", sird.state_names, ["mu", "gamma", "beta"], sird.rhs)
    expected = solve_ode(sird, sird.default_request())
    actual = solve_ode(reordered, sird.default_request())
    np.testing.assert_array_equal(actual.trajectories["I"], expected.trajectories["I"])
