"""
Exercises the figure-comparison tools on synthetic data. No real model is
compared here, so these tests aren't evidence for any model's status.
"""

import json

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd
import pytest

from ode_catalog.validation import (
    compare_to_curve,
    digitization_error,
    load_reference_curve,
    save_overlay,
)

T = np.linspace(0, 10, 201)


def true_curve(t):
    return 5 * t * np.exp(-t / 2)


def write_curve(folder, noise=(0.05, 0.05), x_scale=1.0, y_scale=1.0, **metadata_overrides):
    folder.mkdir()
    rng = np.random.default_rng(0)
    x = np.linspace(0, 10, 41)
    for i, scale in enumerate(noise, start=1):
        y = true_curve(x) + rng.normal(0, scale, len(x))
        pd.DataFrame({"x": x / x_scale, "y": y / y_scale}).to_csv(folder / f"run{i}.csv", index=False)
    metadata = {
        "paper": "Synthetic",
        "figure": "Fig. 0",
        "scenario": "base_case",
        "state": "A",
        "x_unit": "s",
        "y_unit": "mol",
        "x_scale": x_scale,
        "y_scale": y_scale,
        "tool": "None",
        "digitized_by": ["Test"],
        "overlay_reviewed_by": None,
        **metadata_overrides,
    }
    (folder / "metadata.json").write_text(json.dumps(metadata), encoding="utf-8")
    return folder


def test_matching_simulation_passes(tmp_path):
    curve = load_reference_curve(write_curve(tmp_path / "curve"))
    assert compare_to_curve(T, true_curve(T), curve).passed


@pytest.mark.parametrize(
    "wrong, field",
    [
        (lambda t: true_curve(t) * 1.2, "curve_error"),
        (lambda t: true_curve(t) * 1.2, "peak_height_error"),
        (lambda t: true_curve(t - 1.5), "peak_time_error"),
        (lambda t: true_curve(t) + 0.5 * (t > 9), "final_value_error"),
    ],
)
def test_wrong_simulation_fails_on_the_right_measure(tmp_path, wrong, field):
    curve = load_reference_curve(write_curve(tmp_path / "curve"))
    comparison = compare_to_curve(T, wrong(T), curve)
    assert not comparison.passed
    assert getattr(comparison, field) > comparison.tolerance


def test_noisier_digitizing_gives_a_wider_tolerance(tmp_path):
    precise = load_reference_curve(write_curve(tmp_path / "precise", noise=(0.05, 0.05)))
    rough = load_reference_curve(write_curve(tmp_path / "rough", noise=(0.05, 0.4)))
    assert digitization_error(rough) > digitization_error(precise)
    assert (
        compare_to_curve(T, true_curve(T), rough).tolerance
        > compare_to_curve(T, true_curve(T), precise).tolerance
    )


def test_identical_runs_still_get_the_minimum_tolerance(tmp_path):
    curve = load_reference_curve(write_curve(tmp_path / "curve", noise=(0.0, 0.0)))
    assert digitization_error(curve) == 0
    assert compare_to_curve(T, true_curve(T), curve, minimum_tolerance=0.01).tolerance == 0.01


def test_scales_convert_figure_units_to_model_units(tmp_path):
    curve = load_reference_curve(write_curve(tmp_path / "curve", x_scale=60.0, y_scale=1 / 60))
    assert curve.runs[0][-1, 0] == pytest.approx(10.0)
    assert compare_to_curve(T, true_curve(T), curve).passed


def test_one_digitizing_run_is_rejected(tmp_path):
    folder = write_curve(tmp_path / "curve")
    (folder / "run2.csv").unlink()
    with pytest.raises(ValueError, match="at least two"):
        load_reference_curve(folder)


def test_missing_metadata_is_rejected(tmp_path):
    folder = write_curve(tmp_path / "curve")
    metadata = json.loads((folder / "metadata.json").read_text(encoding="utf-8"))
    del metadata["figure"]
    (folder / "metadata.json").write_text(json.dumps(metadata), encoding="utf-8")
    with pytest.raises(ValueError, match="figure"):
        load_reference_curve(folder)


def test_curve_is_unreviewed_until_someone_records_a_review(tmp_path):
    assert load_reference_curve(write_curve(tmp_path / "a")).reviewed_by is None
    reviewed = write_curve(tmp_path / "b", overlay_reviewed_by="A Reviewer")
    assert load_reference_curve(reviewed).reviewed_by == "A Reviewer"


def test_overlay_is_saved(tmp_path):
    curve = load_reference_curve(write_curve(tmp_path / "curve"))
    path = tmp_path / "out" / "overlay.png"
    save_overlay(T, true_curve(T), curve, path)
    assert path.exists()
