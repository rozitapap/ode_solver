"""The experiment scripts must run for any catalog model, and write where they're told."""

import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from ode_catalog import CATALOG

ROOT = Path(__file__).resolve().parent.parent
ENV = {**os.environ, "PYTHONPATH": str(ROOT), "MPLBACKEND": "Agg", "PYTHONDONTWRITEBYTECODE": "1"}


def run_script(name, *args):
    return subprocess.run(
        [sys.executable, str(ROOT / "experiments" / name), *args],
        capture_output=True,
        text=True,
        env=ENV,
        check=False,
    )


@pytest.mark.parametrize("model", CATALOG)
def test_perturbation_study_runs_for_every_model(model, tmp_path):
    done = run_script(
        "perturbation_study.py", "--model", model, "--n-runs", "2", "--output-dir", str(tmp_path)
    )
    assert done.returncode == 0, done.stderr
    assert (tmp_path / "trajectories.png").exists()
    df = pd.read_csv(tmp_path / "trajectories.csv")
    assert df["run"].nunique() == 2
    assert df["success"].all()


@pytest.mark.parametrize("model", CATALOG)
def test_single_run_runs_for_every_model(model, tmp_path):
    done = run_script("single_run.py", "--model", model, "--output-dir", str(tmp_path))
    assert done.returncode == 0, done.stderr
    assert (tmp_path / f"{CATALOG[model].default_scenario}.png").exists()
    assert "run complete" in done.stdout


def test_single_run_accepts_a_scenario(tmp_path):
    done = run_script(
        "single_run.py", "--model", "brain_heart", "--scenario", "setpoint_step",
        "--output-dir", str(tmp_path),
    )
    assert done.returncode == 0, done.stderr
    assert (tmp_path / "setpoint_step.png").exists()


def test_unknown_model_is_rejected(tmp_path):
    done = run_script("single_run.py", "--model", "typo", "--output-dir", str(tmp_path))
    assert done.returncode != 0
    assert "invalid choice" in done.stderr
