# Adding a model

Every model has the same structure, defined by `ModelEntry` in `schema.py`, so the experiments and the website never need to special-case one. An entry is checked by `validate_entry()` when it's created, so an incomplete entry fails at import with a list of everything wrong.

## Steps

1. Copy `_template/` to `ode_catalog/<model_id>/`. The folder name is the model's `id`: lowercase letters, digits and underscores. It's used in URLs and the API, so never change it once published.
2. In `model.py`, replace every value marked `REPLACE`:
   1. Write `rhs(t, y, <parameters...>, <inputs...>)`. It returns derivatives in the same order as `STATES`, and it's the only place the maths lives.
   2. Declare one `StateVariable` per state, in `rhs` order, each with a LaTeX `equation`.
   3. Declare one `Parameter` per `rhs` argument, and one `InputVariable` per time-varying input.
   4. Declare any `Auxiliary` symbols the equations use, such as `N = S + I + R + D`.
   5. Add `Reference`s and set `equations_reference`, or leave both empty if the model isn't from a publication.
   6. Add at least one `Scenario` and set `default_scenario`.
   7. List anything you couldn't verify in `todo`.
3. Write `README.md`. It's the model's long description, both on GitHub and on the website.
4. Write `test_validation.py` with whatever evidence supports the status you choose.
5. Add the entry to `_ENTRIES` in `ode_catalog/__init__.py`.
6. Run `python3 -m ode_catalog.index` to update `INDEX.md`.
7. Run `python3 -m pytest`. The contract, schema and experiment tests pick up the new model automatically.
8. Run `python3 -m ode_catalog.export` if you need a fresh `catalog.json` for the website.

## Required fields

**Identity** (`ModelEntry`)
- [ ] `id`: matches the folder name and `^[a-z][a-z0-9_]*$`.
- [ ] `name`: display name.
- [ ] `version`: `MAJOR.MINOR.PATCH`. Bump the major version when the equations or default values change.
- [ ] `domain`: e.g. `epidemiology`, `bioprocess`, `physiology`, `biochemistry`. Add `tags` if they help.
- [ ] `status`, with `evidence` for anything above `EXECUTABLE`. See below.

**Description**
- [ ] `summary`: one line.
- [ ] `README.md`: context, equations, assumptions and limitations, differences from the paper, and validation. It must contain `**Status:** <status>`.

**Citation**
- [ ] `references`: one `Reference` per publication. `citation` is the full citation. The structured fields (`authors`, `title`, `year`, `venue`, `doi`, `url`) are optional: fill them only from the paper itself, never from memory.
- [ ] `equations_reference`: the `id` of the reference the equations come from, or `None`.
- [ ] The paper's PDF and images are not in the repo. Link to them instead.

**States** (`StateVariable`, one per state, in `rhs` order)
- [ ] `key` (used in scenarios and the API), LaTeX `symbol`, `description`, and `unit` (`None` if it isn't recorded).
- [ ] `initial`: the default initial value.
- [ ] `lower`/`upper`: hard bounds, enforced on every request. For example, `lower=0` for amounts.
- [ ] `ui_range` and `ui_range_basis`: the slider range and why it was chosen. Write `"heuristic: ..."` if it isn't from a source.
- [ ] `equation`: the LaTeX right-hand side of `d(key)/dt`.
- [ ] `conserved_quantities` on the entry, if any set of states sums to a constant. When one covers every state, the engine reports conservation drift.

**Parameters** (`Parameter`, one per `rhs` argument)
- [ ] `key` (the `rhs` argument name), LaTeX `symbol`, `description`, `default`.
- [ ] `unit` as the source prints it, or `None`. Conversions and anything that looks wrong go in `note`. Values and units are never silently corrected.
- [ ] `lower`/`upper` hard bounds. Use `lower_exclusive=True` for strictly positive values.
- [ ] `ui_range`, `ui_range_basis` and `scale` (`Scale.LOG` for rates spanning orders of magnitude; the range must then be positive).
- [ ] `source` (e.g. `"Table I"`) and `reference` (a `Reference.id`) if the value comes from a publication.
- [ ] Rules that involve several values, like "gamma + mu > 0", go in `extra_checks`.

**Inputs** (`InputVariable`, one per time-varying `rhs` argument)
- [ ] `default` is a `Signal`: `Signal.constant(v)`, `Signal.step(before, after, at)`, or `Signal("piecewise_constant", values, times)`. Don't pass Python functions: a website can't show or edit them.

**Equations**
- [ ] Every state has a LaTeX `equation` that uses only declared symbols, `t`, numbers, and `\frac`, `\left`, `\right`, `\cdot`, `\times`, `\exp`, `\ln`, `\log`, `\sqrt`. Every parameter, input and auxiliary symbol appears in at least one equation. This check catches renamed or undeclared symbols, but not wrong maths.
- [ ] If possible, someone other than the person who entered the equations has checked them against the paper.
- [ ] `test_validation.py` has an independent check of the equations.

**Scenarios** (`Scenario`)
- [ ] At least one, with `default_scenario` naming it. `entry.scenario("default")` also returns it.
- [ ] Each has a `description`, `t_span` and `n_points`, plus `initial_overrides` (by state key), `param_overrides`, `inputs` and `solver_overrides` for anything that differs from the defaults.
- [ ] Each has a `source`, or `None` if it isn't from the reference.
- [ ] Solver settings that affect the result are in the scenario or the entry's `solver`, not in a script.

**Solver hints** (`SolverHints`, optional)
- [ ] `method`, `rtol`, `atol`, `max_step`, `stiff`, and `max_points` (the largest `n_points` a backend should accept).

**Reference data** (`ReferenceDataset`, optional)
- [ ] The file or folder is under `reference_data/` and exists.
- [ ] `columns` maps each data column to a state, with a unit for each in `column_units`. `time_column` and `time_unit` describe time.
- [ ] `purpose` is `{DataPurpose.VALIDATION}`, `{DataPurpose.DISPLAY}`, or both.
- [ ] `reference` and `scenario` name what the data came from and what it should be compared with.

## Choosing a status

Each status needs evidence in `test_validation.py` or the contract tests. Start at `EXECUTABLE`, and summarize what supports any higher status in the entry's `evidence` field. A conservation check alone doesn't meet the definition of numerically verified, because it isn't an independent implementation or a trusted reference.

| Status | Definition | Evidence |
|---|---|---|
| `EXECUTABLE` | Simulation runs without crashing | The contract tests pass |
| `NUMERICALLY_VERIFIED` | Results agree across independent implementations or against a trusted reference within defined tolerances | A test comparing against the other implementation or reference, with the tolerance stated |
| `SCIENTIFICALLY_REPRODUCED` | Simulation reproduces a documented, published result | A test comparing against that result, and a person's review of the overlay |

These definitions come from BioModels' empirical quality control findings. Status records validation only. Whether a model comes from a publication is recorded separately by `equations_reference`.

## Comparing against a published figure

Scientifically reproduced needs a documented result from a forward simulation with known inputs. An estimate fitted to measured data isn't a suitable target, because the simulation can't be run under the same conditions.

1. Digitize the curve with a tool like WebPlotDigitizer, then do it again from scratch. Two people is better.
2. Save each attempt in `reference_data/<curve_name>/` as `run1.csv` and `run2.csv`, with columns `x` and `y` in the figure's units.
3. Add `metadata.json` to the same folder:

   ```json
   {
     "paper": "Short citation",
     "figure": "Fig. 2, top panel",
     "scenario": "name_of_the_matching_scenario",
     "state": "P",
     "x_unit": "s",
     "y_unit": "mmHg",
     "x_scale": 1.0,
     "y_scale": 1.0,
     "tool": "WebPlotDigitizer",
     "digitized_by": ["Name", "Name"],
     "overlay_reviewed_by": null
   }
   ```

   `x_scale` and `y_scale` convert the figure's units to the model's, for example `0.016667` from beats/min to beats/s.
4. Declare it in the entry:

   ```python
   ReferenceDataset(
       id="fig2_pressure", path="fig2_pressure", description="Fig. 2, top panel",
       purpose=frozenset({DataPurpose.VALIDATION, DataPurpose.DISPLAY}),
       time_column="x", time_unit="s", columns={"y": "P"}, column_units={"y": "mmHg"},
       scenario="name_of_the_matching_scenario", reference="<reference id>",
   )
   ```
5. In `test_validation.py`, use `load_reference_curve` and `compare_to_curve` from `ode_catalog.validation`. The comparison covers the overall curve, peak height, peak time and final value. Its tolerance comes from how much your two digitizing runs differ.
6. Generate an overlay with `save_overlay`. A person then checks that it's the right figure and conditions, and records their name in `overlay_reviewed_by`.
7. Only then set the status to `SCIENTIFICALLY_REPRODUCED`, and summarize the comparison in `evidence`. The test should assert both `comparison.passed` and `curve.reviewed_by`.

If the comparison fails, the cause could be the paper, the digitizing or the code, so check all three before changing anything.
