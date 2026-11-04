# Adding a model

1. Copy `_template/` to `ode_catalog/<model_name>/`.
2. Replace the equations, parameters, scenarios and README with your model's.
3. Add it to `CATALOG` in `ode_catalog/__init__.py`.
4. Run `python3 -m ode_catalog.index` to update the catalog index.
5. Run `python3 -m pytest`. The contract tests pick up the new model automatically.

## Checklist

**Source**
- [ ] `reference` has the citation, or is `None` if the model isn't from a publication.
- [ ] Citation details come from the paper itself, not from memory.
- [ ] The paper's PDF and images are not in the repo. Link to them instead.

**Equations**
- [ ] Each equation has been checked against the paper by someone other than the person who entered it, if possible.
- [ ] `test_validation.py` has an independent check of the equations.

**Parameters**
- [ ] Every parameter has a value, unit, description and `source`.
- [ ] Conversions, and anything that looks wrong in the paper, are recorded in `note`. Values and units are never silently corrected.

**Scenarios**
- [ ] Every scenario has a `source`, or `None` if it isn't from the reference.
- [ ] Solver settings that affect the result are in the scenario, not in a script.

**README**
- [ ] It explains the biological context, assumptions, and any differences from the paper.

## Choosing a status

Each status needs evidence in `test_validation.py` or the contract tests. Start at `EXECUTABLE`, and summarize what supports any higher status in the entry's `evidence` field. A conservation check alone doesn't meet the definition of numerically verified, because it isn't an independent implementation or a trusted reference.

| Status | Definition | Evidence |
|---|---|---|
| `EXECUTABLE` | Simulation runs without crashing | The contract tests pass |
| `NUMERICALLY_VERIFIED` | Results agree across independent implementations or against a trusted reference within defined tolerances | A test comparing against the other implementation or reference, with the tolerance stated |
| `SCIENTIFICALLY_REPRODUCED` | Simulation reproduces a documented, published result | A test comparing against that result, and a person's review of the overlay |

These definitions come from BioModels' empirical quality control findings. Status records validation only. Whether a model comes from a publication is recorded separately by `reference`.

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
4. In `test_validation.py`, use `load_reference_curve` and `compare_to_curve` from `ode_catalog.validation`. The comparison covers the overall curve, peak height, peak time and final value. Its tolerance comes from how much your two digitizing runs differ.
5. Generate an overlay with `save_overlay`. A person then checks that it's the right figure and conditions, and records their name in `overlay_reviewed_by`.
6. Only then set the status to `SCIENTIFICALLY_REPRODUCED`, and summarize the comparison in `evidence`. The test should assert both `comparison.passed` and `curve.reviewed_by`.

If the comparison fails, the cause could be the paper, the digitizing or the code, so check all three before changing anything.

