# Adding a model

1. Copy `_template/` to `ode_catalog/<model_name>/`.
2. Replace the equations, parameters, scenarios and README with your model's.
3. Add it to `CATALOG` in `ode_catalog/__init__.py`.
4. Run `python3 -m pytest`. The contract tests pick up the new model automatically.

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
