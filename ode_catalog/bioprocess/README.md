# Batch cell-culture kinetics

**Status:** see `model.py`

**Source:** none is recorded. Add a citation here and in `model.py` if these values came from a publication. State and parameter units aren't recorded either. Hours comes only from the time label on the model's earlier plots.

## Biological context
Tracks viable cells in a batch culture, along with the glucose and glutamine they consume and the lactate and ammonia they produce.

## Equations
```
dVCD/dt = mu * VCD
dGlc/dt = -qGlc * VCD
dGln/dt = -qGln * VCD
dLac/dt = qLac * VCD
dAmm/dt = qAmm * VCD
```

## Assumptions and limitations
- Specific rates are constant, with no substrate limitation, so cells grow exponentially without limit.
- Nothing stops glucose or glutamine going negative if they run out, so check them over long runs or with larger rates.

## Differences from the paper
Not applicable.

## Validation
See `test_validation.py`.
