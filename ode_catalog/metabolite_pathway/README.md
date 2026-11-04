# Linear metabolite pathway

**Status:** numerically verified

**Source:** not from a publication. The rate constants are illustrative, and no units are recorded for them.

## Biological context
A generic chain of first-order reactions, where each metabolite converts into the next. It doesn't represent a specific pathway.

## Equations
```
dM1/dt = -k1 * M1
dM2/dt = k1 * M1 - k2 * M2
dM3/dt = k2 * M2 - k3 * M3
dM4/dt = k3 * M3 - k4 * M4
dM5/dt = k4 * M4
```

## Assumptions and limitations
- Every step is irreversible and first-order, with no enzyme saturation or regulation.
- The total amount of material is conserved.

## Differences from the paper
Not applicable.

## Validation
**Numerically verified.** The solver matches the exact Bateman solution within a relative error of 1e-5, and different solver methods agree within the same tolerance. See `test_validation.py`.
