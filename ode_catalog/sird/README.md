# SIRD epidemic model

**Status:** numerically verified

**Source:** not from a publication. The parameter values are loosely representative of an early-outbreak respiratory illness, and aren't fitted to data.

## Biological context
People move from susceptible to infected, then either recover or die. It's a standard way to illustrate how an outbreak grows and declines.

## Equations
```
N = S + I + R + D
dS/dt = -beta * S * I / N
dI/dt = beta * S * I / N - gamma * I - mu * I
dR/dt = gamma * I
dD/dt = mu * I
```

## Assumptions and limitations
- The population is closed and well mixed, with no births, other deaths, or loss of immunity.
- Rates are constant over time.
- Because no source is recorded, the results shouldn't be read as a prediction for any real disease.

## Differences from the paper
Not applicable.

## Validation
**Numerically verified.** Within a relative error of 1e-5, the solver matches the exact relations `S = S0 * exp(-(beta/N)/(gamma+mu) * removed)` and `D/R = mu/gamma`, conserves the total population, and agrees with a separate fixed-step RK4 implementation. Different solver methods agree within the same tolerance. This says the equations are solved correctly, not that they describe a real outbreak. See `test_validation.py`.
