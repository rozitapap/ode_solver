# SIRD epidemic model

**Status:** see `model.py`

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
See `test_validation.py`.
