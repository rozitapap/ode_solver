# Mass-spring-damper system

**Status:** numerically verified

**Source:** V. Krcheva, "Modelling a mass-spring system using a second-order homogeneous linear ordinary differential equation with constant coefficients," Balkan Journal of Applied Mathematics and Informatics (BJAMI), vol. 7, no. 1, pp. 7-18, 2024. [Article page](https://js.ugd.edu.mk/index.php/bjami/en/article/view/5932).

## Context
A body of mass `m` hangs from a spring with constant `k` and is attached to a dashpot with damping coefficient `c`. Its displacement `y(t)` is measured from static equilibrium, positive downwards. Pulled down and released, the body oscillates or creeps back to equilibrium, depending on how much damping there is. The paper uses this to show the four regimes of a second-order homogeneous linear ODE with constant coefficients.

## Equations
Newton's second law with the restoring force `-k y` and damping force `-c y'` gives Eq. (2.2.2), which the paper writes as the first-order system Eq. (3.5):

```
m y'' + c y' + k y = 0

y'  = v
v'  = (-c v - k y) / m
```

Here `v` is the paper's `y2 = y'`. The characteristic equation `λ² + (c/m) λ + k/m = 0` (2.2.3) sets the regime:

| Scenario | c | c² vs 4mk | Regime | Paper |
|---|---|---|---|---|
| `no_damping` | 0 | 0 < 64 | Harmonic oscillation, `ω0 = 4` | Eq. (4.1), Fig. 4 |
| `underdamping` (default) | 2 | 4 < 64 | Decaying oscillation, `ω* = √15` | Eq. (4.3), Fig. 5 |
| `critical_damping` | 8 | 64 = 64 | Double root `λ = -4` | Eq. (4.6), Fig. 6 |
| `overdamping` | 10 | 100 > 64 | Roots `λ = -2, -8` | Eq. (4.9), Fig. 7 |

All four use `m = 1`, `k = 16`, `y(0) = 1`, `y'(0) = 0` on `[0, 20]`, as in Section 3.

## Assumptions and limitations
- The spring and damper are linear, and the spring's mass is neglected.
- Motion is strictly vertical. Gravity only shifts the equilibrium, so it doesn't appear in the equation.
- There's no external forcing.

## Differences from the paper
- **Units:** the paper gives `m` in kg, `c` in Ns/m and `k` in N/m, but doesn't print units for `y`, `y'` or `t`. Seconds follow from its frequencies in Hz, and the SI parameter units imply metres for `y`. The state units are left as not recorded.
- **Printed solutions that don't match the paper's own equations:**
  - **Eq. (4.5)**, underdamping: printed as `y = e^{-t} cos(√15 t)`, with `B = 0`. With `y'(0) = 0`, `B = α/ω* = 1/√15`, so the solution is `e^{-t} (cos √15 t + sin(√15 t)/√15)`. The printed function has `y'(0) = -1`. `test_validation.py` checks both.
  - **Eq. (4.8)**, critical damping: printed as `y = (1+4) e^{-4x}`. With the stated `C1 = 1`, `C2 = 4` and Eq. (2.2.7), it should read `(1 + 4t) e^{-4t}`.
  - Eqs. (4.2), (4.5), (4.8) and (4.11) use `x` for time.
- **Other slips that don't affect the model:** the text calls the roots in (4.7) and (4.10) "complex conjugate", but they're real. The paper names its MATLAB function after the Van der Pol oscillator, but the code is the linear system (3.5). The figure for the critically damped case is unnumbered; this README calls it Fig. 6.
- **Default:** the paper has no single default case, so the default scenario is the underdamped one, which is `c = 2` in `model.py`.

## Validation
`test_validation.py` checks that:
- In each scenario, the solver matches the exact closed-form solution for `y` and `y'` within a relative error of 1e-5.
- The solver agrees with the paper's printed solutions (4.2) and (4.11), and with the corrected (4.5) and (4.8).
- The period `T = 1.57` and frequency `f = 0.64 Hz` reported for the undamped case follow from the parameters.
- Mechanical energy is conserved without damping and decreases with it.

That supports **numerically verified**. **Scientifically reproduced** would also need the paper's Figs. 4-7 digitized, compared, and the overlay reviewed by a person. That hasn't been done.
