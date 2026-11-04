# Autonomic cardiac regulation (delay-free baroreflex model)

**Status:** numerically verified

**Source:** S. N. Sadoun, G. A. D'Inverno, A. Boutin, F. Cottin, and T.-M. Laleg-Kirati, "Physics-Informed Neural Estimation of State and Unknown Input in Autonomic Cardiac Dynamics with Left-Invertibility Constraints," 2026 European Control Conference (ECC), Reykjavík, Iceland, July 7-10, 2026, pp. 1795-1800.

## Biological context
The model describes one part of brain-heart interaction: how the autonomic nervous system regulates heart rate H and mean arterial blood pressure P through the baroreflex. The baroreflex set-point Psp(t) acts as a control input from the brain, shifting the pressure the reflex responds to.

Psp can't be measured directly. The reference treats it as an unknown input and estimates it, along with P, from measured heart rate using a physics-informed neural network.

## Equations
From Eqs. (3) and (4) of the reference, with the parasympathetic and sympathetic activity:

```
Tp = 1 / (1 + exp(-asp * (P - Psp(t))))
Ts = 1 - Tp

dH/dt = betaH * Ts / (1 + gamma * Tp) - VH * Tp + deltaH * (H0 - H)
dP/dt = -P / (Ra * (1 + alpha * Ts) * Ca) + H * dV / Ca
```

The parameter values are from Table I and are recorded with their sources in `model.py`.

## Assumptions and limitations
- This is the delay-free version of Ottesen's 1997 model. The reference notes that the delay is responsible for certain physiological oscillations (Mayer waves), so this model doesn't produce them.
- Parameters are nominal values, not fitted to an individual.

## Differences from the paper
- **Time units:** this implementation uses seconds throughout, so H is in beats/s and `H0 = 100/60` instead of Table I's 100 min^-1.
- **Printed units for `Ra` and `Ca`:** Table I prints `mmHg·g·ml^-1` and `mmHg·g^-1`, but Eq. (4) is only dimensionally consistent with `mmHg·s·ml^-1` and `ml·mmHg^-1`. `model.py` records them as printed, with a note. Check them against Ottesen (1997).
- **`deltaH`:** Table I calls it a relaxation time, but its unit and role are those of a rate.
- **Scenarios aren't from the paper:** the paper estimates Psp from data, so `resting` and `setpoint_step`, including their initial conditions, are our own. In `resting`, P settles from 90 to about 85.7 mmHg, so the first 10 s or so are a transient.

## Validation
`test_validation.py` checks that:
- `rhs` agrees with a separate transcription of Eqs. (3) and (4).
- Simulated trajectories satisfy the paper's algebraic relation, Eq. (13).

- Results agree with a separate fixed-step RK4 implementation of Eqs. (3) and (4), within a relative error of 1e-5, in both scenarios.
- Different solver methods agree within the same tolerance.

That supports **numerically verified**, with one caveat: both implementations were written in this project. They share no integrator and were typed separately, but an error in reading the paper could affect both. An external implementation would be stronger evidence.

**Scientifically reproduced** needs a documented result to compare against, and that isn't straightforward with this paper:
- Fig. 3 shows neural-network estimates from the PASS dataset, not a forward simulation with a known input.
- The paper states there's no ground truth for P and Psp.

Options are to use a forward simulation published in Ottesen (1997) or Ataee et al. (2010), or to reproduce Fig. 3 by rerunning the estimator on the PASS data.
