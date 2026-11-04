# Regression baseline

These CSVs are the numerical output of the code before the library restructuring: `ode_solver` at commit `9cf9bee` and `PINNSLearning` at commit `fff435d`.

They only show that restructuring didn't change any results. They aren't scientific validation, because they were produced by the code they're checking.

`tests/test_baseline.py` compares current output against these files. Don't regenerate them to make a failing test pass unless the numerical change is intended and understood.
