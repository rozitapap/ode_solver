"""
Parameter Perturbation Study
=============================
Runs generate_dataset() n_runs times with a model's base rate parameters
each randomly scaled by +/-10%, then saves all trajectories + parameter
values to a CSV and an overlay figure of all runs per state.

Change MODEL_NAME below to run the same study against a different
ODE system (see models/ for what's registered).
"""

import matplotlib.pyplot as plt

from data_generator import generate_dataset
from models import MODELS

MODEL_NAME = "bioprocess"
model = MODELS[MODEL_NAME]

U_init = model.default_U_init
t_pred = model.default_t_pred
n_runs = 200

df = generate_dataset(
    model,
    U_init,
    t_pred,
    base_params=model.default_params,
    n_runs=n_runs,
    mode="random",
    perturb_frac=0.1,
    seed=0,
)

fig, axes = plt.subplots(len(model.state_keys), 1, figsize=(6, 2.4 * len(model.state_keys)), sharex=True)
for _, run_df in df.groupby("run"):
    for ax, key in zip(axes, model.state_keys):
        ax.plot(run_df["t"], run_df[key], alpha=0.2, color="C0")

for ax, key in zip(axes, model.state_keys):
    ax.set_ylabel(key)
axes[-1].set_xlabel("Time (h)")
axes[0].set_title(f"{n_runs} parameter perturbations (+/-10% of base rates) -- {model.name}")
fig.tight_layout()
fig.savefig("perturbation_results.png", dpi=150)
plt.show()

df.to_csv("perturbation_results.csv", index=False)

print(f"Saved {n_runs} runs -> perturbation_results.csv, perturbation_results.png")
