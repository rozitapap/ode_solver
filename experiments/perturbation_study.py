"""
Parameter perturbation study: run a catalog model many times with each parameter
randomly scaled, then save all trajectories and an overlay figure.

    python3 experiments/perturbation_study.py --model sird

Results go to outputs/perturbation/<model>/ unless --output-dir is given.
"""

import argparse
from pathlib import Path

from ode_catalog import CATALOG
from ode_solver import generate_dataset, plot_sweep

OUTPUTS = Path(__file__).resolve().parent.parent / "outputs"


def main():
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument("--model", default="bioprocess", choices=sorted(CATALOG))
    parser.add_argument("--scenario", help="scenario name; defaults to the model's default scenario")
    parser.add_argument("--n-runs", type=int, default=200)
    parser.add_argument("--perturb-frac", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--show", action="store_true", help="also show the figure")
    args = parser.parse_args()

    entry = CATALOG[args.model]
    output_dir = args.output_dir or OUTPUTS / "perturbation" / args.model
    output_dir.mkdir(parents=True, exist_ok=True)

    df = generate_dataset(
        entry.spec,
        entry.request(args.scenario),
        n_runs=args.n_runs,
        mode="random",
        perturb_frac=args.perturb_frac,
        seed=args.seed,
    )

    csv_path = output_dir / "trajectories.csv"
    png_path = output_dir / "trajectories.png"
    df.to_csv(csv_path, index=False)
    plot_sweep(
        df,
        entry.spec,
        save_path=str(png_path),
        show=args.show,
        title=f"{args.n_runs} runs, parameters +/-{args.perturb_frac:.0%} -- {entry.spec.name}",
        **entry.units,
    )

    n_failed = args.n_runs - df[df["success"]]["run"].nunique()
    print(f"Saved {args.n_runs} runs ({n_failed} failed) -> {csv_path}, {png_path}")


if __name__ == "__main__":
    main()
