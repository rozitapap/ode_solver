"""
Run one scenario of a catalog model, print a summary, and save a figure.

    python3 experiments/single_run.py --model brain_heart --scenario setpoint_step

Results go to outputs/single_run/<model>/ unless --output-dir is given.
"""

import argparse
from pathlib import Path

from ode_catalog import CATALOG
from ode_solver import Pipeline

OUTPUTS = Path(__file__).resolve().parent.parent / "outputs"


def main():
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument("--model", default="bioprocess", choices=sorted(CATALOG))
    parser.add_argument("--scenario", help="scenario name; defaults to the model's default scenario")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--show", action="store_true", help="also show the figure")
    args = parser.parse_args()

    entry = CATALOG[args.model]
    scenario = entry.scenario(args.scenario)
    output_dir = args.output_dir or OUTPUTS / "single_run" / args.model

    run = Pipeline().run(
        entry.spec,
        entry.request(scenario.name),
        plot_path=str(output_dir / f"{scenario.name}.png"),
        show_plot=args.show,
        **entry.units,
    )
    if not run.ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
