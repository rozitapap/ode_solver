"""
Export the whole catalog as JSON for the website.

    python3 -m ode_catalog.export                 # writes catalog.json
    python3 -m ode_catalog.export -o path.json
    python3 -m ode_catalog.export -o -            # prints to stdout

Callables (rhs, extra_checks) are left out; the backend runs models through
CATALOG[id].request(...) and the engine.
"""

import argparse
import json
import sys
from pathlib import Path

from ode_catalog import CATALOG

SCHEMA_VERSION = 1


def catalog_to_dict():
    """The full catalog: what a "get model" call returns, for every model."""
    return {
        "schema_version": SCHEMA_VERSION,
        "models": [CATALOG[model_id].to_dict() for model_id in sorted(CATALOG)],
    }


def catalog_summary():
    """What a "list models" call returns."""
    return [CATALOG[model_id].summary_dict() for model_id in sorted(CATALOG)]


def main(argv=None):
    parser = argparse.ArgumentParser(description="Export the model catalog as JSON.")
    parser.add_argument("-o", "--output", default="catalog.json", help="output path, or - for stdout")
    args = parser.parse_args(argv)

    text = json.dumps(catalog_to_dict(), indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    if args.output == "-":
        sys.stdout.write(text)
    else:
        Path(args.output).write_text(text, encoding="utf-8")
        print(f"Wrote {len(CATALOG)} models to {args.output}")


if __name__ == "__main__":
    main()
