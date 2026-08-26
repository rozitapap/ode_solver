"""
ODE model registry.

Add a new system by creating models/<name>.py with an ODEModel instance
named MODEL, then registering it below.
"""

from . import bioprocess, sird

MODELS = {
    "bioprocess": bioprocess.MODEL,
    "sird": sird.MODEL,
}
