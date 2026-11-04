"""
ODE model registry.

Add a new system by creating models/<name>.py with an ODESystemSpec instance
named MODEL, then registering it below.
"""

from . import bioprocess, brain_heart, metabolite_pathway, sird

MODELS = {
    "bioprocess": bioprocess.MODEL,
    "brain_heart": brain_heart.MODEL,
    "metabolite_pathway": metabolite_pathway.MODEL,
    "sird": sird.MODEL,
}
