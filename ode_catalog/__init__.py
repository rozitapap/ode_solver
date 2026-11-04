"""
Collection of mathematical models with their context, sources and validation.

To add a model, follow ode_catalog/ADDING_A_MODEL.md.
"""

from .schema import ModelEntry, ModelStatus, Parameter, Reference, Scenario
from .bioprocess.model import ENTRY as BIOPROCESS
from .brain_heart.model import ENTRY as BRAIN_HEART
from .metabolite_pathway.model import ENTRY as METABOLITE_PATHWAY
from .sird.model import ENTRY as SIRD

CATALOG = {
    "bioprocess": BIOPROCESS,
    "brain_heart": BRAIN_HEART,
    "metabolite_pathway": METABOLITE_PATHWAY,
    "sird": SIRD,
}

__all__ = ["CATALOG", "ModelEntry", "ModelStatus", "Parameter", "Reference", "Scenario"]
