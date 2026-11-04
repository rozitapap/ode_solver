"""
Collection of mathematical models with their context, sources and validation.

To add a model, follow ode_catalog/ADDING_A_MODEL.md.
"""

from .schema import ModelEntry, ModelStatus, Parameter, Reference, Scenario
from .brain_heart.model import ENTRY as BRAIN_HEART

CATALOG = {
    "brain_heart": BRAIN_HEART,
}

__all__ = ["CATALOG", "ModelEntry", "ModelStatus", "Parameter", "Reference", "Scenario"]
