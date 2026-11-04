"""
Collection of mathematical models with their context, sources and validation.

To add a model, follow ode_catalog/ADDING_A_MODEL.md.
"""

from .schema import ModelEntry, ModelStatus, Parameter, Reference, Scenario

CATALOG = {}

__all__ = ["CATALOG", "ModelEntry", "ModelStatus", "Parameter", "Reference", "Scenario"]
