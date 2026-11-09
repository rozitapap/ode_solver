"""
Collection of mathematical models with their context, sources and validation.

To add a model, follow ode_catalog/ADDING_A_MODEL.md.
"""

from .schema import (
    Auxiliary,
    ConservedQuantity,
    DataPurpose,
    InputVariable,
    ModelEntry,
    ModelStatus,
    Parameter,
    Reference,
    ReferenceDataset,
    Scale,
    Scenario,
    Signal,
    SolverHints,
    StateVariable,
    validate_entry,
)
from .bioprocess.model import ENTRY as BIOPROCESS
from .brain_heart.model import ENTRY as BRAIN_HEART
from .metabolite_pathway.model import ENTRY as METABOLITE_PATHWAY
from .sird.model import ENTRY as SIRD

_ENTRIES = [BIOPROCESS, BRAIN_HEART, METABOLITE_PATHWAY, SIRD]

# Keyed by each entry's id, so the key can't drift from it.
CATALOG = {entry.id: entry for entry in _ENTRIES}

__all__ = [
    "CATALOG",
    "Auxiliary",
    "ConservedQuantity",
    "DataPurpose",
    "InputVariable",
    "ModelEntry",
    "ModelStatus",
    "Parameter",
    "Reference",
    "ReferenceDataset",
    "Scale",
    "Scenario",
    "Signal",
    "SolverHints",
    "StateVariable",
    "validate_entry",
]
