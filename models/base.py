"""
Generic ODE model interface.

Every concrete model (bioprocess, SIRD, ...) builds one ODEModel instance
exposing the same shape, so the solver and dataset generator never need to
know which system they're integrating.
"""

from dataclasses import dataclass
from typing import Callable, Mapping, Sequence


@dataclass
class ODEModel:
    name: str
    state_keys: list[str]
    param_names: tuple[str, ...]
    rhs: Callable          # rhs(X, t, *params) -> list of derivatives, len == len(state_keys)
    default_params: Mapping[str, float]
    default_U_init: Sequence[float]
    default_t_pred: Sequence[float]
