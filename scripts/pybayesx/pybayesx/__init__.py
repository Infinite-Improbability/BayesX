# TODO: Log to stdout

"""A Python wrapper to make using BayesX simpler."""
import importlib.metadata
from logging import NullHandler, getLogger

from .analysis import Analysis, AnalysisConfig
from .data import ARF, RMF, DataConfig, Events, Mask
from .model import (
    DeltaPrior,
    LogNormalPrior,
    LogUniformPrior,
    Model,
    NormalPrior,
    Prior,
    Property,
    UniformPrior,
    nfw_einasto,
    nfw_gnfw,
    polytropic,
)
from .plot import plot

__version__ = importlib.metadata.version("pybayesx")

# Set default logging handler to avoid "No handler found" warnings.
getLogger(__name__).addHandler(NullHandler())
