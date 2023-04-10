# TODO: Log to stdout

"""A Python wrapper to make using BayesX simpler."""
from .analysis import Analysis, AnalysisConfig
from .data import ARF, RMF, DataConfig, Events, Mask, load_all_from_fits
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
