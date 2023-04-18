# TODO: Log to stdout

"""A Python wrapper to make using BayesX simpler."""
import importlib.metadata
from logging import NullHandler, getLogger

# We include redundant as statements to indicate to Pylance that these are
# supposed to be public so it doesn't complain about unused imports
from .analysis import Analysis as Analysis
from .analysis import AnalysisConfig as AnalysisConfig
from .data import ARF as ARF
from .data import RMF as RMF
from .data import DataConfig as DataConfig
from .data import Events as Events
from .data import Mask as Mask
from .model import DeltaPrior as DeltaPrior
from .model import LogNormalPrior as LogNormalPrior
from .model import LogUniformPrior as LogUniformPrior
from .model import Model as Model
from .model import NormalPrior as NormalPrior
from .model import Prior as Prior
from .model import Property as Property
from .model import UniformPrior as UniformPrior
from .model import nfw_einasto as nfw_einasto
from .model import nfw_gnfw as nfw_gnfw
from .model import polytropic as polytropic
from .plot import plot as plot

__version__ = importlib.metadata.version("pybayesx")

# Set default logging handler to avoid "No handler found" warnings.
getLogger(__name__).addHandler(NullHandler())
