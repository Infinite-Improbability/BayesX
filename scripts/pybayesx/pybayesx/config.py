from dataclasses import dataclass
from logging import getLogger
from pathlib import Path

log = getLogger(__name__)


@dataclass
class AnalysisConfig:
    # Data related
    n: int = 64  # Number of steps for discretising r

    # MutliNest parameters
    IS: bool = False  # Do Importance Nested Sampling?
    multimodal: bool = False  # Do mode seperation
    nlive: int = 1000  # Number of live points.
    eff: float = 0.8  # Target efficiency
    tol: float = 0.5  # Tolerance value for convergence
    seed: int = -1  # Seed for RNG, negative value takes seed from system clock
    updint: int = (
        100  # No. of iterations after which the output files should be updated
    )
    maxmodes: int = 20  # Maximum no. of modes (for memory allocation)
    nCdims: int = 2  # no. of parameters on which clustering should be performed if
    # mode separation is enabled, default value: 2

    # Root for MEKAL data files
    filion: Path = Path("data/MEKAL/mekal1.dat")
    filrec: Path = Path("data/MEKAL/mekal2.dat")
    filkar: Path = Path("data/MEKAL/mekal3.dat")
    filcon: Path = Path("data/MEKAL/mekal4.dat")
    fillin: Path = Path("data/MEKAL/mekal5.dat")
    filcaf: Path = Path("data/MEKAL/mekal6.dat")

    mass_function: int = 2  # Type of the mass function (with mass_function = 1, 2 & 3
    # for Evrard, Jenkins & Tinker mass functions respectively)
