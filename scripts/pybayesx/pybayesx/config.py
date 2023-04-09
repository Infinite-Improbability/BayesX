from dataclasses import dataclass
from logging import getLogger
from pathlib import Path

log = getLogger(__name__)


@dataclass
class AnalysisConfig:
    """Dataclass hold configuration tied to general analysis and not a specific dataset.

    :param n: Number of steps for discretising r, defaults to 64
    :type n: int, optional
    :param nlive: Number of live points, defaults to 1000
    :type nlive: int, optional
    :param IS: Do Importance Nested Sampling, defaults to False
    :type IS: bool
    :param multimodal: Do more seperation, defaults to False
    :type multimodal: bool, optional
    :param eff: Target efficency, defaults to 0.8
    :type eff: float, optional
    :param tol: Tolerance value for convergence, defaults to 0.5
    :type tol: float, optional
    :param seed: Seed for RNG, negative value takes seed for from system clock,
     defaults to -1
    :type seed: int, optional
    :param updint: Number of iterations after which the output files should be updated,
     defaults to 100
    :type updint: int, optional
    :param maxmodes: Maximum number of modes for memory allocation, defaults to 20
    :type maxmodes: int, optional
    :param nCdims: Number of parameters on which cluster should be performed if mode
     seperatation is enabled. Defaults to 2
    :type nCdims: int, optional
    :param filion: Path to MEKAL data file
    :type filion: Path, optional
    :param filrec: Path to MEKAL data file
    :type filrec: Path, optional
    :param filkar: Path to MEKAL data file
    :type filkar: Path, optional
    :param filcon: Path to MEKAL data file
    :type filcon: Path, optional
    :param fillin: Path to MEKAL data file
    :type fillin: Path, optional
    :param filcaf: Path to MEKAL data file
    :type filcaf: Path, optional
    :param mass_function: Type of the mass function (with mass_function = 1, 2 & 3 for
     Evrard, Jenkins & Tinker mass functions respectively). Used only for M-z joint
     prior.
    :type mass_function: int, optional
    """

    # Data related
    n: int = 64  # Number of steps for discretising r

    # MutliNest parameters
    nlive: int = 1000  # Number of live points.
    IS: bool = False  # Do Importance Nested Sampling?
    multimodal: bool = False  # Do mode seperation
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
