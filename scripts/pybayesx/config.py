from dataclasses import dataclass
from pathlib import Path


@dataclass
class DataConfig:
    # Input data
    filBG: Path
    filevent: Path
    filmask: Path
    filARF: Path  # in txt format
    filRMF: Path  # in text format

    nx: int  # Number of pixels in x direction
    ny: int  # Number of pixels in y direction
    xrayNbin: int  # Number of energy bins
    xrayNch: int  # Number of energy bins

    xraycell: float  # Spatial pixel size in arcsecond
    xrayEmin: float  # Minimum value of the energy range in keV
    xrayEmax: float  # Maximum value of the energy range in keV
    sexpotime: float  # Source exposure time in second
    bexpotime: float  # Background exposure time in second

    NHcol: float = 4.0e20  # Hydrogen column density in cm^2
    xrayBG_model: float = (
        8.4e-6  # Predicted background rate at each pixel in counts cm^-2 arcmin^-2s^-1,
    )

    rmin: float = 0.01  # Minimum radius, Mpc
    rmax: float = None  # Maximum radius for xray emission and GNFW model, Mpc
    rlimit: float = (
        None  # Used to calculate logr, may need to be slightly higher than rmax, Mpc
    )


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
    nCdims: int = 2  # no. of parameters on which clustering should be performed if mode separation is enabled, default value: 2

    # Root for MEKAL data files
    filion: Path = Path("data/MEKAL/mekal1.dat")
    filrec: Path = Path("data/MEKAL/mekal2.dat")
    filkar: Path = Path("data/MEKAL/mekal3.dat")
    filcon: Path = Path("data/MEKAL/mekal4.dat")
    fillin: Path = Path("data/MEKAL/mekal5.dat")
    filcaf: Path = Path("data/MEKAL/mekal6.dat")
