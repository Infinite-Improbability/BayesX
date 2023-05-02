import logging
from collections.abc import Iterable
from dataclasses import asdict, dataclass
from datetime import datetime
from multiprocessing import cpu_count
from os import makedirs
from pathlib import Path
from subprocess import STDOUT
from subprocess import run as sys_run
from typing import Optional, Union

from .data import DataConfig
from .model import (
    DeltaPrior,
    LogUniformPrior,
    Model,
    NormalPrior,
    Prior,
    Property,
    nfw_gnfw,
)
from .plot import plot

log = logging.getLogger(__name__)


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


class Analysis:
    def __init__(
        self,
        data_config: DataConfig,
        model: Model,
        priors: Iterable[Prior],
        analysis_config: AnalysisConfig = AnalysisConfig(),
        label: str = str(datetime.now().strftime("%Y%m%d%H%M%S")),
        base_path: Union[Path, str] = Path("chains/"),
    ) -> None:
        base_path = Path(base_path)
        if not base_path.is_dir():
            log.warn(
                f"Base path {base_path} does not exist"
                + " and may be created automatically."
            )

        self.data_config = data_config
        self.model = model
        self.priors = priors
        self.analysis_config = analysis_config
        self.label = label
        self.base_path = base_path.joinpath(label)

        self.config_path: Optional[Path] = None
        self.output_path: Optional[Path] = None

    def export_infile(self, path: Optional[Union[Path, str]] = None):
        """
        Export configuration in BayesX infile format.

        :param path: Path to save infile as.
        :type label: Path or string, optional
        :raises ValueError: When missing priors required by model.
        """

        if path is not None:
            self.config_path = Path(path)
        elif self.config_path is None:
            self.config_path = self.base_path.joinpath(f"infile_{self.label}.inp")

        # Append filenames to paths
        self.output_path = self.base_path.joinpath("out")

        # Verify we have all required priors
        if not self.model.check_priors([p.property for p in self.priors]):
            raise ValueError("Missing required priors for model.")

        makedirs(self.config_path.parent, exist_ok=True)

        with open(self.config_path, "w") as f:
            for key, value in (
                asdict(self.data_config) | asdict(self.analysis_config)
            ).items():
                if value is None:
                    log.debug(f"Skipping entry {key} (is None)")
                    continue
                elif isinstance(value, float):
                    value = str(value).replace(
                        "e", "d"
                    )  # Output using fortran double notation 4.0d5 = 4.0 * 10^5
                elif isinstance(value, Path):
                    value = f"'{value}'"
                elif isinstance(value, bool):
                    if value:
                        value = "T"
                    else:
                        value = "F"

                f.write(f"#{key}\n")
                f.write(str(value) + "\n")

                log.debug(f"Written {key}: {value} to infile")

            for prior in self.priors:
                prior: Prior
                f.write(f"#{prior.property.value}\n")
                f.write(prior.export() + "\n")

            f.write("#cluster_model\n")
            f.write(str(self.model.num) + "\n")

            f.write("#root\n")
            f.write(f"'{self.output_path}'\n")

        log.info(f"Exported infile to {self.config_path}")

        return self.config_path

    def run(
        self,
        binary_path: Union[Path, str] = Path("bin/BayesX"),
        processes: int = -1,
        logfile: bool = True,
        flags: list[str] = ["--use-hwthread-cpus"],
    ):
        """Run BayesX on given analysis
        :param binary_path: Path to the Bayes-X binary,defaults to `bin/BayesX`
        :type binary_path: Union[Path, str], optional
        """

        log.info(f"Initalizing run {self.label}.")

        # Convert input data to desired formats
        binary_path = Path(binary_path)

        # TODO: Needs better handling for the case where a filename prefix is included
        if self.output_path is None:
            raise Exception  # Output path should be set by the infile
            # self.output_path = self.base_path.joinpath("out")
        elif self.output_path.is_dir():
            log.info(
                "Output path already exists. BayesX will resume from last checkpoint."
            )  # TODO: Handle this better
        elif self.output_path.parent.is_dir():
            # We're probably writing to a file with a name prefix
            log.info(
                "Output path already exists. BayesX will resume from last checkpoint."
            )  # TODO: Handle this better
        else:
            makedirs(self.output_path.parent)

        if self.config_path is None:
            raise Exception  # TODO: Actual error handling

        if processes < 0:
            processes = cpu_count()

        # write output to file
        run_args = (
            [
                "mpiexec",
                "-n",
                str(processes),
            ]
            + flags
            + [binary_path, self.config_path]
        )

        # if True write BayesX output to a file
        if logfile:
            i = 0
            log_path = self.output_path.parent.joinpath("log.log")
            while log_path.is_file():
                i += 1
                log_path = self.output_path.parent.joinpath(f"log{i}.log")
            lf = log_path.open("w")
        else:
            lf = None

        log.info(f"Launching BayesX for run {self.label}")
        start_time = datetime.now()

        # Actually run BayesXs
        sys_run(run_args, check=True, stdout=lf, stderr=STDOUT)

        log.info(
            f"BayesX finished with run {self.label} "
            f"after {(datetime.now() - start_time).total_seconds()}s"
        )

    def plot(self, plot_priors: list[str | tuple[str, float] | tuple[str, None]] = []):
        plot_path = self.base_path.joinpath(
            "plots", f"plot_tri_{datetime.now().strftime('%Y%m%d%H%M%S')}.svg"
        )
        makedirs(plot_path.parent)

        log.info(f"Creating plot at {plot_path}")

        if not plot_priors:
            # TODO: This logic seems sketchy, verify
            p_count = 1
            for p in self.priors:
                if p.type != 0:
                    # plotting script expects p001, p002, etc
                    p_label = f"p{str(p_count).zfill(3)}"
                    if p.true_value is not None:
                        plot_priors.append((p_label, p.true_value))
                    else:
                        plot_priors.append(p_label)
                    p_count += 1

        if self.output_path is None:
            raise Exception  # TODO: Error handling

        plot(self.output_path, plot_priors, plot_file=plot_path)

    def analyse(self):
        self.export_infile()
        self.run()
        self.plot()


def demo():
    log.info("Running example analysis")

    dc = DataConfig(
        filevent=Path("data/simtestdatafiles/data64by64by32.txt"),
        filBG=Path("data/simtestdatafiles/BG64by64by32.txt"),
        filARF=Path("data/simtestdatafiles/ARF_32by1.txt"),
        filRMF=Path("data/simtestdatafiles/RMF_32by32.txt"),
        nx=64,
        ny=64,
        xrayNbin=32,
        xrayNch=32,
        xraycell=0.492,
        xrayEmin=0.7,
        xrayEmax=7,
        sexpotime=300e3,
        bexpotime=300e3,
    )

    ac = AnalysisConfig(nlive=100)

    ps = [
        DeltaPrior(Property.x, 1, 1),
        DeltaPrior(Property.y, 1, 1),
        LogUniformPrior(Property.M_200, 1e14, 6e15),
        NormalPrior(Property.fg_200, 0.13, 0.02),
        DeltaPrior(Property.a_GNFW, 1.062, 1.062),
        DeltaPrior(Property.b_GNFW, 5.4807, 5.4807),
        DeltaPrior(Property.c_GNFW, 0.3292, 0.3292),
        DeltaPrior(Property.c500_GNFW, 1.156, 1.156),
        DeltaPrior(Property.z, 0.5, 0.5),
    ]

    a = Analysis(dc, nfw_gnfw, ps, ac)
    a.analyse()
