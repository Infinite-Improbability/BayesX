import logging
from collections.abc import Iterable
from dataclasses import asdict
from datetime import datetime
from os import makedirs
from multiprocessing import cpu_count
from pathlib import Path
from subprocess import run as sys_run
from subprocess import STDOUT
from typing import Optional, Union
from sys import stderr, stdout

from pybayesx.config import AnalysisConfig
from pybayesx.input import DataConfig
from pybayesx.model import (
    DeltaPrior,
    LogUniformPrior,
    Model,
    NormalPrior,
    Prior,
    Property,
    nfw_gnfw,
)
from pybayesx.plot import plot

log = logging.getLogger(__name__)


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
        """Run BayesX with a given configuration.

        :param label: Label for the run and related output, defaults to current time in
        format YYYYMMDDHHmmSS
        :type label: str, optional
        :raises ValueError: When missing priors required by model.
        """

        if path is not None:
            self.config_path = Path(path)
        elif self.config_path is None:
            self.config_path = self.base_path.joinpath(f"infile_{self.label}.inp")

        # Append filenames to paths
        self.output_path = self.base_path.joinpath("label", "out_")

        # Verify we have all required priors
        if not self.model.check_priors([p.property for p in self.priors]):
            raise ValueError("Missing required priors for model.")

        with open(self.config_path, "w") as f:
            for key, value in (
                asdict(self.data_config) | asdict(self.analysis_config)
            ).items():
                if value is None:
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
            for prior in self.priors:
                prior: Prior
                f.write(f"#{prior.property.value}\n")
                f.write(prior.export() + "\n")

            f.write("#cluster_model\n")
            f.write(str(self.model.num) + "\n")

            f.write("#root\n")
            f.write(f"'{self.output_path}'\n")

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
            log_path = self.output_path.joinpath("log.log")
            while log_path.is_file():
                i += 1
                log_path = self.output_path.joinpath(f"log{i}.log")
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
        plot_path = self.base_path.joinpath("plots")
        makedirs(plot_path)

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

        plot(self.output_path, plot_priors, plot_path)

    def analyse(self):
        self.export_infile()
        self.run()
        self.plot()


if __name__ == "__main__":
    # Write logging to shell
    log.setLevel(logging.INFO)  # log everything
    infoHandler = logging.StreamHandler(stream=stdout)
    infoHandler.setLevel(logging.INFO)  # write info to stdout
    log.addHandler(infoHandler)
    errorHandler = logging.StreamHandler(stream=stderr)
    errorHandler.setLevel(logging.WARNING)  # warnings or errors go to stderr
    log.addHandler(errorHandler)

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
