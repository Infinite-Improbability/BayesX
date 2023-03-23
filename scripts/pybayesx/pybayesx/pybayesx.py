import logging
from collections.abc import Iterable
from dataclasses import asdict
from datetime import datetime
from os import cpu_count, mkdir
from pathlib import Path
from subprocess import run as sys_run
from typing import Union

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


def run(
    data_config: DataConfig,
    analysis_config: AnalysisConfig,
    model: Model,
    priors: Iterable[Prior],
    binary_path: Union[Path, str],
    base_path: Union[Path, str],
    label: str = str(datetime.now().strftime("%Y%m%d%H%M%S")),
):
    """Run BayesX with a given configuration.

    :param data_config: Input data files and associated configuration
    :type data_config: DataConfig
    :param analysis_config: Data-independent configuration for this analysis
    :type analysis_config: AnalysisConfig
    :param model: Model being used.
    :type model: Model
    :param priors: Prior distributions. Must include all priors required by the model.
    :type priors: Iterable[Prior]
    :param binary_path: Path to the Bayes-X binary
    :type binary_path: Union[Path, str]
    :param base_path: Path to root folder for all output, e.g. chains.
     Actual output will go into a subfolder determined by the label.
    :type base_path: Union[Path, str]
    :param label: Label for the run and related output, defaults to current time in
     format YYYYMMDDHHmmSS
    :type label: str, optional
    :raises ValueError: When missing priors required by model.
    """
    # Convert input data to desired formats
    binary_path = Path(binary_path) if isinstance(binary_path, str) else binary_path
    base_path = Path(base_path) if isinstance(base_path, str) else base_path

    # Generate paths for this run based on base_path and label
    base_path = base_path.joinpath(label)
    config_path = base_path.joinpath(f"infile_{label}.inp")
    output_path = base_path.joinpath("out")
    plot_path = base_path.joinpath("plots")

    # Make output folders
    mkdir(base_path)
    mkdir(output_path)
    mkdir(plot_path)

    # Append filenames to paths
    output_path = output_path.joinpath("out_")
    plot_path = plot_path.joinpath("automatic_tri.svg")

    # Verify we have all required priors
    if not model.check_priors([i.property for i in priors]):
        raise ValueError("Missing required priors for model.")

    with open(config_path, "w") as f:
        for key, value in (asdict(data_config) | asdict(analysis_config)).items():
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
        for prior in priors:
            prior: Prior
            f.write(f"#{prior.property.value}\n")
            f.write(prior.export() + "\n")

        f.write("#cluster_model\n")
        f.write(str(model.num) + "\n")

        f.write("#root\n")
        f.write(f"'{output_path}'\n")

    log.info(f"Launching BayesX for run {label}")
    start_time = datetime.now()

    sys_run(["mpiexec", "-n", str(cpu_count()), binary_path, config_path], check=True)

    log.info(
        f"BayesX finished with run {label} "
        f"after {(datetime.now() - start_time).total_seconds()}s"
    )

    # TODO: This logic seems sketchy, verify
    plot_priors: list[str | tuple[str, float] | tuple[str, None]] = []
    p_count = 1
    for p in priors:
        if p.type != 0:
            p_label = (
                f"p{str(p_count).zfill(3)}"  # plotting script expects p001, p002, etc
            )
            if p.true_value is not None:
                plot_priors.append((p_label, p.true_value))
            else:
                plot_priors.append(p_label)
            p_count += 1

    plot(output_path, plot_priors, plot_path)


if __name__ == "main":
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
        str(datetime.now().strftime("%Y%m%d%H%M%S")),
    ]

    run(dc, ac, nfw_gnfw, ps, Path("bin/BayesX"), Path("chains/"))
