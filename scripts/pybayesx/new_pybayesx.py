import logging
from collections.abc import Iterable
from dataclasses import asdict
from datetime import datetime
from os import cpu_count, mkdir
from pathlib import Path
from subprocess import run as sys_run

from config import AnalysisConfig, DataConfig
from model import Model, nfw_gnfw
from priors import Delta, LogUniform, Normal, Prior, Property

log = logging.getLogger(__name__)

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
    Delta(Property.x, 1, 1),
    Delta(Property.y, 1, 1),
    LogUniform(Property.M_200, 1e14, 6e15),
    Normal(Property.fg_200, 0.13, 0.02),
    Delta(Property.a_GNFW, 1.062, 1.062),
    Delta(Property.b_GNFW, 5.4807, 5.4807),
    Delta(Property.c_GNFW, 0.3292, 0.3292),
    Delta(Property.c500_GNFW, 1.156, 1.156),
    Delta(Property.z, 0.5, 0.5),
]


def run(
    data_config: DataConfig,
    analysis_config: AnalysisConfig,
    priors: Iterable[Prior],
    model: Model,
    binary_path: Path,
    base_path: Path,
    label: str = str(datetime.now().strftime("%Y%m%d%H%M%S")),
):
    base_path = base_path.joinpath(label)
    config_path = base_path.joinpath(f"infile_{label}.inp")
    output_path = base_path.joinpath("out")
    plot_path = base_path.joinpath("plots")

    mkdir(base_path)
    mkdir(output_path)
    mkdir(plot_path)

    output_path = output_path.joinpath("out_")  # add filename

    if not model.check_priors([i.property for i in priors]):
        raise Exception("Missing required priors for model.")

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
        f"BayesX finished with run {label} after {(datetime.now() - start_time).total_seconds()}s"
    )

    plot_priors = []
    true_priors = ["--true"]
    p_count = 1
    for p in priors:
        if p.type != 0:
            plot_priors.append(
                f"p{str(p_count).zfill(3)}"
            )  # plotting script expects p001, p002, etc
            if p.true_value is not None:
                true_priors.append(str(p.true_value))
            p_count += (
                1  # TODO: Have check in plot script to ensure it plots correct priors
            )

    if len(true_priors) > 1:
        plot_priors += true_priors

    sys_run(
        ["python3", "scripts/auto_plot_tri.py", output_path] + plot_priors,
        check=True,
    )


run(dc, ac, ps, nfw_gnfw, Path("bin/BayesX"), Path("chains/"))
