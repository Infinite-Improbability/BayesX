"""
BayesX
 - Config
   - make_config()
 - Data
   - import()
 - Output
   - plot()
   - run_report()
"""

from pathlib import Path
from typing import Union
import numpy as np
from datetime import datetime

rng: np.Generator = np.random.default_rng()


class BayesX:
    """Wrap everything needed for running BayesX on a dataset.
    Includes config, input data and results.
    
    :param label: A human-readable name used to identify this run, such as naming files. Defaults to current time in ISO8601 format.
    :type label: str

    :param path: Path to write all files to. Defaults to 'chains/{self.label}/'
    :type path: pathlib.Path

    :param config: Infile key/value pairs, excluding priors
    :type config: dict[str, any]

    :param priors: Infile key/value pairs for priors
    :type priors: dict[str, Prior]

    """

    def __init__(
        self,
        label: str = str(datetime.now().isoformat()),
        chain_root: Union[str, bytes, Path] = "chains",
    ) -> None:
        self.label: str = label
        self.path: Path = Path.join(chain_root, self.label)

        self.config: dict[str] = None
        self.priors: dict[str, Prior] = None

    def set_config(self, **kwargs) -> None:
        config: None = BayesX.default_config()

        config["root"] = self.output_path
        # TODO: Data paths
        # config["filBG"] = ''  # The root for background file
        # config["filevent"] = ''  # The root for the event file
        # config["filmask"] = ''
        # TODO: Dynamic stuff
        for k, v in kwargs.items():
            config[k] = v
        # TODO: Sanity checks - sensible values, actually exist

        # Call priors here?

        self.config = config

    def set_priors(self, *kwargs) -> None:
        priors: dict = {}
        priors["x_prior"] = Prior(0, 1, 1)
        priors["y_prior"] = Prior(0, 1, 1)
        priors["m200_prior"] = Prior(2, 1e14, 6e15, 5e15)
        priors["fgas200_prior"] = Prior(3, 0.13, 0.02, 0.131)
        priors["a_GNFW_prior"] = Prior(0, 1.062, 1.062)
        priors["b_GNFW_prior"] = Prior(0, 5.4807, 5.4807)
        priors["c_GNFW_prior"] = Prior(0, 0.3292, 0.3292)
        priors["c500_GNFW_prior"] = Prior(0, 1.156, 1.156)
        priors["alpha_model2_prior"] = Prior(1, 0.05, 5)
        priors["gamma0_poly_prior"] = Prior(1, 0.1, 0.4)
        priors["gammaR_poly_prior"] = Prior(1, 0.1, 0.4)
        priors["t0_poly_prior"] = Prior(1, 1, 5)
        priors["z_Prior"] = Prior(0, 0.5, 0.5)

        # TODO: Sanity checking

        for k, v in kwargs.items():
            priors[k] = v

        self.priors = priors

    @classmethod
    def default_config(cls) -> dict:
        # All parameters will be defined in this dict
        # Except priors, which have their own dict.
        params: dict = {}

        # The root for MEKAL data files:
        params["filion"] = Path("data/MEKAL/mekal1.dat")
        params["filrec"] = Path("data/MEKAL/mekal2.dat")
        params["filkar"] = Path("data/MEKAL/mekal3.dat")
        params["filcon"] = Path("data/MEKAL/mekal4.dat")
        params["fillin"] = Path("data/MEKAL/mekal5.dat")
        params["filcaf"] = Path("data/MEKAL/mekal6.dat")

        # The telescope
        params[
            "XrayTelescope"
        ] = "CHA"  # The first three letters of the X-ray telescope name
        params["filARF"] = Path(
            "data/simtestdatafiles/ARF_32by1.txt"
        )  # Root for ARF telescope file in .txt format.
        params["filRMF"] = Path(
            "data/simtestdatafiles/RMF_32by32.txt"
        )  # Root for RMF telescope file in .txt format.

        # Misc params
        # TODO: Set based on data
        params["n"] = 64  # Number of steps for discretising r
        params["nx"] = 64  # Number of pixels in x direction
        params["ny"] = 64  # Number of pixels in y direction
        params["xrayNbin"] = 32  # Number of energy bins
        params["xrayNch"] = 32  # Number of energy bins
        params["Aeffave"] = 250  # Average effective area of the telescope in cm^{2}
        params["xraycell"] = 0.492
        params["xrayEmin"] = 0.7
        params["xrayEmax"] = 7.0
        params["sexpotime"] = "120d3"
        params["bexpotime"] = "120d3"
        params["NHcol"] = "2.2d20"
        params["xrayBG_model"] = "8.4d-6"
        params["rmin"] = 0.01
        params["rmax"] = 0.3
        params["rlimit"] = 0.3

        params["eff"] = 0.8
        params["tol"] = 0.5
        params["nlive"] = 100
        params["seed"] = -1

        params["nCdims"] = "2"

        return params

    def export_config(self, path: Path = None) -> None:
        if not path:
            path: Path = self.path.joinpath(f"infile-{self.label}.inp")

        with open(self.path, "w") as f:
            for key, value in self.config.items():
                f.writelines(f"#{key}\n")
                f.write(str(value) + "\n")
            for key, value in self.priors.items():
                f.write(f"#{key}\n")
                f.write(value.free() + "\n")

    def load_data(self, path: Path, file_type: str = None) -> None:
        ext = Path.suffix
        if file_type:
            ext = file_type if file_type[0] == "." else "." + file_type
        if ext == ".txt":
            self.load_data_from_txt(self, path)
        elif ext == ".fits":
            self.load_data_from_fits(self, path)
        else:
            raise ValueError("Unrecognised data file extension.")

    def load_data_from_txt(self, path: Path) -> None:
        pass

    def load_data_from_fits(self, path: Path) -> None:
        pass


class Prior:
    """
    Define prior.

    True value is selected from the distribution. See readme for BayesX
    for information on type, param1 and param2.
    """

    def __init__(
        self, type_: int, param1: float, param2: float, value: float = None
    ) -> None:
        """
        Parameters
        ----------
        type_ : int
            Prior type as specified in BayesX's readme
        param1 : float
            Parameter for prior
        param2 : float
            Parameter for prior
        value : float, optional
            Override random assignment of true value
        """
        self.type: int = type_
        self.param1: float = param1
        self.param2: float = param2
        self.value: float

        if value:
            self.value = value
        else:
            self.value = self.generate_value()

    def generate_value(self) -> float:
        if self.type == 0:
            assert self.param1 == self.param2
            value = self.param1
        elif self.type == 1:
            # Uniform
            value = rng.uniform(self.param1, self.param2)
        elif self.type == 2:
            # LogUniform
            lx1 = np.log10(self.param1)
            lx2 = np.log10(self.param2)
            r = rng.random()
            value = 10 ** (lx1 + r * (lx2 - lx1))
        elif self.type == 3:
            # Gaussian
            value = rng.normal(self.param1, self.param2)
        elif self.type == 4:
            # LogGaussian
            value = rng.lognormal(self.param1, self.param2)
        return value

    def free(self) -> str:
        return f"{self.type}    {self.param1}  {self.param2}"

    def fixed(self) -> str:
        return f"0    {self.value}  {self.value}"
