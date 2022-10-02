"""
BayesX
 - Config
   - make_config()
 - Data
   - import()
   - Way to trim data to certain channels, spatial ranges
 - Output
   - plot()
   - run_report()
"""
# TODO: Check compatability with earlier version of Python 3

from datetime import datetime
import logging
from pathlib import Path
from typing import Any, Union

import numpy as np
from astropy.io import fits

from .binning import bin

log = logging.getLogger(__name__)
rng = np.random.default_rng()


class Prior:
    """
    A prior as defined in BayesX.

    :param type_: Prior type as specified in BayesX's readme
    :type type_: int
    :param param1: Parameter for prior as specified in BayesX's readme for prior type.
    :type param1: float
    :param param1: See param1.
    :type param1: float
    :param true_value: True value of prior
    :type true_value: float, optional
    """

    def __init__(
        self,
        name: str,
        type_: int,
        param1: float,
        param2: float,
        true_value: float = None,
        randomise: bool = False,
    ) -> None:
        """
        :param name: Prior name as used by BayesX infile, excluding hash
        :type name: str
        :param type_: Prior type as specified in BayesX's readme
        :type type_: int
        :param param1: Parameter for prior
        :type param1: float
        :param param1: Parameter for prior
        :type param1: float
        :param true_value: True value of prior
        :type true_value: float, optional
        :param randomise: If true, pick true value from prior.
        :type randomise: bool, optional
        """
        self.name: str = name
        self.type: int = type_
        self.param1: float = param1
        self.param2: float = param2
        self.true_value: float = true_value

    def sample_prior(self) -> float:
        """Randomly select a value from prior, using the same methods as in BayesX.

        :return: Random value selected from prior.
        :rtype: float
        """
        value: float
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
        """Return prior as defined.

        :return: String of prior params formatted as in BayesX infile.
        :rtype: str
        """
        return f"{self.type}    {self.param1}  {self.param2}"

    def fixed(self) -> str:
        """Return true value of prior as fixed prior

        :return: String of fixed prior at true value formatted as in BayesX infile.
        :rtype: str
        """
        if self.true_value is None:
            raise AttributeError(
                "Unable to export fixed prior. No true value has been set."
            )
        return f"0    {self.true_value}  {self.true_value}"


class BayesX:
    """Wrap everything needed for running BayesX on a dataset.
    Includes config, input data and results.

    :param label: A human-readable name used to identify this run, such as naming files.
    Defaults to current time in ISO8601 format.
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
        """
        :param label: Human readable name for related output, defaults to
        current time in ISO8601 format.
        :type label: str, optional
        :param chain_root: Root path for output of any BayesX run, defaults to "chains"
        :type chain_root: Union[str, bytes, Path], optional
        """
        self.label: str = label
        self.path: Path = Path.join(chain_root, self.label)

        self.config: dict[str] = None
        self.priors: dict[str, Prior] = None

    def make_config(self, **kwargs) -> None:
        """
        Create config dictionary from default and dynamic values, with optional
        overrides applied through keyword arguments. For details on keys and values
        see the readme for BayesX. This overwrites any existing config.

        Validation is applied where possible to ensure reasonable inputs.
        """

        config: dict[str, Any] = BayesX.default_config()

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

        if self.config:
            log.warn("Overwriting config.")

        self.config = config

    def set_priors(self, **kwargs: Prior) -> None:
        """
        Create prior dictionary from default, with optional overrides applied through
        keyword arguments. For informnation on prior names and specifications see the
        readme for BayesX.
        """
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
    def default_config(cls) -> "dict[str, Any]":
        """Return a config dictonary with some default values. Does not include any
        dynamic values created based on the data.

        :return: Dictonary of key-value pairs as used in BayesX infile.
        :rtype: dict[str, Any]
        """
        params: dict = {}

        # The root for MEKAL data files:
        params["filion"] = Path("data/MEKAL/mekal1.dat")
        params["filrec"] = Path("data/MEKAL/mekal2.dat")
        params["filkar"] = Path("data/MEKAL/mekal3.dat")
        params["filcon"] = Path("data/MEKAL/mekal4.dat")
        params["fillin"] = Path("data/MEKAL/mekal5.dat")
        params["filcaf"] = Path("data/MEKAL/mekal6.dat")

        # The telescope
        params["XrayTelescope"] = "CHA"  # The 1st 3 letters of the X-ray telescope name
        params["filARF"] = Path(
            "data/simtestdatafiles/ARF_32by1.txt"
        )  # Root for ARF telescope file in .txt format.
        params["filRMF"] = Path(
            "data/simtestdatafiles/RMF_32by32.txt"
        )  # Root for RMF telescope file in .txt format.

        # Misc params
        # TODO: Set based on data in set_config()
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
        """Export self.config and self.priors to file

        :param path: Path to destination file. Defaults to infile-{self.label}.inp in
        the folder referenced by self.path
        :type path: Path, optional
        """
        if not path:
            path: Path = self.path.joinpath(f"infile-{self.label}.inp")

        if not (self.config and self.priors):
            log.warn("Exporting config when config or priors not set.")

        with open(self.path, "w") as f:
            for key, value in self.config.items():
                f.write(f"#{key}\n")
                f.write(str(value) + "\n")
            for key in self.model.priors():
                prior: Prior = self.priors[key]
                f.write(f"#{prior.name}\n")
                f.write(prior.free() + "\n")

        log.info(f"Exported config to {path}")

    def load_events(self, path: Path, file_type: str = None) -> None:
        """Try to intelligently handle events input from a number of file types.

        :param path: Path to events file
        :type path: Path
        :param file_type: Override type detection by providing a file extension.
        :type file_type: str, optional
        :raises ValueError: If the file type cannot be recognised.
        """
        ext = Path.suffix
        if file_type:
            ext = file_type if file_type[0] == "." else "." + file_type
        if ext == ".txt":
            self.load_events_from_txt(self, path)
        elif ext == ".fits":
            self.load_data_from_fits(self, path)
        else:
            raise ValueError("Unrecognised data file extension.")

    def load_events_from_txt(self, path: Path) -> None:
        """Load events from a BayesX-style text file.
         The event file must contains the counts per pixel per energy channel as 1-D
         array of dimension of (nx times ny times xrayNch) with each element on a
         new line.

        :param path: Path to events txt file
        :type path: Path
        """
        self.config["filevent"] = path

    def load_events_from_fits(
        self,
        path: Path,
        x_key: str = "X",
        y_key: str = "Y",
        channel_key: str = "PI",
        **kwargs,
    ) -> None:
        """Load events from a fits file.
        The file will be converted to the text format used by BayesX.

        :param path: Path to events fits file
        :type path: Path
        :param x_key: _description_, defaults to "X"
        :type x_key: str, optional
        :param y_key: _description_, defaults to "Y"
        :type y_key: str, optional
        :param energy_key: _description_, defaults to "ENERGY"
        :type energy_key: str, optional
        """
        with fits.open(path) as f:
            outpath: Path = self.path.joinpath(f"{self.label}-events.inp")
            bin(f[x_key], f[y_key], f[channel_key], outfile=outpath**kwargs)
            # TODO: Get exposure time, etc
