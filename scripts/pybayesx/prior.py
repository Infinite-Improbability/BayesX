from enum import Enum
import logging
import numpy as np

log = logging.getLogger(__name__)
rng = np.random.default_rng()


class Property(Enum):
    x = "x_prior"  # prior on x coordinate on the sky in arcsec of the cluster center
    y = "y_prior"  # prior on y coordinate on the sky in arcsec of the cluster center
    M_200 = "m200_prior"  # prior on cluster mass within the overdensity radius of r_{200} (M_{200}) in M_{\sun}
    fg_200 = "fgas200_prior"  # prior on cluster gas mass fraction within the overdensity radius of r_{200} (f_{g,200})
    a_GNFW = "a_GNFW_prior"  # prior on slope parameter "a" in GNFW pressure profile
    b_GNFW = "b_GNFW_prior"  # prior on slope parameter "b" in GNFW pressure profile
    c_GNFW = "c_GNFW_prior"  # prior on slope parameter "c" in GNFW pressure profile
    c500_GNFW = "c500_GNFW_prior"  # prior on gas concentration parameter, c_{500}
    alpha_model2 = "alpha_model2_prior"  # prior on Einasto shape parameter, alpha
    gamma0_poly = (
        "gamma0_poly_prior"  # prior gamma0, free parameter in polytropic model
    )
    gammaR_poly = (
        "gammaR_poly_prior"  # prior gammaR, free parameter in polytropic model
    )
    t0_poly = "t0_poly_prior"  # prior T0, free parameter in polytropic model
    z = "z_Prior"  # prior on cluster redshift


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
        property: Property,
        type_: int,
        param1: float,
        param2: float,
        true_value: float = None,
        randomise: bool = False,
    ) -> None:
        """
        :param property: Property described by this prior
        :type property: Property
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
        self.property: Property = property
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
        else:
            raise Exception("Invalid prior type")
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


class Uniform(Prior):
    def __init__(
        self,
        property: Property,
        type_: int,
        param1: float,
        param2: float,
        true_value: float = None,
        randomise: bool = False,
    ) -> None:
        super().__init__(property, type_, param1, param2, true_value, randomise)
