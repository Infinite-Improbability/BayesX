from abc import ABC, abstractmethod
from collections.abc import Iterable
from enum import Enum
from logging import getLogger
from typing import Optional

import numpy as np

log = getLogger(__name__)
rng = np.random.default_rng()


class Property(Enum):
    x = "x_prior"  # prior on x coordinate on the sky in arcsec of the cluster center
    y = "y_prior"  # prior on y coordinate on the sky in arcsec of the cluster center
    M_200 = "m200_prior"  # prior on cluster mass within the overdensity radius of
    # r_{200} (M_{200}) in M_{\sun}
    fg_200 = "fgas200_prior"  # prior on cluster gas mass fraction within the
    # overdensity radius of r_{200} (f_{g,200})
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


class Prior(ABC):
    """
    A prior as defined in BayesX.

    This class is abstract and so cannot be directly instantiated.
    Subclasses must implement the `sample()` method.

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
        true_value: Optional[float] = None,
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
        self.true_value: Optional[float] = true_value

    @abstractmethod
    def sample(self) -> float:
        """Randomly select a value from prior, using the same methods as in BayesX.

        :return: Random value selected from prior.
        :rtype: float
        """

    def export(self, fixed=False) -> str:
        """Return prior as defined.

        :param fixed: Export true value as delta prior, defaults to False
        :type fixed: bool, optional
        :raises AttributeError: When attempting to export fixed prior without having set
         the true value
        :return: String describing prior as used in BAYES-X infile.
        :rtype: str
        """
        if fixed:
            if self.true_value is None:
                raise AttributeError(
                    "Unable to export fixed prior because no true value has been set."
                )
            return f"0    {self.true_value}  {self.true_value}"
        else:
            return f"{self.type}    {self.param1}  {self.param2}"


class DeltaPrior(Prior):
    def __init__(
        self,
        property: Property,
        value: float,
        true_value: Optional[float] = None,
        randomise: bool = False,
    ) -> None:
        super().__init__(property, 0, value, value, true_value, randomise)

    def sample(self) -> float:
        return self.param1


class UniformPrior(Prior):
    def __init__(
        self,
        property: Property,
        min: float,
        max: float,
        true_value: Optional[float] = None,
        randomise: bool = False,
    ) -> None:
        super().__init__(property, 1, min, max, true_value, randomise)

    def sample(self) -> float:
        return rng.uniform(self.param1, self.param2)


class LogUniformPrior(Prior):
    def __init__(
        self,
        property: Property,
        min: float,
        max: float,
        true_value: Optional[float] = None,
        randomise: bool = False,
    ) -> None:
        super().__init__(property, 2, min, max, true_value, randomise)

    def sample(self) -> float:
        lx1 = np.log10(self.param1)
        lx2 = np.log10(self.param2)
        r = rng.random()
        return 10 ** (lx1 + r * (lx2 - lx1))


class NormalPrior(Prior):
    def __init__(
        self,
        property: Property,
        mean: float,
        standard_deviation: float,
        true_value: Optional[float] = None,
        randomise: bool = False,
    ) -> None:
        super().__init__(property, 3, mean, standard_deviation, true_value, randomise)

    def sample(self) -> float:
        return rng.normal(self.param1, self.param2)


class LogNormalPrior(Prior):
    def __init__(
        self,
        property: Property,
        mean: float,
        width: float,
        true_value: Optional[float] = None,
        randomise: bool = False,
    ) -> None:
        super().__init__(property, 4, mean, width, true_value, randomise)

    def sample(self) -> float:
        return rng.lognormal(self.param1, self.param2)


# TODO: Missing joint prior (see readme)


class Model:
    def __init__(
        self, num: int, required_priors: Iterable[Property], name: str = ""
    ) -> None:
        self.name = name
        self.num = num
        self.required_priors: Iterable[Property] = required_priors

    def check_priors(self, input_priors: Iterable[Property]) -> bool:
        for prop in self.required_priors:
            if prop not in input_priors:
                return False
        return True


nfw_gnfw = Model(
    num=1,
    required_priors=[
        Property.x,
        Property.y,
        Property.M_200,
        Property.fg_200,
        Property.a_GNFW,
        Property.b_GNFW,
        Property.c_GNFW,
        Property.c500_GNFW,
        Property.z,
    ],
    name="NFW-GNFW",
)

nfw_einasto = Model(
    num=2,
    required_priors=[
        Property.x,
        Property.y,
        Property.M_200,
        Property.fg_200,
        Property.a_GNFW,
        Property.b_GNFW,
        Property.c_GNFW,
        Property.c500_GNFW,
        Property.alpha_model2,
        Property.z,
    ],
    name="Einasto",
)

polytropic = Model(
    num=3,
    required_priors=[
        Property.x,
        Property.y,
        Property.M_200,
        Property.gamma0_poly,
        Property.gammaR_poly,
        Property.t0_poly,
        Property.z,
    ],
    name="NFW-GNFW",
)
