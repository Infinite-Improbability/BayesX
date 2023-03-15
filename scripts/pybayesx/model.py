import logging
from enum import Enum
from typing import Iterable

from .priors import Property

log = logging.getLogger(__name__)


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
        Property.alpha_model2,
        Property.z,
    ],
    name="Einasto",
)

polytropic = Model(
    num=1,
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
