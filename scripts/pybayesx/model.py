from typing import Iterable
from .priors import Property
import logging

log = logging.getLogger(__name__)


class Model:
    def __init__(
        self,
        name: str,
        num: int,
        required_priors: Iterable[Property],
    ) -> None:
        self.name = name
        self.num = num
        self.required_priors: Iterable[Property] = required_priors

    def check_priors(self, input_priors: Iterable[Property]) -> bool:
        for prop in self.required_priors:
            if prop not in input_priors:
                return False
        return True
