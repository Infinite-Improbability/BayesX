from logging import getLogger
from pathlib import Path
from typing import Iterable

import numpy as np

# TODO: Properly integrate with Mask class

log = getLogger(__name__)


class Ellipse:
    """Ellipse definition."""

    def __init__(
        self, x: float, y: float, radius1: float, radius2: float, angle: float
    ):
        super().__init__()
        self.origin_x = x
        self.origin_y = y
        self.radius1 = radius1  # width if angle = 0
        self.radius2 = radius2  # height if angle = 0
        self.angle = angle  # Radians, counterclockwise from x axis

    def ellipse_equation(self, x, y):
        t1 = (
            np.cos(self.angle) * (x - self.origin_x)
            + np.sin(self.angle) * (y - self.origin_y)
        ) ** 2

        t2 = (
            np.sin(self.angle) * (x - self.origin_x)
            - np.cos(self.angle) * (y - self.origin_y)
        ) ** 2
        return t1 / (self.radius1**2) + t2 / (self.radius2**2)

    def contains_point(self, x, y):
        # TODO: Reimplement overdraw as tolerance here
        if self.ellipse_equation(x, y) <= 1:
            return True
        else:
            return False

    @property
    def semimajor_axis(self) -> float:
        return np.max((self.radius1, self.radius2))

    def __str__(self) -> str:
        return f"ellipse: x0={self.origin_x} y0={self.origin_y}, r1={self.radius1}, r2={self.radius2}, angle={self.angle}"


def mask(
    xMin: float,
    xMax: float,
    yMin: float,
    yMax: float,
    maskFiles: Iterable[Path],
):
    # Convert lines of ellipse definitions in input file to a list of Ellipse objects.
    log.info("Loading ellipses.")
    ellipses: list[Ellipse] = []
    for filePath in maskFiles:
        with open(filePath, "r") as f:
            for line in f:
                line = line.strip()
                if "ellipse" not in line:
                    continue
                # Currently assuming pixel coordinates. In theory this isn't guaranteed.
                x, y, r1, r2, angle = [
                    float(el.strip(" )")) for el in line.split("(")[1].split(",")
                ]
                ellipses.append(Ellipse(x, y, r1, r2, angle * np.pi / 180))

    x_list = []
    y_list = []

    for ellipse in ellipses:
        log.debug(f"Processing ellipse {ellipse}")

        min_x = ellipse.origin_x - ellipse.semimajor_axis
        min_y = ellipse.origin_y - ellipse.semimajor_axis
        max_x = ellipse.origin_x + ellipse.semimajor_axis
        max_y = ellipse.origin_y + ellipse.semimajor_axis

        if min_x < xMin:
            min_x = xMin
        if min_y < yMin:
            min_y = yMin
        if max_x > xMin:
            max_x = xMax
        if max_y > yMin:
            max_y = yMax

        for x in range(int(np.floor(min_x)), int(np.ceil(max_x)) + 1):
            for y in range(int(np.floor(min_y)), int(np.ceil(max_y)) + 1):
                if ellipse.contains_point(x, y):
                    x_list.append(x)
                    y_list.append(y)

    ones = np.ones(len(x_list))

    return np.column_stack((x_list, y_list, ones))
