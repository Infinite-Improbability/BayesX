from logging import getLogger
from pathlib import Path
from typing import Iterable, Optional

import numpy as np

# TODO: Properly integrate with Mask class

log = getLogger(__name__)


class Ellipse:
    """Ellipse definition."""

    def __init__(self, x, y, radius1, radius2, angle):
        super().__init__()
        self.x = x
        self.y = y
        self.radius1 = radius1
        self.radius2 = radius2
        self.angle = angle  # Radians

    def generate_matrix(self, overdraw=0):
        """Create a square array of 0 and 1, with 1 representing the shape of the
        ellipse."""

        # Add one to ensure we don't miss anything on the edge.
        radius = max(self.radius1, self.radius2) + 1
        # Add 1 to account for Python's end behaviour (last element of 1:10 is 9, etc)
        r = np.arange(-radius, radius + 1, 1)

        x = np.tile(r, (len(r), 1))  # Create matrix of x values for each point
        y = x.T  # Create same matrix for y values

        # For every point check if in ellipse. If so, set to 1.
        return (
            x * np.cos(self.angle) + y * np.sin(self.angle)
        ) ** 2 / self.radius1**2 + (
            x * np.sin(self.angle) - y * np.cos(self.angle)
        ) ** 2 / self.radius2**2 <= (
            1 + overdraw
        )


def mask(
    xMin: float,
    xMax: float,
    yMin: float,
    yMax: float,
    maskFiles: Iterable[Path],
    overdraw: float,
    output_file: Optional[Path] = None,
    display: bool = True,
):
    def coord2index(x, y, referenceCoordinates=(xMin, yMin)):
        """
        Convert a coordinate to index.

        Intended to handle input coordinates that don't align with matrix coordinates
        due to origin position.

        Aligns referenceCoordinates to index [0,0]. Reference coordinates should be
        (xMin, yMin.)

        x-axis points down.
        """

        # Adds zero index to each coordinate and rounds
        # TODO: Needs confirmation of correct rounding behaviour.
        return [
            int(np.ceil(coord - ref))
            for coord, ref in zip((x, y), referenceCoordinates)
        ]

    size = (int(yMax - yMin), int(xMax - xMin))
    # TODO: Custom step size support. Currently assumes step of 1.

    # Convert lines of ellipse definitions in input file to a list of Ellipse objects.
    log.info("Loading ellipses.")
    ellipses = []
    for filePath in maskFiles:
        with open(filePath, "r") as f:
            for line in f:
                line = line.strip()
                if line[:8] != "ellipse":
                    continue
                # Currently assuming pixel coordinates. In theory this isn't guaranteed.
                x, y, r1, r2, angle = [
                    float(el.strip(" )")) for el in line.split("(")[1].split(",")
                ]
                ellipses.append(Ellipse(x, y, r1, r2, angle * np.pi / 180))

    # Print ellipses on a matrix
    log.info("Generating mask array.")
    # Generate matrix of zeros (no masking) covering the entire region of data
    mask = np.full(size, False)

    for el in ellipses:
        # Convert ellipse into matrix
        elM = el.generate_matrix(overdraw)

        # Get indices to place the ellipse on the mask matrix
        radius = elM.shape[0] / 2  # elM is square
        x0, y0 = coord2index(el.x - radius, el.y - radius)
        x1, y1 = coord2index(el.x + radius, el.y + radius)

        # If the ellipse goes off the edge of the region crop it
        if x0 < 0:  # Left
            elM = elM[:, -x0:]
            x0 = 0
        if y0 < 0:  # Top
            elM = elM[-y0:, :]
            y0 = 0
        if x1 > size[1]:  # Right
            elM = elM[:, : -(x1 - size[1])]
            x1 = size[1]
        if y1 > size[0]:  # Bottom
            elM = elM[: -(y1 - size[0]), :]
            y1 = size[0]
        if 0 in elM.shape:
            continue

        # Add matrix on mask by numpy.logical_or
        mask[y0:y1, x0:x1] = mask[y0:y1, x0:x1] | elM

    log.info("Preparing array for export.")
    # Rearrange for export
    x = np.tile(np.arange(xMin, xMax), size[0])  # [1,2,3,1,2,3,...]
    y = np.repeat(np.arange(yMin, yMax), size[1])  # [1,1,1,2,2,2,...]
    mask1D = np.column_stack((x, y, np.ravel(mask)))
    # Don't need to export umasked regions
    mask1D = mask1D[mask1D[:, 2] == 1]

    if output_file is not None:
        filename_parts = output_file.name.strip().split(".")
        file_type = filename_parts[-1]

        if file_type == "gz":
            file_type = filename_parts[-2]
        if file_type == "txt":
            log.info("Exporting to .txt file.")
            np.savetxt(output_file, mask1D, fmt=["%.10g", "%.10g", "%.1d"])
        elif file_type == "npy":
            log.info("Exporting to .npy file.")
            np.save(output_file, mask1D)
        elif file_type == "matxt":
            log.info("Exporting entire mask matrix to text file.")
            np.savetxt(output_file, mask, fmt="%.1d")
        else:
            log.info("Exporting to compressed .npz file.")
            np.savez_compressed(output_file, mask1D)
        log.info("File saved as {}.".format(output_file.name))

    if display:
        # Visual output
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            log.warn("Couldn't import matplotlib, exiting without plotting.")
            exit()

        # Trim to masked region
        for i in range(0, size[0]):
            if mask[i, :].any():
                mask = mask[i:, :]
                break
        for j in range(1, size[0]):
            if mask[-j, :].any():
                j -= 1
                mask = mask[:-j, :]
                break
        for k in range(0, size[1]):
            if mask[:, k].any():
                mask = mask[:, k:]
                break
        for n in range(1, size[1]):
            if mask[:, -n].any():
                mask = mask[:, :-n]
                break

        # Convert indices to coordinates
        yMin += i  # type: ignore
        yMax -= j  # type: ignore
        xMin += k  # type: ignore
        xMax -= n  # type: ignore

        # Get ellipse centres
        ox = []
        oy = []

        for el in ellipses:
            ox.append(el.x)
            oy.append(el.y)

        fig, ax = plt.subplots()

        ax.imshow(mask, extent=(xMin, xMax, yMin, yMax), origin="lower")
        ax.scatter(ox, oy, marker="+")  # type: ignore

        plt.show()

    return mask1D
