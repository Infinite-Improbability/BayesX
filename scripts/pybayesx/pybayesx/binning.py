from logging import getLogger
from pathlib import Path
from typing import Any, Optional

import numpy as np
from numpy.typing import ArrayLike

log = getLogger(__name__)


def bin(
    x: ArrayLike,
    y: ArrayLike,
    channel: ArrayLike,
    n_bins: int,
    cellsize: float,
    x0: Optional[float] = None,
    y0: Optional[float] = None,
    n_channels: Optional[int] = None,
    outfile: Optional[Path] = None,
    mask: bool = False,
) -> "np.ndarray[Any, np.dtype[np.float64]]":
    """Bin data given spatial coordinates and channel coordinates of events.
    Bins are only spatial and do not reduce channel count.
    This will crop points if they are not within `n_bins/2` bins of the center.
    If `mask = True` then instead bin passed on mask status.

    :param x: 1D sequence of x coordinates of events. Should have an entry for every point.
    :type x: numpy.typing.ArrayLike
    :param y: Sequence of y coordinates of events, assumed to use the same units and order as x.
    :type y: numpy.typing.ArrayLike
    :param channel: Sequence of channel of events. If `mask` is True then this value
     should be 1 for all masked spatial coordinates and 0 otherwise.
    :type channel: numpy.typing.ArrayLike
    :param nbins: The number of bins along each spatial axis
    :type nbins: int
    :param cellsize: The size of a bin, in units matching x and y
    :type cellsize: float
    :param x0: Central (origin) point along x axis in data, defaults to None. If None
     then the midpoint of x is used.
    :type x0: float, optional
    :param y0: See `x0`, defaults to None
    :type y0: float, optional
    :param n_channels: Number of possible channels, defaults to None. If None, then the
     maximum value of `channel` is used.
     If `mask` is True then a positive integer value is required.
    :type n_channels: int, optional
    :param outfile: Path of file to export binned data to, defaults to None. If None
     export is skipped.
    :type outfile: Path, optional
    :param mask: If True then the channel sequence is treated as a 1/0 (T/F) boolean
     definining masked regions. Defaults to False.
    :type mask: bool, optional

    :raises ValueError: If the value of max_chan is <=0 or (if in mask mode) if not set.
    :return: A 1D array of counts, of length `n_bins*n_bins*chan_max`. Nested as
     x_bin_index(y_bin_index)

    :rtype: np.ndarray[Any, np.dtype[np.float64]]
    """

    # Coerce input
    x = np.asanyarray(x)
    y = np.asanyarray(y)
    channel = np.asanyarray(channel)

    # Do some input checking
    if len(x) != len(y):
        raise ValueError("Coordinate lists x and y have different lengths.")
    if x.ndim != 1 or y.ndim != 1:
        raise ValueError("x or y array is not one dimensional.")

    # Set no. of channels based on largest channel number present
    if mask and (n_channels is None or n_channels <= 0):
        raise ValueError("If binning a mask, max_chan is a required argument.")
    elif n_channels is None:
        n_channels = int(np.max(channel))

    # If origins are not provided default to centre of data
    if x0 is None:
        x0 = (np.max(x) + np.min(x)) / 2
    if y0 is None:
        y0 = (np.max(y) + np.min(y)) / 2

    # Set relative coordinates
    dx = x - x0
    dy = y - y0

    # Initalize output array
    counts = np.zeros((n_bins, n_bins, n_channels))

    # Set center coordinates in bin basis
    i0 = n_bins / 2
    j0 = n_bins / 2

    # Loop over every cell
    # The x and y arrays contain an entry for every point so we don't need a nested iteration over dy
    for k in range(len(dx)):
        # Get cell coordinates in original basis, relative to centre
        u = dx[k]
        v = dy[k]

        # Get bin indices
        signx = 1.0
        if u < 0.0:
            signx = -1
        signy = 1.0
        if v < 0.0:
            signy = -1

        # Convert coordinates to bin basis
        i = int(np.floor(u / cellsize + signx * 0.5) + i0)
        j = int(np.floor(v / cellsize + signy * 0.5) + j0)

        # If bin indices out of range then abort
        if (i < 0) | (i >= n_bins):
            continue
        if (j < 0) | (j >= n_bins):
            continue

        # Increments count for bin
        if mask:
            counts[i, j, :] = np.maximum(counts[i, j, :], channel[k])
            # This should cause the mask file to match the structure of the binned data.
            # The file is chan_max times longer, but we avoid having to expand
            # the mask across all channels in BayesX.
        else:
            # Get energy channel and convert to index by subtracting 1
            vr = int(channel[k]) - 1
            counts[i, j, vr] += 1

    counts_1d = counts.ravel()

    if outfile:
        np.savetxt(outfile, counts_1d, "%5d")

    return counts_1d
