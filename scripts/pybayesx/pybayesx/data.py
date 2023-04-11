from __future__ import annotations  # python 3.7 and up

from dataclasses import dataclass
from logging import getLogger
from pathlib import Path
from typing import Any, Optional, Union

import numpy as np
from astropy.io import fits
from astropy.io.fits.hdu import PrimaryHDU
from numpy.typing import ArrayLike

log = getLogger(__name__)

# TODO: Refactor into BinnableData subclassS


class Data:
    """Abstract base class for source data"""

    def __init__(self, data: ArrayLike) -> None:
        """Basic initialiser

        :param data: Array on input data
        :type data: ArrayLike
        """
        self.data = np.array(data)
        self.path = None  # path to data in format ready for BayesX

    def bin(
        self, n_bins: int, cellsize: float, outfile: Optional[Path] = None, **kwargs
    ):
        """Bin arbitary data with three dimensions in first two dimensions.

        :param n_bins: Amount of bins in one dimension
        :type n_bins: int
        :param cellsize: Size of a single bin, in units of source data
        :type cellsize: float
        :param outfile: Path to write binned data to, defaults to None
        :type outfile: Optional[Path], optional
        :return: 1D array of binned data
        :rtype: np.ndarray[Any, np.dtype[np.float64]]
        """
        self.path = outfile
        b = self._bin(
            self.data[:, 0],
            self.data[:, 1],
            self.data[:, 2],
            n_bins,
            cellsize,
            outfile=outfile,
            **kwargs,
        )
        self.nx, self.ny, self.n_channels = b.shape
        return b

    def _bin(
        self,
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

        The current implementation of this function leaves a cross artefact on the
        binned data.

        :param x: 1D sequence of x coordinates of events. Should have an entry for every
        point.
        :type x: numpy.typing.ArrayLike
        :param y: Sequence of y coordinates of events, assumed to use the same units and
        order as x.
        :type y: numpy.typing.ArrayLike
        :param channel: Sequence of channels of events, assumed to use the same order
         as y. If `mask` is True then this value  should be 1 for all masked spatial
         coordinates and 0 otherwise.
        :type channel: numpy.typing.ArrayLike
        :param nbins: The number of bins along each spatial axis
        :type nbins: int
        :param cellsize: The size of a bin, in units matching x and y
        :type cellsize: float
        :param x0: Central (origin) point along x axis in data, defaults to None. If
         None then the midpoint of x is used.
        :type x0: float, optional
        :param y0: See `x0`, defaults to None
        :type y0: float, optional
        :param n_channels: Number of possible channels, defaults to None. If None, then
         the maximum value of `channel` is used. If `mask` is True then a positive
         integer value is required.
        :type n_channels: int, optional
        :param outfile: Path of file to export binned data to, defaults to None. If None
        export is skipped.
        :type outfile: Path, optional
        :param mask: If True then the channel sequence is treated as a 1/0 (T/F) boolean
        definining masked regions. Defaults to False.
        :type mask: bool, optional

        :raises ValueError: If the value of max_chan is <=0 or (if in mask mode) if not
         set.

        :return: A 1D array of counts, of length `n_bins*n_bins*chan_max`. Nested as
        x_bin_index(y_bin_index)
        :rtype: np.ndarray[Any, np.dtype[np.float64]]
        """

        # TODO: Rewrite for being class method
        # And fix the cross problem

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
        # The x and y arrays contain an entry for every point so we don't need a nested
        # iteration over dy
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
                # This should cause the mask file to match the structure of the binned
                # data.
                # The file is chan_max times longer, but we avoid having to expand
                # the mask across all channels in BayesX.
            else:
                # Get energy channel and convert to index by subtracting 1
                vr = int(channel[k]) - 1
                counts[i, j, vr] += 1

        counts_1d = counts.ravel()

        if outfile:
            np.savetxt(outfile, counts_1d, "%5d")

        return counts

    def __str__(self) -> str:
        if self.path is None:
            raise ValueError("Data class has no file set. Please export to file.")
        return str(self.path)

    # @abstractmethod
    # def export():
    #     raise NotImplementedError


class Events(Data):
    def __init__(self, data: ArrayLike, background: bool, exposure_time: float) -> None:
        """X-ray event data

        :param data: Source data array. Three columns, with each row containing (x, y, channel).
        :type data: ArrayLike
        :param background: True if data is for the X-ray background
        :type background: bool
        :param exposure_time: Observation exposure time (live).
        :type exposure_time: float
        """
        # TODO: Figure out consistent format for input data in docstring
        super().__init__(data)
        self.background = background
        self.exposure_time = exposure_time

    @classmethod
    def load_txt(
        cls, path: Path, background: bool, exposure_time: float, **kwargs
    ) -> Events:
        """
        Load events from a text file with x, y and ch columns.

        :param path: Path to file
        :type path: Path
        :param background: True if data is for the X-ray background
        :type background: bool
        :param exposure_time: Observation exposure time (live).
        :type exposure_time: float

        Additional keyword arguments are passed to `numpy.loadtxt`.
        """
        data = np.loadtxt(path, **kwargs)
        return cls(data, background, exposure_time)

    @classmethod
    def load_from_fits(
        cls,
        path: Path,
        background: bool,
        x_key: str = "X",
        y_key: str = "Y",
        channel_key: str = "PI",
        du_index=1,  # TODO: Confirm correct
    ) -> Events:
        """Load events from a fits file.
        The file will be converted to the text format used by BayesX.

        :param path: Path to events fits file
        :type path: Path
        :param x_key: Key of 'x' column in fits file, defaults to "X"
        :type x_key: str, optional
        :param y_key: Key of 'y' column in fits file, defaults to "Y"
        :type y_key: str, optional
        :param channel_key: Key of channel column in fits file, defaults to "PI"
        :type channel_key: str, optional
        :param du_index: List index of data unit in HDUList with events (0-indexed)
        :type du_index: int
        :param mode: `'evts'` for events, `bg` for background.
        :type mode: str
        """
        with fits.open(path) as fi:
            assert du_index < len(fi)

            f: PrimaryHDU = fi[du_index]  # type: ignore

            if f.header["extname"] != "EVENTS":
                log.warn(
                    "Trying to load events from a data unit that lacks events extension"
                )

            data = np.column_stack((f.data[x_key], f.data[y_key], f.data[channel_key]))  # type: ignore
            exposure_time: float = f.header["livetime"]  # type: ignore

            return cls(data, background, exposure_time)


class Mask(Data):
    def __init__(self, data: ArrayLike) -> None:
        super().__init__(data)
        assert self.data.ndim == 2

    @classmethod
    def load_from_npz(cls, path: Path):
        data = np.load(path)["arr_0"]
        return cls(data)

    @classmethod
    def load_from_reg(
        cls,
        path: Path,
        **kwargs,
    ) -> Mask:
        # TODO: Load using masking script
        raise NotImplementedError

    def bin(
        self,
        n_bins: int,
        cellsize: float,
        n_channels: int,
        outfile: Optional[Path] = None,
        **kwargs,
    ):
        kwargs["mask"] = True
        super().bin(n_bins, cellsize, outfile, n_channels=n_channels, **kwargs)


class ARF(Data):
    def __init__(self, data: ArrayLike) -> None:
        super().__init__(data)
        assert self.data.ndim == 1  # TODO: Better verification
        self.xrayNbins = len(self.data)  # TODO: Verify correctness

    @classmethod
    def load_from_txt(cls, path: Path, **kwargs) -> ARF:
        data = np.loadtxt(path)  # TODO: Update fornat for tabular text export

        return cls(data)

    @classmethod
    def load_from_fits(cls, path: Path) -> ARF:
        with fits.open(path) as f:
            data = f[1].data["specresp"]  # type: ignore

            # self.xrayEmin = np.min(f[1].data["energ_lo"])
            # self.xrayEmax = np.max(f[1].data["energ_hi"])

        return cls(data)

    def export(self, outfile: Path):
        self.path = outfile
        np.savetxt(outfile, self.data)

    def bin(self):
        raise NotImplementedError("ARF is not binnable")


class RMF(Data):
    def __init__(self, data: ArrayLike) -> None:
        super().__init__(data)
        assert self.data.ndim == 2  # TODO: Better verification
        self.xrayNbins, self.xrayNch = self.data.shape  # TODO: Verify correctness

    @classmethod
    def load_from_txt(cls, path: Path) -> RMF:
        data = np.loadtxt(path)  # TODO: Update fornat for tabular text export
        # reshaped = np.reshape(
        #     data, (nBins, n_channels), order="C"
        # )  # TODO: Verify correctness
        return cls(data)

    @classmethod
    def load_from_fits(cls, path: Path) -> RMF:
        with fits.open(path) as f:
            rmf = f[1].data["matrix"]  # type: ignore
            xrayNch = len(rmf[-1])
            xrayNbin = len(rmf)
            mat = np.zeros((xrayNbin, xrayNch))

            for i in range(0, len(rmf)):
                mat[i, : len(rmf[i])] = rmf[i]

        return cls(mat)

    def export(self, outfile: Path):
        self.path = outfile
        np.savetxt(outfile, np.ravel(self.data))

    def bin(self):
        raise NotImplementedError("RMF is not binnable")


@dataclass
class DataConfig:
    """Configuration options relevant to the input data."""

    # Input data
    filBG: Union[Path, Events]
    filevent: Union[Path, Events]
    filARF: Path  # in txt format
    filRMF: Path  # in text format

    nx: int  # Number of pixels in x direction
    ny: int  # Number of pixels in y direction
    xrayNbin: int  # Number of energy bins
    xrayNch: int  # Number of energy bins

    xraycell: float  # Spatial pixel size in arcsecond
    xrayEmin: float  # Minimum value of the energy range in keV
    xrayEmax: float  # Maximum value of the energy range in keV
    sexpotime: float  # Source exposure time in second
    bexpotime: float  # Background exposure time in second

    Aeffave: float = 250  # Average effective area of the telescope in cm^{2}

    filmask: Optional[Union[Path, Mask]] = None

    NHcol: float = 2.20e20  # Hydrogen column density in cm^2
    xrayBG_model: float = (
        8.4e-6  # Predicted background rate at each pixel in counts cm^-2 arcmin^-2s^-1,
    )

    rauto: bool = True  # Automatic radius calculation
    rmin: Optional[float] = None  # Minimum radius, Mpc
    rmax: Optional[float] = None  # Maximum radius for xray emission and GNFW model, Mpc
    rlimit: Optional[
        float
    ] = None  # Used to calculate logr, may need to be slightly higher than rmax, Mpc

    @classmethod
    def generate(
        cls,
        evts: Events,
        bg: Events,
        arf: ARF,
        rmf: RMF,
        out_path: Union[str, Path],
        bin_cellsize: int,
        nbins: int,
        energy_min: float,
        energy_max: float,
        mask: Optional[Mask] = None,
    ):
        # TODO: Automatic energy range

        out_path = Path(out_path)

        log.info("Binning events")
        evts.bin(nbins, bin_cellsize, out_path.joinpath("evts.txt"))
        log.info("Binning background")
        bg.bin(nbins, bin_cellsize, out_path.joinpath("bg.txt"))

        log.info(f"Events have dimensions ({evts.nx}, {evts.ny}, {evts.n_channels})")

        if (evts.nx, evts.ny, evts.n_channels) != (bg.nx, bg.ny, bg.n_channels):
            raise ValueError("Mismatched binned datasets.")

        mask_path = None
        if mask is not None:
            log.info("Binning mask")
            mask_path = out_path.joinpath("mask.txt")
            mask.bin(nbins, bin_cellsize, outfile=mask_path, n_channels=evts.n_channels)

        rmf.export(out_path.joinpath("rmf.txt"))
        arf.export(out_path.joinpath("arf.txt"))

        return cls(
            filBG=bg.path,  # type: ignore
            filevent=evts.path,  # type: ignore
            filARF=arf.path,
            filRMF=rmf.path,  # type: ignore
            nx=evts.nx,
            ny=evts.ny,
            xrayNbin=rmf.data.shape[0],
            xrayNch=rmf.data.shape[1],
            xraycell=bin_cellsize * 0.492,
            xrayEmin=energy_min,
            xrayEmax=energy_max,
            sexpotime=evts.exposure_time,
            bexpotime=bg.exposure_time,
            filmask=mask_path,
        )
