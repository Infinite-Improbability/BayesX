from __future__ import annotations  # python 3.7 and up

from dataclasses import dataclass
from logging import getLogger
from os import makedirs
from pathlib import Path
from typing import Optional, Sequence, Union

import numpy as np
from astropy.io import fits
from astropy.io.fits.hdu import PrimaryHDU
from numpy.typing import ArrayLike, NDArray
from scipy.stats import binned_statistic_dd

from .mask import mask

log = getLogger(__name__)


class Data:
    """Abstract base class for source data"""

    def __init__(self, data: ArrayLike) -> None:
        """Basic initialiser

        :param data: Array on input data
        :type data: ArrayLike
        """
        self.data = np.array(data)
        self.path = None  # path to data in format ready for BayesX

    def __str__(self) -> str:
        if self.path is None:
            raise ValueError("Data class has no file set. Please export to file.")
        return str(self.path)


class BinnableData(Data):
    def __init__(self, data: ArrayLike) -> None:
        super().__init__(data)

    def bin(
        self,
        cellsize: float,
        outfile: Optional[Union[Path, str]] = None,
        **kwargs,
    ):
        """Bin arbitary data with three dimensions in first two dimensions.

        :param n_bins: Amount of bins in one dimension
        :type n_bins: int
        :param cellsize: Size of a single bin, in units of source data
        :type cellsize: float
        :param outfile: Path to write binned data to, defaults to None
        :type outfile: Path | str, optional
        :return: 1D array of binned data
        :rtype: np.ndarray[Any, np.dtype[np.float64]]
        """
        if outfile is not None:
            outfile = Path(outfile)
            self.path = outfile

        binned, new_edges = self._bin(
            self.data[:, 0],
            self.data[:, 1],
            self.data[:, 2],
            cellsize,
            outfile=outfile,
            **kwargs,
        )
        self.binned = binned
        self.nx, self.ny, self.n_channels = binned.shape
        return binned, new_edges

    def _bin(
        self,
        x: ArrayLike,
        y: ArrayLike,
        channels: ArrayLike,
        cellsize: float,
        edges: Optional[Sequence[ArrayLike]] = None,
        n_bins: Optional[tuple[int, int]] = None,
        origin: Optional[tuple[float, float]] = None,
        n_channels: Optional[int] = None,
        outfile: Optional[Path] = None,
        mask: bool = False,
    ) -> tuple[NDArray, list[NDArray]]:
        """Bin data given spatial coordinates and channel coordinates of events.
        Bins are only spatial and do not reduce channel count.
        This will crop points if they are not within `n_bins/2` bins of the center.
        If `mask = True` then instead bin passed on mask status.

        The current implementation of this function leaves a cross artefact on the
        binned data.

        :param x: 1D array of x coordinates of events. Should have an entry for every
        point.
        :type x: numpy.typing.ArrayLike
        :param y: 1D array of y coordinates of events, assumed to use the same units and
        order as x.
        :type y: numpy.typing.ArrayLike
        :param channel: 1D array of channels of events, assumed to use the same order
         as x. If `mask` is True then this value is ignored.
        :type channel: numpy.typing.ArrayLike
        :param cellsize: The size of a bin, in units matching x and y
        :type cellsize: float
        :param nbins: The number of bins along each spatial axis, defaults to None. If excluded
         bin number is determined by `floor((max-min) / cellsize)` on each axis.
        :type nbins: [int, int], optional
        :param origin: Central (origin) point (x0, y0) along in data, defaults to None. If
         None then the midpoint is determined from the range of the data.
        :type origin: [float, float], optional
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
        # TODO: Rewrite as class method

        # Coerce input
        x = np.asarray(x)
        y = np.asarray(y)
        channels = np.asarray(channels)

        # Set no. of channels based on largest channel number present
        if mask and (n_channels is None or n_channels < 1):
            raise ValueError("If binning a mask, max_chan is a required argument.")
        elif n_channels is None:
            n_channels = int(np.ptp(channels)) + 1  # max - min

        if mask:
            log.debug(f"n_channels is given as {n_channels}")
            channels = np.linspace(0, n_channels, n_channels)
            channels = np.tile(channels, x.size)
            x = x.repeat(n_channels)
            y = y.repeat(n_channels)

        # Do some input checking
        # We do this now to ensure mask modifications are included
        if len(x) != len(y) != len(channels):
            raise ValueError("Coordinate lists x and y have different lengths.")
        if x.ndim != 1 or y.ndim != 1 or channels.ndim != 1:
            raise ValueError("An input array is not one dimensional.")

        # If origins are not provided default to centre of data
        if origin is None:
            x0 = (np.max(x) + np.min(x)) / 2
            y0 = (np.max(y) + np.min(y)) / 2
        else:
            x0 = origin[0]
            y0 = origin[1]

        # Get edges of data
        x_min = np.min(x)
        y_min = np.min(y)
        x_max = np.max(x)
        y_max = np.max(y)
        ch_min = np.min(channels)
        ch_max = np.max(channels)

        # Calculate the number of bins.
        # We don't want partial bins so we use integer division.
        n_x_bins = int((x_max - x_min) // cellsize)
        n_y_bins = int((y_max - y_min) // cellsize)
        n_ch_bins = int(n_channels) + 1

        # If the user has provided bins check their input against our own
        # Then do what they ask because sometimes it is valid.
        if n_bins:
            if n_x_bins < n_bins[0] or n_y_bins < n_bins[1]:
                log.warning(
                    "Requesting more bins than the data spans. This may be due to bad input, or there may simply be no events at the edge of the observation."
                )
            n_x_bins, n_y_bins = n_bins

        if (x_max - x_min) % cellsize != 0:
            log.info("Trimmed x data to have whole number of bins")
            x_min = x0 - (cellsize * n_x_bins / 2)
            x_max = x0 + (cellsize * n_x_bins / 2)
        if (y_max - y_min) % cellsize != 0:
            log.info("Trimmed y data to have whole number of bins")
            y_min = y0 - (cellsize * n_y_bins / 2)
            y_max = y0 + (cellsize * n_y_bins / 2)

        # It would be easier to just get scipy to do all the edge generation for us but
        # this gives us more control and repeatibility
        if edges is not None:
            x_edges, y_edges, _ = edges
        else:
            x_edges = np.linspace(x_min, x_max, n_x_bins)
            y_edges = np.linspace(y_min, y_max, n_y_bins)
        ch_edges = np.linspace(ch_min, ch_max, n_ch_bins)

        data = np.column_stack((x, y, channels))

        log.debug(f"Binning on range x={x_edges.min()}:{x_edges.max()}, y={y_edges.min()}:{y_edges.max()}, ch={ch_edges.min()}:{ch_edges.max()}")  # type: ignore

        counts: NDArray
        new_edges: list[NDArray]
        counts, new_edges, bin_numbers = binned_statistic_dd(
            data,
            values=None,
            statistic="count",
            bins=[x_edges, y_edges, ch_edges],  # type: ignore
        )

        # log.debug(f"Binning done with edges {edges}")

        counts_1d = counts.ravel()

        if outfile:
            np.savetxt(outfile, counts_1d, "%5d")

        return counts, new_edges

    def bin_plot(self, channel_index=None):
        if not hasattr(self, "binned"):
            raise ValueError("Data not binned")

        try:
            import matplotlib.pyplot as plt
            from matplotlib.colors import LogNorm

            getLogger("matplotlib").setLevel("WARNING")
            if channel_index is None:
                # Sum over counts
                binned = np.sum(self.binned, 2)
            else:
                binned = self.binned[:, :, channel_index]

            plt.imshow(binned, norm=LogNorm(), origin="lower", interpolation="none")
            plt.colorbar()
        except ImportError:
            log.warn("Matplotlib not installed")
            return


class Events(BinnableData):
    def __init__(
        self,
        data: ArrayLike,
        exposure_time: float,
        background: bool,
        energy_range: Optional[Sequence[float]] = None,
    ) -> None:
        """X-ray event data

        :param data: Source data array. Three columns, with each row containing (x, y, channel).
        :type data: ArrayLike
        :param exposure_time: Observation exposure time (live).
        :type exposure_time: float
        :param background: True if data is for the X-ray background
        :type background: bool
        :param energy_range: Range of input energies, lower limit followed by upper limit in keV, defaults to None
        :type energy_range: Sequence[float], optional
        """
        # TODO: Figure out consistent format for input data in docstring
        super().__init__(data)
        self.background = background
        self.exposure_time = exposure_time
        self.energy_range = energy_range

    @classmethod
    def load_txt(
        cls,
        path: Union[Path, str],
        exposure_time: float,
        background: bool = False,
        **kwargs,
    ) -> Events:
        """
        Load events from a text file with x, y and ch columns.

        :param path: Path to file
        :type path: Path | str
        :param background: True if data is for the X-ray background, defaults to False.
        :type background: bool
        :param exposure_time: Observation exposure time (live).
        :type exposure_time: float

        Additional keyword arguments are passed to `numpy.loadtxt`.
        """
        path = Path(path)
        data = np.loadtxt(path, **kwargs)
        return cls(data, exposure_time, background)

    @classmethod
    def load_fits(
        cls,
        path: Union[Path, str],
        background: bool = False,
        x_key: str = "X",
        y_key: str = "Y",
        channel_key: str = "PI",
        energy_key: str = "energy",
        du_index=1,  # TODO: Confirm correct
    ) -> Events:
        """Load events from a fits file.
        The file will be converted to the text format used by BayesX.

        :param path: Path to events fits file
        :type path: Path | str
        :param background: True if data is for the X-ray background, defaults to False.
        :type background: bool
        :param x_key: Key of 'x' column (assumes pixels) in fits file, defaults to "X"
        :type x_key: str, optional
        :param y_key: Key of 'y' column (assumes pixels) in fits file, defaults to "Y"
        :type y_key: str, optional
        :param channel_key: Key of channel column in fits file, defaults to "PI"
        :type channel_key: str, optional
        :param energy_key: Key of energy column (assumes eV) in fits file, defaults to "energy"
        :type energy_key: str, optional
        :param du_index: List index of data unit in HDUList with events (0-indexed)
        :type du_index: int
        :param mode: `'evts'` for events, `bg` for background.
        :type mode: str
        """
        path = Path(path)

        with fits.open(path) as fi:
            assert du_index < len(fi)

            f: PrimaryHDU = fi[du_index]  # type: ignore

            if f.header["extname"] != "EVENTS":
                log.warn(
                    "Trying to load events from a data unit that lacks events extension"
                )

            data = np.column_stack((f.data[x_key], f.data[y_key], f.data[channel_key]))  # type: ignore
            exposure_time: float = f.header["livetime"]  # type: ignore

            # Assuming eV
            energy_min = np.min(f.data[energy_key])  # type: ignore
            energy_max = np.max(f.data[energy_key])  # type: ignore

            log.info(
                f"Detected energy range {energy_min / 1000}:{energy_max / 1000} keV"
            )

            # Convert from eV to 10eV then convert to keV
            energy_min = np.floor(energy_min / 10) / 100
            energy_max = np.ceil(energy_max / 10) / 100

            log.info(f"Rounded energy range to {energy_min}:{energy_max} keV")
            log.info("Events loading complete")

            return cls(data, exposure_time, background, (energy_min, energy_max))


class Mask(BinnableData):
    def __init__(self, data: ArrayLike) -> None:
        super().__init__(data)
        assert self.data.ndim == 2

    @classmethod
    def load_npz(cls, path: Union[Path, str]):
        """Load mask as boolean array from numpy export.

        :param path: Path to array file
        :type path: Path | str
        :return: Returns a new Mask object
        :rtype: Mask
        """
        path = Path(path)
        data = np.load(path)["arr_0"]
        return cls(data)

    @classmethod
    def load_reg(
        cls,
        path: Union[Path, str],
        xMin: float,
        xMax: float,
        yMin: float,
        yMax: float,
    ):
        """Load mask as boolean array from ds9 reg. Only supports ellipses at present

        :param path: Path to array file
        :type path: Union[Path, str]
        :return: Returns a new Mask object
        :rtype: Mask

        For other parameters see `Mask.mask()`
        """
        path = Path(path)
        data = mask(xMin, xMax, yMin, yMax, [path])

        log.info("Region loading complete")

        return cls(data)

    def bin(
        self,
        cellsize: float,
        n_channels: int,
        outfile: Optional[Union[Path, str]] = None,
        **kwargs,
    ):
        kwargs["mask"] = True
        super().bin(cellsize, outfile, n_channels=n_channels, **kwargs)


class ARF(Data):
    def __init__(
        self,
        data: ArrayLike,
    ) -> None:
        super().__init__(data)
        assert self.data.ndim == 1  # TODO: Better verification
        self.xrayNbins = len(self.data)  # TODO: Verify correctness

    @classmethod
    def load_txt(cls, path: Union[Path, str], **kwargs) -> ARF:
        path = Path(path)
        data = np.loadtxt(path)  # TODO: Update fornat for tabular text export

        return cls(data)

    @classmethod
    def load_fits(cls, path: Union[Path, str]) -> ARF:
        path = Path(path)
        with fits.open(path) as f:
            data = f[1].data["specresp"]  # type: ignore

        log.info("ARF loading complete")

        return cls(data)

    def export(self, outfile: Union[Path, str]):
        outfile = Path(outfile)
        self.path = outfile

        np.savetxt(outfile, self.data)


class RMF(Data):
    def __init__(self, data: ArrayLike) -> None:
        super().__init__(data)
        assert self.data.ndim == 2  # TODO: Better verification
        self.xrayNbins, self.xrayNch = self.data.shape  # TODO: Verify correctness
        self.energy_bounds: Optional[NDArray] = None

    @classmethod
    def load_txt(cls, path: Union[Path, str]) -> RMF:
        path = Path(path)
        data = np.loadtxt(path)  # TODO: Update fornat for tabular text export
        # reshaped = np.reshape(
        #     data, (nBins, n_channels), order="C"
        # )  # TODO: Verify correctness
        return cls(data)

    @classmethod
    def load_fits(cls, path: Union[Path, str]) -> RMF:
        path = Path(path)

        with fits.open(path) as f:
            data = f[1].data  # type: ignore

            rmf = data["MATRIX"]

            first_channel = np.concatenate(data["F_CHAN"])
            n_channels = np.concatenate(data["N_CHAN"])
            last_channel = first_channel + n_channels - 1

            xrayNch = np.max(last_channel)
            xrayNbin = len(rmf)

            log.debug(f"RMF has {xrayNbin} bins and {xrayNch} channels")

            mat = np.zeros((xrayNbin, xrayNch))

            for i in range(0, xrayNbin):
                # correct for zero indexing
                first_index = first_channel[i] - 1
                # We don't need zero indexing correction because it's cancelled out by Python's exclusive slice endpoint
                last_index = last_channel[i]

                mat[i, first_index:last_index] = rmf[i]

            new = cls(mat)
            new.energy_bounds = f[2].data  # type: ignore

        log.info("RMF loading complete")

        return new

    def export(
        self, outfile: Union[Path, str], energy_range: Optional[Sequence[float]] = None
    ):
        outfile = Path(outfile)
        self.path = outfile

        data = self.data

        # Trim bins
        min_index = 0
        max_index = self.data.size
        if energy_range is not None:
            if self.energy_bounds is not None:
                # Very unoptimised

                for i, e in enumerate(self.energy_bounds["E_MAX"]):
                    # if e_max for channel greater than specified minimum energy
                    if e > energy_range[0]:
                        # We are assuming Ebounds order matches channel order, ascending
                        min_index = i
                        log.debug(f"Min index found ({i}) with max energy of {e}")
                        break
                else:
                    log.warn("Unable to find minimum index using energy bounds")

                for i, e in enumerate(self.energy_bounds["E_MIN"]):
                    # if e_min for channel greater than specified maximum energy
                    if e > energy_range[1]:
                        max_index = i
                        log.debug(f"Max index found ({i}) with min energy of {e}")
                        break
                else:
                    log.warn("Unable to find maximum index using energy bounds")

                if min_index >= max_index:
                    raise ValueError("Energy bounding got inverted bounds.")

                data = self.data[:, min_index:max_index]
            else:
                log.warn("Unable to constrain energy range on RMF")

        self.n_bins = data.shape[0]
        self.n_channels = data.shape[1]
        np.savetxt(outfile, np.ravel(data))


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

    comments: Optional[str] = None  # User comments

    @classmethod
    def generate(
        cls,
        evts: Events,
        bg: Events,
        arf: ARF,
        rmf: RMF,
        out_path: Union[str, Path],
        bin_size: int,
        energy_range: Optional[Sequence[float]] = None,
        mask: Optional[Mask] = None,
        cell_size: float = 0.492,
    ):
        """Generate DataConfig from a collection of Data objects

        :param evts: Source events
        :type evts: Events
        :param bg: Background events
        :type bg: Events
        :param arf: Ancillary Response File
        :type arf: ARF
        :param rmf: Response/Redistribution Matrix File
        :type rmf: RMF
        :param out_path: Folder in which to save exported data
        :type out_path: Union[str, Path]
        :param bin_size: Size of bin (side length) in same units as spatial coordinates
        :type bin_size: int
        :param nbins: Number of bins kept (side length)
        :type nbins: int
        :param energy_range: Range of input energies, lower limit followed by upper limit in keV, defaults to None
        :type energy_range: Sequence[float], optional
        :param mask: Mask file of regions, defaults to None
        :type mask: Mask, optional
        :param cell_size: Size (side length) of pixel in arcseconds before binning, defaults to 0.492 (as for Chandra ACIS)
        :type cell_size: float, optional
        :raises ValueError: Raised if events and background have different shapes after binning.
        :return: A new DataConfig object
        :rtype: DataConfig
        """

        out_path = Path(out_path)

        makedirs(out_path, exist_ok=True)

        # Set energy range from source data
        if energy_range is None:
            if evts.energy_range is not None:
                energy_range = evts.energy_range
                if evts.energy_range != bg.energy_range:
                    log.warn(
                        "Background and source energy ranges don't match, using source range."
                    )
            elif bg.energy_range is not None:
                energy_range = bg.energy_range
            else:
                raise ValueError("Unable to detect energy range, please specify.")

        # Bin
        log.info("Binning source events...")
        _, edges = evts.bin(bin_size, out_path.joinpath("evts.txt"))
        log.info("Binning background events...")
        bg.bin(bin_size, out_path.joinpath("bg.txt"), edges=edges)

        # Validate binning before proceeding
        if (evts.nx, evts.ny, evts.n_channels) != (bg.nx, bg.ny, bg.n_channels):
            log.info(
                f"Source events have dimensions ({evts.nx}, {evts.ny}, {evts.n_channels})"
            )
            log.info(
                f"Background events have dimensions ({bg.nx}, {bg.ny}, {bg.n_channels})"
            )
            raise ValueError(
                "Source and background datasets have different shapes after binning."
            )
        else:
            log.info(
                f"Events have dimensions ({evts.nx}, {evts.ny}, {evts.n_channels})"
            )

        # Bin mask
        mask_path = None
        if mask is not None:
            log.info("Binning mask...")
            mask_path = out_path.joinpath("mask.txt")
            mask.bin(
                bin_size, outfile=mask_path, n_channels=evts.n_channels, edges=edges
            )

            log.info(f"Mask has dimensions ({mask.nx}, {mask.ny}, {mask.n_channels})")

            # Validate binning before proceeding
            if (evts.nx, evts.ny, evts.n_channels) != (
                mask.nx,
                mask.ny,
                mask.n_channels,
            ):
                raise ValueError(
                    f"Source and mask datasets have different shapes after binning."
                )

        # Export RMF and ARF
        rmf.export(out_path.joinpath("rmf.txt"), energy_range=energy_range)
        arf.export(out_path.joinpath("arf.txt"))

        log.info(f"ARF has {arf.data.size} energy bins.")
        log.info(
            f"RMF has {rmf.n_bins} energy bins and {rmf.n_channels} energy channels (in export)."
        )

        # More data validation
        if arf.data.size != rmf.n_bins:
            raise ValueError(
                f"ARF and RMF have different numbers of energy bins ({arf.data.size}, {rmf.n_bins})"
            )
        elif rmf.n_channels != evts.n_channels:
            raise ValueError(
                f"RMF and source data have different numbers of energy channels ({rmf.n_channels}, {evts.n_channels})"
            )

        return cls(
            filBG=bg.path,
            filevent=evts.path,
            filARF=arf.path,
            filRMF=rmf.path,
            nx=evts.nx,
            ny=evts.ny,
            xrayNbin=rmf.n_bins,
            xrayNch=rmf.n_channels,
            xraycell=bin_size * cell_size,
            xrayEmin=energy_range[0],
            xrayEmax=energy_range[1],
            sexpotime=evts.exposure_time,
            bexpotime=bg.exposure_time,
            filmask=mask_path,
        )
