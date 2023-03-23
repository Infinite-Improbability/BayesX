from __future__ import annotations  # python 3.7 and up

from abc import ABC, abstractmethod
from dataclasses import dataclass
from logging import getLogger
from pathlib import Path
from typing import Optional, Union

import numpy as np
from astropy.io import fits
from astropy.io.fits.hdu import PrimaryHDU
from numpy.typing import ArrayLike
from pybayesx.binning import bin

log = getLogger(__name__)


class Data(ABC):
    def __init__(self, data: ArrayLike) -> None:
        self.data = np.array(data)
        self.path = None  # path to data in format ready for BayesX
        pass

    @classmethod
    def load_from_file(cls, path: Path, file_type: Optional[str] = None):
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
            return cls.load_from_txt(path)
        elif ext == ".fits":
            return cls.load_from_fits(path)
        else:
            raise ValueError("Unrecognised data file extension.")

    @classmethod
    @abstractmethod
    def load_from_txt(cls, path: Path) -> Data:
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def load_from_fits(cls, path: Path) -> Data:
        raise NotImplementedError

    def bin(
        self, n_bins: int, cellsize: float, outfile: Optional[Path] = None, **kwargs
    ):
        if outfile is not None:
            self.path = outfile
        return bin(
            self.data[:, 0],
            self.data[:, 1],
            self.data[:, 2],
            n_bins,
            cellsize,
            outfile=outfile,
            **kwargs,
        )

    def __str__(self) -> str:
        if self.path is None:
            raise ValueError("Data class has no file set. Please export to file.")
        return str(self.path)


class Events(Data):
    def __init__(self, data: ArrayLike, background: bool) -> None:
        super().__init__(data)
        assert self.data.ndim == 3
        self.background = background
        self.nx, self.ny, self.nChannels = self.data.shape  # TODO: Verify correctness

    @classmethod
    def load_from_txt(
        cls, path: Path, background: bool, nx: int, ny: int, nChannels: int, **kwargs
    ) -> Events:
        data = np.loadtxt(path)  # TODO: Update fornat for tabular text export
        reshaped = np.reshape(
            data, (nx, ny, nChannels), order="C"
        )  # TODO: Verify correctness
        """
        Used order C, where
            "read / write the elements using C-like index order,
            with the last axis index changing fastest, back to the first
            axis index changing slowest"
        based on the following from the Fortran code:

        DO xrayxpix = 1, xraynx
                DO xrayypix = 1, xrayny
                DO i = 1, xrayNch
                    xLENi = xLENi + 1
                    xrayCpred(xLENi) = xrayCmap(i, xrayxpix, xrayypix) + xrayBG(xLENi)
                END DO
                END DO
            END DO
        """
        return Events(reshaped, background)

    def load_from_fits(
        self,
        path: Path,
        background: bool,
        x_key: str = "X",
        y_key: str = "Y",
        channel_key: str = "PI",
        du_index=1,  # TODO: Confirm correct
        mode="evts",
        **kwargs,
    ) -> Events:
        """Load events from a fits file.
        The file will be converted to the text format used by BayesX.

        :param path: Path to events fits file
        :type path: Path
        :param x_key: _description_, defaults to "X"
        :type x_key: str, optional
        :param y_key: _description_, defaults to "Y"
        :type y_key: str, optional
        :param channel_key: _description_, defaults to "PI"
        :type channel_key: str, optional
        :param du_index: List index of data unit with events (0-indexed)
        :type du_index: int
        :param mode: `'evts'` for events, `bg` for background.
        :type mode: str
        """
        raise NotImplementedError  # needs verification

        with fits.open(path) as fi:
            assert du_index < len(fi)

            f: PrimaryHDU = fi[du_index]  # type: ignore

            if f.header["extname"] != "EVENTS":
                log.warn(
                    "Trying to load events from a data unit that lacks events extension"
                )

            if mode == "evts":
                self.events = np.column_stack((f[x_key], f[y_key], f[channel_key]))
                self.config["sexpotime"] = f.header["livetime"]
            elif mode == "bg":
                self.bg = np.column_stack((f[x_key], f[y_key], f[channel_key]))
                self.config["bexpotime"] = f.header["livetime"]


class Mask(Data):
    def __init__(self, data: ArrayLike, background: bool) -> None:
        super().__init__(data)
        assert self.data.ndim == 3
        self.background = background
        self.nx, self.ny, self.nChannels = self.data.shape  # TODO: Verify correctness

    @classmethod
    def load_from_txt(
        cls, path: Path, background: bool, nx: int, ny: int, nChannels: int, **kwargs
    ) -> Mask:
        data = np.loadtxt(path)
        reshaped = np.reshape(
            data, (nx, ny, nChannels), order="C"
        )  # TODO: Verify correctness
        return Mask(reshaped, background)

    def load_from_fits(
        self,
        path: Path,
        **kwargs,
    ) -> Mask:
        raise NotImplementedError

    def bin(
        self, n_bins: int, cellsize: float, outfile: Optional[Path] = None, **kwargs
    ):
        kwargs["mask"] = True
        super().bin(n_bins, cellsize, outfile, **kwargs)


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

    rmin: float = 0.01  # Minimum radius, Mpc
    rmax: float = 10  # Maximum radius for xray emission and GNFW model, Mpc
    rlimit: float = (
        10  # Used to calculate logr, may need to be slightly higher than rmax, Mpc
    )

    def load_all(self) -> tuple[Events, Events, Mask | None]:
        # This function is getting a little heavy on the error handling.
        # Can we make it more elegant?
        events_file = (
            self.filevent.path if isinstance(self.filevent, Events) else self.filevent
        )
        bg_file = self.filBG.path if isinstance(self.filBG, Events) else self.filBG

        if events_file is None or bg_file is None:
            raise ValueError(
                "Cannot load events or bg data: path to exported data not set."
            )

        events = Events.load_from_file(events_file)
        bg = Events.load_from_file(bg_file)

        mask = None
        if self.filmask is not None:
            mask_file = (
                self.filmask.path if isinstance(self.filmask, Mask) else self.filmask
            )
            if mask_file is None:
                raise ValueError(
                    "Cannot load mask data: path to exported data not set."
                )
            mask = Mask.load_from_file(mask_file)

        return events, bg, mask  # type: ignore

    def bin_all(self, nbins: int, cellsize: int, outdir: Path, **kwargs):
        (
            events,
            bg,
            mask,
        ) = self.load_all()  # should we load and discard each to lower memory use?

        filevent = outdir.joinpath(f"events_binned{nbins}of{cellsize}.txt")
        events.bin(
            nbins,
            cellsize,
            outfile=filevent,
            **kwargs,
        )
        self.filevent = (
            filevent  # run after so it doesn't get changed if the binning fails
        )

        filBG = outdir.joinpath(f"background_binned{nbins}of{cellsize}.txt")
        bg.bin(
            nbins,
            cellsize,
            outfile=filBG,
            **kwargs,
        )
        self.filBG = filBG

        filmask = outdir.joinpath(f"mask_binned{nbins}of{cellsize}.txt")
        if mask:
            mask.bin(
                nbins,
                cellsize,
                outfile=filmask,
                **kwargs,
            )
        self.filmask = filmask

        self.xraycell *= cellsize
        self.nx = nbins
        self.ny = nbins
