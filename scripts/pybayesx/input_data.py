from __future__ import annotations

import logging  # Python 3.7 and up
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Optional

import numpy as np
from astropy.io import fits
from astropy.io.fits.hdu import PrimaryHDU
from numpy.typing import ArrayLike

from .binning import bin

log = logging.getLogger(__name__)


class Data:
    def __init__(self, data: ArrayLike) -> None:
        self.data = np.array(data)
        pass

    @classmethod
    def load_from_file(cls, path: Path, file_type: Optional[str] = None) -> Data:
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
        return NotImplementedError

    @classmethod
    @abstractmethod
    def load_from_fits(cls, path: Path) -> Data:
        return NotImplementedError

    @abstractmethod
    def bin_and_export(
        self, n_bins: int, cellsize: float, outfile: Path, n_channels: int, **kwargs
    ):
        pass


class Events(Data):
    def __init__(self, data: ArrayLike, background: bool) -> None:
        super().__init__(data)
        assert self.data.ndim == 3
        self.background = background
        self.nx, self.ny, self.nChannels = self.data.shape  # TODO: Verify correctness

    @classmethod
    def load_from_txt(
        cls, path: Path, background: bool, nx: int, ny: int, nChannels: int, **kwargs
    ) -> Data:
        data = np.loadtxt(path)
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

    def load_events_from_fits(
        self,
        path: Path,
        background: bool,
        x_key: str = "X",
        y_key: str = "Y",
        channel_key: str = "PI",
        du_index=1,  # TODO: Confirm correct
        mode="evts",
        **kwargs,
    ) -> Data:
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
        return NotImplementedError  # type: ignore

        with fits.open(path) as fi:
            assert du_index < len(fi) + 1

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
