import numpy as np
from astropy.constants import k_B, h
from scipy.integrate import simpson
from ntsz_utils import get_planck_bandpass, get_data_path
from collections import deque
from ntsz import ntsz_single_interpol, ntsz_broken_interpol
import os


class BandpassCorrector:
    """
    Compute bandpass-corrected (non-thermal / kinetic) Sunyaev-Zeldovich (SZ)
    spectra for Planck frequency channels.

    Parameters
    ----------
    path_to_grid : str, optional
        Path to the ntSZ spectra grid FITS file. Defaults to the data folder
        provided with this package if none specified.
    alpha : float, optional
        Index for the single power-law describing the non-thermal
        electrons. Default is 3.61.
    alpha2 : float, optional
        Index for the high-energy part of the broken power-law.
        Default is 3.61.
    pmin : float, optional
        Minimum momentum for non-thermal electrons. Default is 1.0.
    pbreak : float, optional
        Break momentum for the broken model. Required if using the broken
        power-law momentum distribution for the non-thermal electrons.

    Notes
    -----
    Public Methods
    --------------
    - compute(model="single"):
        Compute bandpass-corrected SZ spectra for the chosen model. One can
        choose between "single", "broken", or "ksz" if choosing to compute
        ntSZ from single or broken power-law models for non-thermal electrons
        or kSZ effect.

    Private Methods
    ---------------
    - _validate_grid_path():
        Checks the validity of the provided grid path.
    - _bandpass_weight_single():
        Calculates weights for the single power-law model.
    - _bandpass_weight_broken():
        Calculates weights for the broken power-law model.
    - _bandpass_weight_ksz():
        Calculates weights for the kSZ model.
    - _get_band_data():
        Collects Planck band frequencies and transmission data from IMO.
    """

    def __init__(
        self, path_to_grid=None, alpha=3.61, alpha2=3.61, pmin=1.0,
        pbreak=None, model='single'
            ):
        self.T_CMB = 2.7255  # K
        self.alpha = alpha
        self.alpha2 = alpha2
        self.pmin = pmin
        self.pbreak = pbreak
        self.path_to_grid = path_to_grid
        self.model = model
        self.planck_freq = {
            "LFI": [70],
            "HFI": [100, 143, 217, 353, 545, 857]
        }
        if self.path_to_grid is None:
            self.path_to_grid = get_data_path(
                subfolder="data"
            )/f"ntsz_{model}_grid.fits"
        # Validate grid path
        self._validate_grid_path()

    def _validate_grid_path(self):
        if os.path.isfile(self.path_to_grid):
            print(f"Using {self.path_to_grid} as the ntSZ spectra grid\
                 for interpolation.")
        elif os.path.isdir(self.path_to_grid):
            raise IsADirectoryError(
                f"Expected a FITS file for the ntSZ grid, but found a\
                    directory at: {self.path_to_grid}\n"
                f"Please provide the full path to the FITS file."
            )
        else:
            raise FileNotFoundError(
                f"Missing grid FITS file or directory at:\
                    {self.path_to_grid}\n"
                f"Expected something like '<path-to-dir>/ntsz_grid.fits'"
            )

    @staticmethod
    def _bandpass_weight_single(
        freq, freq_nominal, transmission, pmin, path_to_grid, alpha
            ):
        h_x = ntsz_single_interpol(pmin, freq, path_to_grid, alpha)
        weight = simpson(h_x * transmission, freq) / (
            simpson(transmission * (freq_nominal / freq / 1e9), freq)
        )
        return weight

    @staticmethod
    def _bandpass_weight_broken(
        freq, freq_nominal, transmission, pbreak, path_to_grid, pmin, alpha2
            ):
        h_x = ntsz_broken_interpol(pbreak, freq, path_to_grid, pmin, alpha2)
        weight = simpson(h_x * transmission, freq) / (
            simpson(transmission * (freq_nominal / freq / 1e9), freq)
        )
        return weight

    def _bandpass_weight_ksz(self, freq, freq_nominal, transmission):
        x = freq * 1e9 * h / k_B / self.T_CMB
        h_x = (x ** 4) * np.exp(x) / (np.exp(x) - 1) ** 2
        weight = simpson(h_x * transmission, freq) / (
            simpson(transmission * (freq_nominal / freq / 1e9), freq)
        )
        return weight

    def _get_band_data(self):
        """Collect Planck band frequencies and transmissions."""
        frequencies = deque()
        transmissions = deque()

        for k, chans in self.planck_freq.items():
            for chan in chans:
                freq_ch, trans_ch = get_planck_bandpass(frequency=chan)
                frequencies.append(freq_ch)
                transmissions.append(trans_ch)

        freq_nominal = np.array((70, 100, 143, 217, 353, 545, 857))
        return np.array(frequencies, dtype=object), \
            np.array(transmissions, dtype=object), freq_nominal

    def compute(self, model="single"):
        """Compute bandpass-corrected SZ weights for the specified model."""
        frequencies, transmissions, freq_nominal = self._get_band_data()
        est_weights = np.zeros(len(frequencies))

        if model == "single":
            for i in range(len(frequencies)):
                est_weights[i] = self._bandpass_weight_single(
                    np.array(frequencies[i])[1:],
                    freq_nominal[i] * 1e9,
                    transmissions[i][1:],
                    self.pmin,
                    self.path_to_grid,
                    self.alpha
                )

        elif model == "broken":
            if self.pbreak is None:
                raise ValueError("pbreak must be set for 'broken' model.")
            for i in range(len(frequencies)):
                est_weights[i] = self._bandpass_weight_broken(
                    np.array(frequencies[i])[1:],
                    freq_nominal[i] * 1e9,
                    transmissions[i][1:],
                    self.pbreak,
                    self.path_to_grid,
                    self.pmin,
                    self.alpha2
                )

        elif model == "ksz":
            for i in range(len(frequencies)):
                est_weights[i] = self._bandpass_weight_ksz(
                    np.array(frequencies[i])[1:],
                    freq_nominal[i] * 1e9,
                    transmissions[i][1:]
                )

        else:
            raise ValueError(f"Unknown model '{model}'. \
            Choose from 'single', 'broken', or 'ksz'.")

        return est_weights


if __name__ == "__main__":
    corrector = BandpassCorrector()
    weights = corrector.compute(model="single")
    print("Computed bandpass-corrected SZ weights:", weights)
