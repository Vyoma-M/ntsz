"""
The classes in this module can be used to generate grids of ntSZ spectra
where the axes can be imagined to be the observing frequency and the pmin
(pbreak) in case of a power-law (broken power-law) describing the momentum
distribution of the scattering non-thermal electrons.
"""
import numpy as np
from astropy.io import fits
from astropy.constants import pc, e, m_e, k_B, c, sigma_T, h
from scipy.integrate import quad
from multiprocessing import Pool
from ntsz.ntsz import photon_scat_kernel, i_x


class NTSZGridGenerator_single:
    """
    Compute bandpass-corrected (non-thermal / kinetic) Sunyaev-Zeldovich (SZ)
    spectra for Planck frequency channels.

    Parameters
    ----------
    alpha : float
        Index for the single power-law describing the non-thermal
        electrons. Default is 3.61.
    nproc : int
        Number of processes for parallel computation
    filename : str
        Name of file along with path for where the grid needs to be saved.
    pmax : float, optional
        Maximum momentum for non-thermal electrons. Default is 5e5 (5*10^5).

    Notes
    -----
    Public Methods
    --------------
    - create_grid():
        Compute ntSZ spectra and save them into FITS file with <filename>. If
        method is assigned to a variable, variable will also contain the data.

    Private Methods
    ---------------
    - _powerlaw_int():
        Compute ntSZ spectrum for power-law distribution of non-thermal
        electrons.
    - _save_grid():
        Utility function that saves computed grid into a FITS file with
        user-defined path_to_grid='<path-to-file->/filename'.
    """
    def __init__(
        self, alpha=3.61, nproc=32, filename='ntsz_grid_single.fits', pmax=5e5
            ):
        self.alpha = alpha
        self.nproc = nproc
        self.filename = filename
        self.pmax = pmax

        # Physical constants
        self.pc = pc.si.value
        self.sigma_T = sigma_T.si.value
        self.c = c.si.value
        self.m_e = m_e.si.value
        self.e = e.si.value
        self.h = h.si.value
        self.k_B = k_B.si.value
        self.T_CMB = 2.7255  # K

        # Frequency setup
        self.frequencies = np.geomspace(1e10, 2.8e12, 1000)
        self.xf = (self.h * self.frequencies) / (self.k_B * self.T_CMB)
        self.i0 = 2 * (self.k_B * self.T_CMB) ** 3 / (self.h * self.c) ** 2\
                    * 1e20

    def _powerlaw_int(self, pmin):
        alpha = self.alpha
        pmax = self.pmax
        f = self.frequencies
        xx = (self.h * f) / (self.k_B * self.T_CMB)
        spec = np.zeros(len(xx))
        norm = (alpha - 1) / (pmin ** (1 - alpha) - pmax ** (1 - alpha))

        def integrand(p, alpha, x):
            return photon_scat_kernel(p, x) * p ** (-alpha)

        for i in range(len(xx)):
            spec[i] = quad(integrand, pmin, pmax, args=(alpha, xx[i]))[0]

        return norm * spec

    def _save_grid(self, grid_data):
        hdr = fits.Header()
        hdr['COMMENT'] = 'This is a file containing non-thermal SZ spectra due'
        hdr['COMMENT'] = 'to non-thermal population of electrons with momenta'
        hdr['COMMENT'] = 'distribution corresponding to a power-law as'
        hdr['COMMENT'] = 'a function of pmin and observing frequency. It is '
        hdr['COMMENT'] = 'a grid where the axes are pmin and frequency in Hz'
        hdr['COMMENT'] = 'and each point corresponds to a value of the corresponding'
        hdr['COMMENT'] = 'distortion in specific intensity due to ntSZ.'
        primary_hdu = fits.PrimaryHDU(header=hdr)

        hdr = fits.Header()
        hdr['NAME'] = 'ntSZ spectra'
        hdr['COMMENT'] = f'alpha = {self.alpha}'
        hdr['COMMENT'] = f'pmin = {self.pmin}'
        hdr['COMMENT'] = 'pmin = np.linspace(-1, +3.0, 200)'
        hdr['COMMENT'] = 'gridf = np.geomspace(1e10, 2.8e12, 1000)'
        hdr['COMMENT'] = f'pmax = {self.pmax}'
        col = fits.ImageHDU(data=grid_data, header=hdr)
        hdul = fits.HDUList([primary_hdu, col])
        hdul.writeto(self.filename, overwrite=True)

    def create_grid(self):
        pm1 = np.linspace(-1, 3.0, 200)
        pm = 10 ** pm1

        print("Starting integration to create a grid...")
        pool = Pool(processes=self.nproc)
        data = pool.map(self._powerlaw_int, [pmin for pmin in pm])
        pool.close()
        parr = np.array(data)
        parr = (parr - i_x(self.xf)) * self.i0

        self._save_grid(self.filename, parr)
        print(f"Grid saved to {self.filename}")

        return parr


class NTSZGridGenerator_broken:
    """
    Compute ntSZ on a grid with the observing frequency and pbreak as the
    axes. Such a grid can be used for fast computation of ntSZ spectra
    for a broken power-law of momenta distribution for non-thermal electrons
    as a function of break momentum and observing frequency while keeping
    other variables of the electron distribution fixed.

    Parameters
    ----------
    alpha : float
        Index for the flat part of broken power-law describing the non-thermal
        electrons. Default is 0.05. User is recommended to use >0. for
        computational purposes.
    alpha2 : float
        Index for the single power-law describing the non-thermal
        electrons. Default is 3.61.
    pmin : float, optional
        Minimum momentum for non-thermal electrons. Default is 1.0.
    nproc : int
        Number of processes for parallel computation
    filename : str
        Name of file along with path for where the grid needs to be saved.
    pmax : float, optional
        Maximum momentum for non-thermal electrons. Default is 5e5 (5*10^5).

    Notes
    -----
    Public Methods
    --------------
    - create_grid():
        Compute ntSZ spectra and save them into FITS file with <filename>. If
        method is assigned to a variable, variable will also contain the data.

    Private Methods
    ---------------
    - _broken_powerlaw_int():
        Compute ntSZ spectrum for broken power-law distribution of non-thermal
        electrons.
    - _save_grid():
        Utility function that saves computed grid into a FITS file with
        user-defined path_to_grid='<path-to-file->/filename'.
    """
    def __init__(
        self, alpha=0.05, pmin=1.0, alpha2=3.61, nproc=32,
        filename='ntsz_grid_broken.fits', pmax=5e5
            ):
        self.alpha = alpha
        self.alpha2 = alpha2
        self.pmin = pmin
        self.nproc = nproc
        self.filename = filename
        self.pmax = pmax

        # Physical constants
        self.pc = pc.si.value
        self.sigma_T = sigma_T.si.value
        self.c = c.si.value
        self.m_e = m_e.si.value
        self.e = e.si.value
        self.h = h.si.value
        self.k_B = k_B.si.value
        self.T_CMB = 2.7255  # K

        # Frequency setup
        self.frequencies = np.geomspace(1e10, 2.8e12, 1000)
        self.xf = (self.h * self.frequencies) / (self.k_B * self.T_CMB)
        self.i0 = 2 * (self.k_B * self.T_CMB) ** 3 / (self.h * self.c) ** 2\
                    * 1e20

    def _brokenpowerlaw_int(self, pbreak):
        alpha1 = self.alpha
        alpha2 = self.alpha2
        pmin = self.pmin
        pmax = self.pmax
        xx = self.xf
        spec = np.zeros(len(xx))
        norm = ((pmin**(1-alpha1)-pbreak**(1-alpha1))/(alpha1-1)
                + pbreak ** (-alpha1+alpha2)/(alpha2-1)
                * (pbreak**(1-alpha2)-pmax**(1-alpha2)))**(-1)

        def integrand1(p, alpha, x):
            return photon_scat_kernel(p, x) * p**(-alpha)

        def integrand2(p, alpha, alpha2, x, pbreak):
            return photon_scat_kernel(p, x)\
                * pbreak ** (-alpha1+alpha2) * p ** (-alpha2)
        for i in np.arange(len(xx)):
            spec[i] = quad(integrand1, pmin, pbreak, args=(alpha1, xx[i]))[0]\
                + norm * quad(
                    integrand2, pbreak, pmax,
                    args=(alpha1, alpha2, xx[i], pbreak))[0]
        return norm * spec

    def _save_grid(self, grid_data):
        hdr = fits.Header()
        hdr['COMMENT'] = 'This is a file containing non-thermal SZ spectra due'
        hdr['COMMENT'] = ' to non-thermal population of electrons with momenta'
        hdr['COMMENT'] = 'distribution corresponding to a broken power-law as'
        hdr['COMMENT'] = 'a function of pmin and observing frequency. It is a'
        hdr['COMMENT'] = 'grid where the axes are pmin and frequency in Hz'
        hdr['COMMENT'] = 'and each point corresponds to a value of the corresponding'
        hdr['COMMENT'] = 'distortion in specific intensity due to ntSZ.'
        primary_hdu = fits.PrimaryHDU(header=hdr)

        hdr = fits.Header()
        hdr['NAME'] = 'ntSZ spectra'
        hdr['COMMENT'] = 'alpha1 = 0.05'
        hdr['COMMENT'] = f'alpha2 = {self.alpha2}'
        hdr['COMMENT'] = f'pmin = {self.pmin}'
        hdr['COMMENT'] = 'gridpbr = np.geomspace(100, 1200, 100)'
        hdr['COMMENT'] = 'gridf = np.geomspace(1e10, 2.8e12, 1000)'
        hdr['COMMENT'] = f'pmax = {self.pmax}'
        col = fits.ImageHDU(data=grid_data, header=hdr)
        hdul = fits.HDUList([primary_hdu, col])
        hdul.writeto(self.filename, overwrite=True)

    def create_grid(self):
        pbr = np.geomspace(100, 1200, 100)
        print("Starting integration to create a grid..")

        pool = Pool(processes=self.nproc)
        data = pool.map(self._brokenpowerlaw_int, [pcr for pcr in pbr])
        pool.close()
        doubleparr = np.array(data)
        doubleparr = (doubleparr - i_x(self.xf)) * self.i0

        self._save_grid(self.filename, doubleparr)
        print(f"Grid saved to {self.filename}")

        return doubleparr
