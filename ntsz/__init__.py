# -*- coding: utf-8 -*-
#
#  This file is part of ntsz.
#
#  ntsz is free software; you can redistribute it and/or modify
#  it under the terms of the MIT License.
#
#  ntsz is distributed in the hope that it will be useful,but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See
#  the provided copy of the MIT License for more details.

"""
ntsz Package
A package for calculating non-thermal Sunyaev-Zel'dovich (ntSZ) effect,
the associated distortion in the specific intensity of the CMB and
related Planck bandpass corrections.

Modules:
--------
- ntsz: Core functions to compute the following:
        1) Photon scattering kernel
        2) Pseudo-temperature (non-thermal equivalent of thermal kinetic temp.)
        3) ntSZ spectra for single/broken power-law distributions of
           the momenta of non-thermal electrons.
        4) Thermal SZ effect with relativistic corrections for a given temp.
- ntsz_utils: Utility and file-handling functions.
- planck_bandpass: Estimation of Planck bandpass-corrected spectra.
- ntsz_grid_generator: Compute grids of ntSZ spectra in order to save
        computation time when computing multiple spectra for variable pmin for
        single power-law or variable pbreak for broken power-law models. User
        is recommended to generate these grids before computing
        bandpass-corrected spectra.

Data products:
-------------
The following data products are included with this package:
- ntsz_grid_single.fits: This is a grid of the estimated ntSZ spectra for
                    different pmin values of the power-law.
- ntsz_grid_double.fits: This is a grid of the estimated ntSZ spectra for
                    different pbreak values of the broken power-law.
- planck-bandpass_fits: A FITS file containing the wavenumbers/frequencies at
                        which the transmissions are measured. This file is
                        created by extracting information from the
                        HFI_RIMO_3.31.fits and LFI_RIMO_3.30.fits files that
                        are publicly available in the Planck Legacy Archive.
"""

from .ntsz import (
    freq_to_x,
    x_to_freq,
    intensity,
    i_x,
    photon_scat_kernel,
    thermal_int,
    f_e,
    powerlaw_int,
    broken_fe,
    brokenpowerlaw_int,
    ksz,
    ntsz_single_interpol,
    ntsz_broken_interpol
    )

from .planck_bandpass import (
    BandpassCorrector
)

from .ntsz_utils import (
    get_data_path,
    get_path,
    writefile,
    create_fits,
    read_fits
)
from .ntsz_grid_generator import (
    NTSZGridGenerator_broken,
    NTSZGridGenerator_single
)

__all__ = [
    # ntsz functions
    "freq_to_x",
    "x_to_freq",
    "intensity",
    "i_x",
    "photon_scat_kernel",
    "thermal_int",
    "f_e",
    "powerlaw_int",
    "broken_fe",
    "brokenpowerlaw_int",
    "ksz",
    "ntsz_single_interpol",
    "ntsz_broken_interpol",
    # bandpass_corrections functions
    "BandpassCorrector",
    # ntsz utility functions
    "get_data_path",
    "get_path",
    "writefile",
    "create_fits",
    "read_fits",
    # ntSZ spectra grid generator classes
    "NTSZGridGenerator_broken",
    "NTSZGridGenerator_single"
]
