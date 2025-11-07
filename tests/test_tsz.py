"""
A unit test to compute rSZ using the Itoh et al. 1998 method and compare it
with the thermal SZ spectrum computed using the ntsz package's thermal_int.
"""
import unittest
from ntsz.ntsz import thermal_int, freq_to_x, i_x
from astropy.constants import k_B, c, h
import numpy as np

# Define constants
c = c.si.value
h = h.si.value
k_B = k_B.si.value
T_CMB = 2.7255  # K


class TestThermalSZ_rel(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        self.temp = 5.0  # Example temperature in keV
        self.i0 = 2 * (k_B * T_CMB) ** 3 / (h * c) ** 2 * 1e20  # in Jy sr^-1
        self.f = np.array([217 * 1e9])
        self.xf = freq_to_x(self.f)
        super().__init__(*args, **kwargs)

    def _itoh_rsz(self):
        sbar = self.xf / np.sinh(self.xf / 2)
        xbar = self.xf * np.cosh(self.xf / 2) / np.sinh(self.xf / 2)
        no = 1 / (np.exp(self.xf) - 1)
        taurel = (9.10938356e-31 * 299792458.0**2)\
            / (1.38064852e-23 * 1.159 * self.temp * 1e7)
        Y0 = -4 + xbar
        Y1 = (
            -10
            + ((47 / 2) * xbar)
            - ((42 / 5) * xbar**2)
            + ((7 / 10) * xbar**3)
            + (sbar**2 * ((-21 / 5) + ((7 / 5) * xbar)))
        )
        Y2 = (
            (-15 / 2)
            + (1023 / 8 * xbar)
            - (868 / 5 * xbar**2)
            + (329 / 5 * xbar**3)
            - (44 / 5 * xbar**4)
            + (11 / 30 * xbar**5)
            + (
                sbar**2
                * ((-434 / 5) + (658 / 5 * xbar)
                    - (242 / 5 * xbar**2) + (143 / 30 * xbar**3))
            )
            + (sbar**4 * ((-44 / 5) + (187 / 60 * xbar)))
        )
        dn = (
            self.xf
            * np.exp(self.xf)
            * (Y0 + Y1 / taurel + Y2 / taurel ** 2)
            * no
            * self.xf**3
            / (np.exp(self.xf) - 1)
            / taurel
        ) * self.i0
        return dn

    def test_thermal_int(self):
        tsz_spec = np.zeros(len(self.xf))
        # Compute thermal SZ spectrum
        for i in range(len(self.xf)):
            tsz_spec[i] = (thermal_int(pmin=1e-8, pmax=0.9,
                                       x=self.xf[i], t=self.temp)
                           - i_x(self.xf[i])) * self.i0
            rsz = self._itoh_rsz()
        self.assertEqual(np.around(tsz_spec, 2), np.around(rsz, 2))


if __name__ == "__main_":
    unittest.main()
