import numpy as np
from astropy.constants import m_e, k_B, c, h, e
from scipy.integrate import quad, simpson
from scipy.special import kv
from astropy.io import fits
from scipy import interpolate

# Define constants
c = c.si.value
m_e = m_e.si.value
h = h.si.value
k_B = k_B.si.value
e = e.si.value
T_CMB = 2.7255  # K


def freq_to_x(f):
    """Convert frequency f (Hz) to dimensionless x = h * f / (k_B * T_{CMB}).
    Parameter
    ---------
    f: float
       Frequency in Hz
    Returns
    ---------
    x: float
        Dimensionless frequency
    """
    return (h * f) / (k_B * T_CMB)


def x_to_freq(x):
    """Convert dimensionless frequency to frequency in GHz
    f = x * (k_B * T_{CMB}) / h.
    Parameter
    ---------
    x: float
       Dimensionless frequency
    Returns
    ---------
    x: float
        Frequency in GHz
    """
    return (k_B * T_CMB) * x / h / 1e9


def intensity(x, y_0):
    I_0 = 2 * (k_B*T_CMB)**3 / (h*c)**2 * 1e20
    return (
        I_0
        * y_0
        * ((x**4 * np.exp(x)) / (np.exp(x) - 1) ** 2)
        * ((x * (np.exp(x) + 1) / (np.exp(x) - 1)) - 4)
    )


def i_x(x):
    """Compute the blackbody intensity spectrum I(x) ~ x^3 / (e^x - 1).
    Parameters
    ----------
    x: float
       Dimensionless frequency

    Returns
    -------
    I: float
       Specific intensity of a blackbody with T used to estimate the
       dimensionless frequency
    """
    return x**3 / (np.exp(x) - 1)


def photon_scat_kernel(p, x):
    """Compute the photon scattering kernel for a given normalised momentum p
    and dimensionless frequency x.

    Parameters
    ----------
    p: float
        normalised physical momentum of an electron where p = beta_e * gamma_e
    x: float
        Dimensionless frequency

    Returns
    -------
    integrand: float
        The value of the integrated photon scattering kernel
    """
    smax = 2 * np.arcsinh(p)
    s = np.linspace(-smax, smax, 500)
    t = np.exp(s)
    f = []
    for i in np.arange(len(t)):
        func = (
            (
                -3
                * np.absolute(1 - t[i])
                / (32 * p**6 * t[i])
                * (1 + (10 + 8 * p**2 + 4 * p**4) * t[i] + t[i] ** 2)
            )
            + (
                3
                * (1 + t[i])
                / (8 * p**5)
                * (
                    (3 + 3 * p**2 + p**4) / np.sqrt(1 + p**2)
                    - (3 + 2 * p**2)
                    / (2 * p)
                    * (2 * np.arcsinh(p) - np.absolute(np.log(t[i])))
                )
            )
        ) * t[i]
        func = func * ((x / t[i]) ** 3 / (np.exp(x / t[i]) - 1))
        f.append(func)
    return simpson(np.array(f), s)


def thermal_int(x, t, pmin=0.05, pmax=0.9):
    """Integrate over thermal electron momentum distribution.
    This computes the thermal SZ with relativistic correstions by assuming a
    Maxwell-Juettner distribution to describe a blackbody with temperature T.

    Parameters
    ----------
    x: float
       Dimensionless frequency estimated with T_{CMB}
    t: float
       Kinetic temperature in keV of the scattering thermal electrons
    pmin: float, optional
       Minimum normalised momentum of the electron distribution. Default: 0.05
    pmax: float, optional
       Maximum normalised momentum of the electron distribution. Default: 0.9

    Returns
    --------
    jx: float
       Computes the intensity of photons scattered away from dimensionless
       frequency x. The distortion in specific intensity of CMB due to this
       relativistic SZ is then computed by subtracting i_x(x) and multiplying
       by i0.
    """
    psi = (
        1.38064852e-23 * 1.159 * t * 1e7
        ) / (
            9.10938356e-31 * 299792458.0**2
            )

    def f(p, psi, x):
        return (
            photon_scat_kernel(p, x)
            * p**2
            * np.exp(-np.sqrt(1 + p**2) / psi)
            / psi
            / kv(2, (1 / psi))
        )

    return quad(f, pmin, pmax, args=(psi, x))[0]


def f_e(alpha, pmin, pmax):
    """Compute electron power-law distribution.
    Parameters
    ----------
    alpha: float
           Index of the power-law describing the momentum
           distribution of non-thermal electrons
    pmin: float
          Minimum momentum of the power-law describing the momentum
          distribution of non-thermal electrons
    pmax: float
          Maximum momentum of the power-law describing the momentum
          distribution of non-thermal electrons

    Returns
    ---------
    edist: float array
          Non-thermal electron momentum distribution
    p: float array
          Array with corresponding momenta
    """
    p = np.linspace(pmin, pmax, 100)
    edist = (alpha - 1) * p ** (-alpha) / (
        pmin ** (1 - alpha) - pmax ** (1 - alpha)
        )
    return edist, p


def powerlaw_int(alpha, pmin, pmax, x):
    """Integrate power-law electron distribution * photon scattering kernel to
    estimate ntSZ effect due to a non-thermal electron distribution described
    by a single power-law.

    Parameters:
    -----------
    alpha: float
           Index of the power-law describing the momentum
           distribution of non-thermal electrons.
    pmin: float
          Minimum momentum of the power-law describing the momentum
          distribution of non-thermal electrons.
    pmax: float
          Maximum momentum of the power-law describing the momentum
          distribution of non-thermal electrons.
    x: float
       Dimensionless frequency estimated with T_{CMB}. User can use the
       freq_to_x(f) function to compute this for frequency f in Hz.

    Returns
    -------
    jx: float
        The intensity of the photons scattered away from dimensionless
        frequency x due to non-thermal electrons with power-law distribution.
        To get distortion in specific intensity of the CMB, compute
        (jx - i_x(x)) *i0 * ynth.
    """
    norm = (alpha - 1) / (pmin ** (1 - alpha) - pmax ** (1 - alpha))

    def f(p, alpha, x):
        return photon_scat_kernel(p, x) * p ** (-alpha)

    return norm * quad(f, pmin, pmax, args=(alpha, x))[0]


def broken_fe_dist(p, alpha1, alpha2, pmin, pcr, pmax):
    if p < pcr:
        return p ** (-alpha1)
    else:
        return pcr ** (-alpha1 + alpha2) * p ** (-alpha2)


def broken_fe(alpha1, alpha2, pmin, pcr, pmax):
    """Compute broken power-law distribution.
    Parameters:
    -----------
    alpha1: float
           Index of the flat part of the broken power-law describing
           momentum distribution of non-thermal electrons.
    alpha2: float
           Index of the power-law part of the broken power-law describing
           momentum distribution of non-thermal electrons. alpha2 is negative.
    pmin: float
          Minimum momentum of the broken power-law describing the momentum
          distribution of non-thermal electrons.
    pcr: float
          Break momentum of the broken power-law describing the momentum
          distribution of non-thermal electrons.
    pmax: float
          Maximum momentum of the broken power-law describing the momentum
          distribution of non-thermal electrons. Note: pmin < pbreak < pmax.
    x: float
       Dimensionless frequency estimated with T_{CMB}. User can use the
       freq_to_x(f) function to compute this for frequency f in Hz.
    Returns
    ---------
    edist: float array
          Non-thermal electron momentum distribution
    p: float array
          Array with corresponding momenta
    """
    p = np.linspace(pmin, pmax, 500)
    norm = (
        (pmin ** (1 - alpha1) - pcr ** (1 - alpha1)) / (alpha1 - 1)
        + pcr ** (-alpha1 + alpha2)
        / (alpha2 - 1)
        * (pcr ** (1 - alpha2) - pmax ** (1 - alpha2))
    ) ** (-1)
    dist = np.zeros(len(p))
    for i in np.arange(len(p)):
        dist[i] = broken_fe_dist(p[i], alpha1, alpha2, pmin, pcr, pmax)
    return norm * dist, p


def brokenpowerlaw_int(alpha1, alpha2, pmin, pcr, pmax, x):
    """Integrate broken power-law electron distribution * photon scattering
    kernel to estimate ntSZ effect due to a non-thermal electron distribution
    described by a broken power-law.

    Parameters:
    -----------
    alpha1: float
           Index of the flat part of the broken power-law describing
           momentum distribution of non-thermal electrons.
    alpha2: float
           Index of the power-law part of the broken power-law describing
           momentum distribution of non-thermal electrons. alpha2 is negative.
    pmin: float
          Minimum momentum of the broken power-law describing the momentum
          distribution of non-thermal electrons.
    pcr: float
          Break momentum of the broken power-law describing the momentum
          distribution of non-thermal electrons.
    pmax: float
          Maximum momentum of the broken power-law describing the momentum
          distribution of non-thermal electrons. Note: pmin < pbreak < pmax.
    x: float
       Dimensionless frequency estimated with T_{CMB}. User can use the
       freq_to_x(f) function to compute this for frequency f in Hz.

    Returns
    -------
    jx: float
        The intensity of the photons scattered away from dimensionless
        frequency x due to non-thermal electrons with broken power-law
        distribution. To get distortion in specific intensity of the CMB,
        compute (jx - i_x(x)) * i0 * ynth.
    """
    norm = (
        (pmin ** (1 - alpha1) - pcr ** (1 - alpha1)) / (alpha1 - 1)
        + pcr ** (-alpha1 + alpha2)
        / (alpha2 - 1)
        * (pcr ** (1 - alpha2) - pmax ** (1 - alpha2))
    ) ** (-1)

    def f1(p, alpha_1, x):
        return photon_scat_kernel(p, x) * p ** (-alpha_1)

    def f2(p, alpha_1, alpha_2, x, p_cr):
        return photon_scat_kernel(p, x) * p_cr ** (-alpha_1 + alpha_2) * p\
            ** (-alpha_2)

    return (
        norm * quad(f1, pmin, pcr, args=(alpha1, x))[0]
        + norm * quad(f2, pcr, pmax, args=(alpha1, alpha2, x, pcr))[0]
    )


def incompletebeta(x, a, b):
    """Incomplete beta function."""
    def func(t, a, b):
        return t ** (a - 1) * (1 - t) ** (b - 1)

    return quad(func, 0, x, args=(a, b))[0]


def pe(alpha, p1, p2):
    """Compute pressure for a power-law distribution of non-thermal
    electrons."""
    factor = m_e * c**2 * (alpha - 1) / (
        6 * (p1 ** (1 - alpha) - p2 ** (1 - alpha))
        )
    x = [(1 + p1**2) ** (-1), (1 + p2**2) ** (-1)]
    aa = (alpha - 2) / 2
    bb = (3 - alpha) / 2
    integrand = incompletebeta(x[0], aa, bb) - incompletebeta(x[1], aa, bb)
    return factor * integrand


def brokenpe(alpha1, alpha2, pmin, pcr, pmax):
    """Compute pressure for a broken power-law distribution of non-thermal
    electrons."""
    norm = (
        (pmin ** (1 - alpha1) - pcr ** (1 - alpha1)) / (alpha1 - 1)
        + (
            pcr ** (-alpha1 + alpha2)
            / (alpha2 - 1)
            * (pcr ** (1 - alpha2) - pmax ** (1 - alpha2))
        )
    ) ** (-1)
    factor = norm * m_e * c**2 / 6
    x = [(1 + pmin**2) ** (-1), (1 + pcr**2) ** (-1), (1 + pmax**2) ** (-1)]
    aa = [(alpha1 - 2) / 2, (alpha2 - 2) / 2]
    bb = [(3 - alpha1) / 2, (3 - alpha2) / 2]
    integrand = (
        incompletebeta(x[0], aa[0], bb[0])
        - incompletebeta(x[1], aa[0], bb[0])
        + (
            (incompletebeta(x[1], aa[1], bb[1])
                - incompletebeta(x[2], aa[1], bb[1]))
            * pcr ** (-alpha1 + alpha2)
        )
    )
    return factor * integrand


def ksz(ytsz, Te, x):
    """Compute kinetic SZ contribution. The peculiar velocity of a cluster
    is drawn from a Normal distribution with dispersion 100 km/s.
    Parameters
    ----------
    ytsz: float
        Compton-y parameter of thermal SZ.
    Te: float
        Kinetic temperature of the scattering electrons in keV.
    x: float
       Dimensionless frequency estimated with T_{CMB}. User can use the
       freq_to_x(f) function to compute this for frequency f in Hz.

    Returns
    -------
    jx: float
        The contribution of kSZ effect. To estimate the distortion in
        specific intensity of the CMB due to kSZ, compute
        (jx - i_x(x)) * i0 * ytSZ.
    """
    vpec = np.random.normal(scale=100)
    h_x = (x**4)*np.exp(x)/(np.exp(x)-1)**2
    tau = ytsz * m_e * c**2 / (Te * 1e3 * e)
    ksz = (-1) * tau * vpec * 1000 / c * h_x
    return ksz


# The following functions estimate ntSZ spectra through interpolation
def ntsz_single_interpol(pmin, freq, path_to_grid, alpha=3.61):
    """Interpolate non-thermal SZ spectra for single power-law momentum
    distribution of non-thermal electrons. This function can be used to compute
    multiple spectra while keeping all variables of the single power-law fixed
    except for pmin. This function interpolates the the spectra from a
    pre-computed grid specified by path_to_grid. If user doesn't provide a
    grid, this function assumes a pre-computed grid with an alpha=3.61 and
    pmax=5e5 describing the power-law.

    Parameters
    ----------
    pmin: float
          Minimum momentum of the power-law describing the momentum
          distribution of non-thermal electrons.
    freq: float
          Observing frequency in GHz.
    path_to_grid: str, optional
          Path to ntSZ grid FITS file computed by user before using this
          function. The class NTSZGridGenerator_single() can be used to
          compute such a grid for single power-law describing the
          non-thermal electrons.
    alpha: float, optional
          Power of the power-law describing momenta distribution of
          non-thermal electrons if user providing their own ntSZ grid.
          Default: 3.61.

    Returns:
    --------
    h_x: float
         Intensity contribution due to ntSZ effect. The distortion in the
         specific intensity of the CMB is computed by
         h_x * ynth.
    """
    hdu = fits.open(path_to_grid)
    ntsz_finegrid = np.transpose(hdu[1].data)
    hdu.close()
    pm1 = np.linspace(-1, +3.0, 200)
    pm = 10**pm1
    gridf = np.geomspace(1e10, 2.8e12, 1000)
    ntsz_interpol = interpolate.RectBivariateSpline(
        gridf, pm, ntsz_finegrid, kx=1, ky=1
    )
    pp1 = pe(alpha, pmin, 5e5)
    h_x = ntsz_interpol(freq * 1e9, pmin)[:, 0] * m_e * c**2 / pp1
    return h_x


def ntsz_broken_interpol(
    pbreak, freq, path_to_grid, pmin=1.0, alpha2=3.61
        ):
    """Interpolate non-thermal SZ spectra for broken power-law momentum
    distribution of non-thermal electrons. This function can be used to compute
    multiple spectra while keeping all variables of the broken power-law fixed
    except for pbreak. This function interpolates the the spectra from a
    pre-computed grid specified by path_to_grid. If user doesn't provide a
    grid, this function assumes a pre-computed grid with an alpha2=3.61 and
    pmax=5e5 describing the power-law.

    Parameters
    ----------
    pbreak: float
          Break momentum of the broken power-law describing the momentum
          distribution of non-thermal electrons.
    freq: float
          Observing frequency in GHz.
    path_to_grid: str, optional
          Path to ntSZ grid FITS file computed by user before using this
          function. The class NTSZGridGenerator_broken() can be used to
          compute such a grid for single power-law describing the
          non-thermal electrons.
    pmin: float, optional
          Minimum momentum of the broken power-law describing the momentum
          distribution of non-thermal electrons. Default: 1.0
    alpha2: float, optional
          Power of the broken power-law describing momenta distribution of
          non-thermal electrons if user providing their own ntSZ grid.
          Default: 3.61.

    Returns:
    --------
    h_x: float
         Intensity contribution due to ntSZ effect. The distortion in the
         specific intensity of the CMB is then computed by
         h_x * ynth.
    """
    hdu = fits.open(path_to_grid)
    ntsz_finegrid = np.transpose(hdu[1].data)
    hdu.close()
    gridf = np.geomspace(1e10, 2.8e12, 1000)
    gridpbr = np.geomspace(100, 1200, 100)
    ntsz_interpol = interpolate.RectBivariateSpline(
        gridf, gridpbr, ntsz_finegrid, kx=1, ky=1
    )
    pp1 = brokenpe(0.05, alpha2, pmin, pbreak, 5e5)
    h_x = ntsz_interpol(freq * 1e9, pbreak)[:, 0] * m_e * c**2 / pp1
    return h_x
