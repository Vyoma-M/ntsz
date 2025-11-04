# ntsz 
This package enables the compution of the spectral distortion in the specific intensity of the CMB due to the Sunyaev-Zeldovich effect for the following models describing the momenta distribution of the scattering electrons:
1) Maxwell-Juettner distribution with a given kinetic temperature. 
2) Single power-law with a minimum and maximum momenta, and an index.
3) Broken power-law with a flat spectrum from minimum to break momentum and then a power-law with a negative index to maximum momentum.

While the first model results in the thermal SZ spectrum with relativistic corrections (referred to as rSZ), the power-law models result in the non-thermal SZ (ntSZ) effect.

## Installation
With ```pip```:

` pip install ntsz `

The following modules need to be installed to use this package:

```
healpy
numpy
astropy
multiprocessing
scipy
```

## Usage
The [example_notebook.ipynb](example_notebook.ipynb) can be used to familiarise oneself with the package.
The user is referred to documentation for further information on other tools available with this package.

## References
The computation of the SZ spectra follows the formalism described here:

[1] Ensslin & Kaiser (2000): Comptonization of the Cosmic Microwave Background by Relativistic Plasma, [arXiv:0001429v2](https://arxiv.org/abs/astro-ph/0001429v2).

The kSZ effect is estimated using the formalism presented in:

[2] Mroczkowski, T., et. al. (2019): Astrophysics with the Spatially and Spectrally Resolved Sunyaev-Zeldovich Effects: A Millimetre/Submillimetre Probe of the Warm and Hot Universe, [arXiv:1811.02310](https://arxiv.org/abs/1811.02310).

## Citation
If you use this package for your work, please cite the following papers:

[1] Vyoma Muralidhara & Kaustuv Basu (2024): Constraining the average magnetic field in galaxy clusters with current and upcoming CMB surveys, [DOI: 10.1088/1475-7516/2024/11/010](https://doi.org/10.1088/1475-7516/2024/11/010).

[2] Ensslin & Kaiser (2000): Comptonization of the Cosmic Microwave Background by Relativistic Plasma, [arXiv:0001429v2](https://arxiv.org/abs/astro-ph/0001429v2).
