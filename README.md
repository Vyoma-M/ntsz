# ntsz 
This package enables the compution of the spectral distortion in the specific intensity of the CMB due to the Sunyaev-Zeldovich effect for the following models describing the momenta distribution of the scattering electrons:
1) Maxwell-Juettner distribution with a given kinetic temperature. 
2) Single power-law with a minimum and maximum momenta, and an index.
3) Broken power-law with a flat spectrum from minimum to break momentum and then a power-law with a negative index to maximum momentum.

While the first model results in the thermal SZ spectrum with relativistic corrections (referred to as rSZ), the power-law models result in the non-thermal SZ (ntSZ) effect.

![Static Badge](https://img.shields.io/badge/GitHub-Vyoma--M-blue?link=https%3A%2F%2Fgithub.com%2FVyoma-M)
![Static Badge](https://img.shields.io/badge/arXiv-2402.17445-green?link=https%3A%2F%2Farxiv.org%2Fabs%2F2402.17445)
![Static Badge](https://img.shields.io/badge/License-MIT-red?link=https%3A%2F%2Fgithub.com%2FVyoma-M%2Fntsz%2Fblob%2Fmain%2FLICENSE)


## Installation
With ```pip```:

``` pip install ntsz ```

## Test
Test your installation with:

```python -m unittest discover```

## Usage
The [example_notebook.ipynb](example_notebook.ipynb) can be used to familiarise oneself with the package.
The user is referred to documentation for further information on other tools available with this package.

## References
The computation of the SZ spectra follows the formalism described here:

[1] Ensslin & Kaiser (2000): Comptonization of the Cosmic Microwave Background by Relativistic Plasma, [arXiv:0001429v2](https://arxiv.org/abs/astro-ph/0001429v2).

The kSZ effect is estimated using the formalism presented in:

[2] Mroczkowski, T., et. al. (2019): Astrophysics with the Spatially and Spectrally Resolved Sunyaev-Zeldovich Effects: A Millimetre/Submillimetre Probe of the Warm and Hot Universe, [arXiv:1811.02310](https://arxiv.org/abs/1811.02310).

## Acknowledgement
If you found this package useful for your work, please cite the following papers:

[1] Vyoma Muralidhara & Kaustuv Basu (2024): Constraining the average magnetic field in galaxy clusters with current and upcoming CMB surveys, [DOI: 10.1088/1475-7516/2024/11/010](https://doi.org/10.1088/1475-7516/2024/11/010). BibTex entry:
```
@article{Muralidhara:2024ipg,
    author = "Muralidhara, Vyoma and Basu, Kaustuv",
    title = "{Constraining the average magnetic field in galaxy clusters with current and upcoming CMB surveys}",
    eprint = "2402.17445",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.CO",
    doi = "10.1088/1475-7516/2024/11/010",
    journal = "JCAP",
    volume = "11",
    pages = "010",
    year = "2024"
}
```

[2] Ensslin & Kaiser (2000): Comptonization of the Cosmic Microwave Background by Relativistic Plasma, [arXiv:0001429v2](https://arxiv.org/abs/astro-ph/0001429v2). BibTex entry:
```
@article{Ensslin:2000mk,
    author = "Ensslin, Torsten A. and Kaiser, Christian R.",
    title = "{Comptonization of the cosmic microwave background by relativistic plasma}",
    eprint = "astro-ph/0001429",
    archivePrefix = "arXiv",
    journal = "Astron. Astrophys.",
    volume = "360",
    pages = "417--430",
    year = "2000"
}
```

## License

Copyright 2025 Vyoma Muralidhara.

ntsz is free software made available under the MIT License. For details see the LICENSE file.
