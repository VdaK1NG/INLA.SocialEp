# INLA.SocialEp

  

Package that includes a series of functions acting as a wrapper for
several disease mapping methods applied to social epidemiology on INLA,
as described by [Escobar-Hernández *et al.*
(2025)](https://doi.org/10.48550/arXiv.2410.21227). See the
[website](https://vdak1ng.github.io/INLA.SocialEp) for more information,
documentation, and examples.

## Installation

You can install the latest development version from
[GitHub](https://github.com/VdaK1NG/INLA.SocialEp):

``` r
install.packages("remotes")
remotes::install_github("VdaK1NG/INLA.SocialEp")
```

### Dependencies

**-Imports**: dplyr, ggraph, spdep, ggplot2, gridExtra, broom, sf,
tidyr, scales, SUMMER

------------------------------------------------------------------------

**-Suggests**: assertthat, doFuture, forcats

------------------------------------------------------------------------

## Usage

Check out the [Example
vignette](https://vdak1ng.github.io/INLA.SocialEp/articles/example.html)
for a quick start tutorial. If you are interested in how the datasets
were simulated, feel free to check the in-depth description at the
[Simulation Study - Data Creation
vignette](https://vdak1ng.github.io/INLA.SocialEp/articles/simulation.html).
For a more in-depth discussion, read and/or take a look at the
[reference
documentation](https://vdak1ng.github.io/INLA.SocialEp/reference/index.html).

## Help & Contributing

If you come across a bug, [open an
issue](https://github.com/VdaK1NG/INLA.SocialEp/issues) and include a
[minimal reproducible example](https://tidyverse.org/help/).

## License

The `INLA.SocialEp` package is licensed under [the MIT
license](https://github.com/VdaK1NG/INLA.SocialEp/LICENSE.md). Text and
images included in this repository, including the INLA.SocialEp logo,
are licensed under the [CC BY 4.0
license](https://creativecommons.org/licenses/by/4.0/).

## Citation

To cite `INLA.SocialEp` in publications, use:

A BibTeX entry for LaTeX users is:

``` R
 @Manual{,
  title = {INLA.SocialEp: Implementation of Advanced Social Disease Mapping Techniques in
INLA},
  author = {Pablo Escobar-Hernández and Francisco Palmí-Perales and Antonio López-Quílez},
  year = {2025},
  note = {R package version 1.1, commit 0edbfae61f99d7d65045fd9be9dccfd2327445bc},
  url = {https://github.com/VdaK1NG/INLA.SocialEp},
} 
```
