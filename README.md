
<!-- README.md is generated from README.Rmd. Please edit that file -->

# INLA.SocialEp <a href='https://github.com/VdaK1NG/INLA.SocialEp'><img src='man/figures/logo.png' align="right" height="120" /></a>

Package that includes a series of functions acting as a wrapper for
several disease mapping methods applied to social epidemiology on INLA.

<!-- badges: start -->

[![license](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/SchlossLab/mikropml/blob/main/LICENSE.md)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03073/status.svg)](https://doi.org/10.21105/joss.03073)
<!-- badges: end -->

An interface to build machine learning models for classification and
regression problems. `INLA.SocialEp` implements the ML pipeline
described by [Escobar-Hernández *et al.*
(2025)](https://doi.org/10.48550/arXiv.2410.21227) with reasonable
default options for data preprocessing, hyperparameter tuning,
cross-validation, testing, model evaluation, and interpretation steps.
See the [website](http://www.schlosslab.org/mikropml/) for more
information, documentation, and examples.

## Installation

You can install the latest development version from
[GitHub](https://github.com/VdaK1NG/INLA.SocialEp):

``` r
install.packages("remotes")
remotes::install_github("VdaK1NG/INLA.SocialEp")
```

### Dependencies

**-Imports**: dplyr,ggraph, spdep, ggplot2, gridExtra, broom, sf, tidyr,
scales, SUMMER

------------------------------------------------------------------------

**-Suggests**: assertthat, doFuture, forcats

------------------------------------------------------------------------

## Usage

Check out the [introductory
vignette](http://www.schlosslab.org/mikropml/articles/introduction.html)
for a quick start tutorial. For a more in-depth discussion, read [all
the vignettes](http://www.schlosslab.org/mikropml/articles/index.html)
and/or take a look at the [reference
documentation](http://www.schlosslab.org/mikropml/reference/index.html).

## Help & Contributing

If you come across a bug, [open an
issue](https://github.com/SchlossLab/mikropml/issues) and include a
[minimal reproducible example](https://tidyverse.org/help/).

## License

The `INLA.SocialEp` package is licensed under [the MIT
license](https://github.com/VdaK1NG/INLA.SocialEp/LICENSE.md). Text and
images included in this repository, including the mikropml logo, are
licensed under the [CC BY 4.0
license](https://creativecommons.org/licenses/by/4.0/).

## Citation

To cite `INLA.SocialEp` in publications, use:

A BibTeX entry for LaTeX users is:

     @Manual{,
      title = {INLA.SocialEp: Implementation of Advanced Social Disease Mapping Techniques in
    INLA},
      author = {Pablo Escobar-Hernández and Francisco Palmí-Perales and Antonio López-Quílez},
      year = {2025},
      note = {R package version 1.1, commit 0edbfae61f99d7d65045fd9be9dccfd2327445bc},
      url = {https://github.com/VdaK1NG/INLA.SocialEp},
    } 
