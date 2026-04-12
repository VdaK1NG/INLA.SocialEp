# Usage Example

## Introduction

``` r
install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
#> Installing package into '/home/runner/work/_temp/Library'
#> (as 'lib' is unspecified)
#> Warning: dependencies 'Rgraphviz', 'graph' are not available
#> also installing the dependencies 'colorspace', 'fracdiff', 'timeDate', 'urca', 'forecast', 'modelr', 'microbenchmark', 'doBy', 'nloptr', 'reformulas', 'RcppEigen', 'carData', 'abind', 'pbkrtest', 'lme4', 'DEoptimR', 'RcppArmadillo', 'xopen', 'corrplot', 'car', 'rbibutils', 'SparseM', 'tensorA', 'robustbase', 'bayesm', 'ellipsis', 'miniUI', 'profvis', 'rcmdcheck', 'rversions', 'urlchecker', 'ggsci', 'cowplot', 'ggsignif', 'polynom', 'rstatix', 'litedown', 'dfidx', 'Formula', 'zoo', 'lmtest', 'statmod', 'Rdpack', 'coda', 'httpuv', 'later', 'otel', 'promises', 'sourcetools', 'xtable', 'mnormt', 'quantreg', 'INLAtools', 'fmesher', 'MatrixModels', 'compositions', 'Deriv', 'Ecdat', 'devtools', 'doParallel', 'evd', 'fastGHQuad', 'ggpubr', 'gsl', 'markdown', 'matrixStats', 'mlogit', 'mvtnorm', 'pixmap', 'rgl', 'runjags', 'shiny', 'sn', 'splancs', 'tidyterra', 'INLAspacetime'
library("INLA.SocialEp")
```

## Understanding the inputs

**DISCLAIMER:** Observed and expected values have to be given in an
specific order. Consider n the number of areas, the first n values (1:n)
of the obs (exp) should be the ones belonging to the FIRST level (the
first position of the **lev.fac1** vector argument) of the FIRST factor
(the first position of the **fac.names** vector argument) and to the
FIRST level (the first position of the **lev.fac2** vector argument) of
the SECOND factor (the second position of the **fac.names** vector
argument). The n following values ((n+1):2n) of the obs (exp) should be
the ones belonging to the FIRST level (the first position of the
**lev.fac1** vector argument) of the FIRST factor (the first position of
the fac.names vector argument) and to the SECOND level (the second
position of the **lev.fac2** vector argument) of the SECOND factor (the
second position of the **fac.names** vector argument).  
  
The n following values ((2n+1):3n) of the obs (exp) should be the ones
belonging to the SECOND level (the second position of the **lev.fac1**
vector argument) of the FIRST factor (the first position of the
**fac.names** vector argument) and to the FIRST level (the first
position of the **lev.fac2** vector argument) of the SECOND factor (the
second position of the **fac.names** vector argument).  
  
The n following values ((3n+1):4n) of the obs (exp) should be the ones
belonging to the SECOND level (the second position of the **lev.fac1**
vector argument) of the FIRST factor (the first position of the
**fac.names** vector argument) and to the SECOND level (the second
position of the **lev.fac2** vector argument) of the SECOND factor (the
second position of the **fac.names** vector argument).  
  
The first n values are O1, the second n values are O2, the third n
values are O3 and the last n values are O4.

``` r
sec.nb <- spdep::poly2nb(sp_object_sim, snap=0.0000005)
spdep::nb2INLA("MapGraph", sec.nb)
graph_states <- INLA::inla.read.graph("MapGraph")
```

## Customizing parameters

There are three mainly aspects regarding the modelization and
post-processing that can be modified by the user:

- **Prior Specification**:
- **Shared Copies in INLA**:
- **Saving Options**:

### Prior specification

### Shared Copies in INLA

### Saving Inputs

## Modelling

### Understanding the output

### High Dimension Datasets

In case that you are working with a high dimensional datasets where the
number of areas is high, we advise running the `inla.SpANOVA.2x2` with
the option `save.mod.data=TRUE`, which will save all the necessary
information to run the desired specification afterwards, and all the
saving inputs options as `FALSE`. This will allow you to run through all
the possible configurations without saving the complete results, which
should speed things up and avoid saving undesired results.

``` r
data_mod <- data.frame(
  "obs" = c(sp_object_sim$OBS_M6_SIM_G1_v3, sp_object_sim$OBS_M6_SIM_G2_v3, 
          sp_object_sim$OBS_M6_SIM_G3_v3, sp_object_sim$OBS_M6_SIM_G4_v3),
  "exp" = c(sp_object_sim$EXP, sp_object_sim$EXP, sp_object_sim$EXP, sp_object_sim$EXP),
  "lev.fac1" = c(rep("F1L1", graph_states$n*2), rep("F1L2", graph_states$n*2)), 
  "lev.fac2" = c(rep("F2L1", graph_states$n), rep("F2L2", graph_states$n), rep("F2L1", graph_states$n), rep("F2L2", graph_states$n))
)


ResMod <- inla.SpANOVA.2x2(
  data = data_mod,
  gr = graph_states,
  fac.names = c("F1", "F2"),
  lev.fac1 = c("F1L1", "F1L2"),
  lev.fac2 = c("F2L1", "F2L2"),
  scale.mod = FALSE,
  sp.prior = "pc.prec",
  pc.prec.val = c(1, 0.01),
  sp.copy.fixed = FALSE,
  save.res=FALSE, 
  save.random=FALSE,
  save.hyper=FALSE, 
  save.mod.data=TRUE
  )
#> -/-/-/-/-/-/-/-/-/- M0 finished -/-/-/-/-/-/-/-/-/- 
#> -/-/-/-/-/-/-/-/-/- M1 finished -/-/-/-/-/-/-/-/-/- 
#> -/-/-/-/-/-/-/-/-/- M2 finished -/-/-/-/-/-/-/-/-/- 
#> -/-/-/-/-/-/-/-/-/- M3 finished -/-/-/-/-/-/-/-/-/- 
#> -/-/-/-/-/-/-/-/-/- M4 finished -/-/-/-/-/-/-/-/-/- 
#>   |                                                          |                                                  |   0%  |                                                          |============                                      |  25%  |                                                          |=========================                         |  50%  |                                                          |======================================            |  75%  |                                                          |==================================================| 100%
#> -/-/-/-/-/-/-/-/-/- M5 finished -/-/-/-/-/-/-/-/-/-
#>   |                                                          |                                                  |   0%  |                                                          |======                                            |  12%  |                                                          |============                                      |  25%  |                                                          |===================                               |  38%  |                                                          |=========================                         |  50%  |                                                          |===============================                   |  62%  |                                                          |======================================            |  75%  |                                                          |============================================      |  88%  |                                                          |==================================================| 100%
#> -/-/-/-/-/-/-/-/-/- M6 finished -/-/-/-/-/-/-/-/-/-
#>    NUMBER                  MODEL      DIC     WAIC        CPU LOOCV LOGCV.3
#> 1       1                     M0 92182.37 89865.26   4.598696  4.56    4.56
#> 2       2                     M1 87736.15 85749.54  65.048164  3.93    3.91
#> 3       3      M2-ind(F1L1-F2L1) 91115.00 88700.55  29.747030  4.38    4.39
#> 4       4      M2-ind(F1L2-F2L1) 91116.86 88701.68  21.517702  4.38    4.39
#> 5       5      M2-ind(F2L1-F1L2) 91117.27 88701.81  27.819517  4.38    4.39
#> 6       6      M2-ind(F2L2-F1L2) 91117.40 88702.72  22.109857  4.38    4.39
#> 7       7           M3-F1.(F1L1) 90446.04 88097.93  69.526960  4.27    4.27
#> 8       8           M3-F1.(F1L2) 89684.16 87308.31  68.850652  4.16    4.17
#> 9       9           M4-F2.(F2L1) 90521.51 88118.29  84.396181  4.28    4.28
#> 10     10           M4-F2.(F2L2) 89761.28 87407.14  71.692321  4.16    4.18
#> 11     11 M5-F1.(F1L1)+F2.(F2L1) 88393.50 86095.59 199.699392  4.03    4.02
#> 12     12 M5-F1.(F1L2)+F2.(F2L1) 88102.19 85848.25 168.984106  3.99    3.99
#> 13     13 M5-F1.(F1L1)+F2.(F2L2) 88127.22 85861.23 213.742012  3.99    4.00
#> 14     14 M5-F1.(F1L2)+F2.(F2L2) 88975.12 86618.34 231.311736  4.07    4.07
#> 15     15 M6-F1.(F1L1)*F2.(F2L1) 86041.82 83972.49 442.670633  3.75    3.76
#> 16     16 M6-F1.(F1L2)*F2.(F2L1) 85781.50 83780.15 373.912864  3.72    3.77
#> 17     17 M6-F1.(F1L1)*F2.(F2L2) 85890.56 83943.71 452.746824  3.73    3.78
#> 18     18 M6-F1.(F1L2)*F2.(F2L2) 86825.43 84816.19 545.709757  3.84    3.80
#> 19     19 M6-F2.(F2L1)*F1.(F1L1) 86043.25 83974.04 407.655375  3.75    3.76
#> 20     20 M6-F2.(F2L2)*F1.(F1L1) 85788.89 83792.78 401.011778  3.72    3.77
#> 21     21 M6-F2.(F2L1)*F1.(F1L2) 85822.91 83844.67 351.595838  3.73    3.77
#> 22     22 M6-F2.(F2L2)*F1.(F1L2) 86761.41 84713.34 536.913172  3.83    3.81
#>    LOGCV.5 LOGCV.10
#> 1     4.56     4.56
#> 2     3.91     3.91
#> 3     4.39     4.39
#> 4     4.39     4.39
#> 5     4.39     4.39
#> 6     4.39     4.39
#> 7     4.27     4.27
#> 8     4.17     4.17
#> 9     4.28     4.28
#> 10    4.17     4.17
#> 11    4.01     4.01
#> 12    3.99     3.99
#> 13    3.99     3.99
#> 14    4.06     4.06
#> 15    3.73     3.72
#> 16    3.74     3.73
#> 17    3.75     3.74
#> 18    3.79     3.79
#> 19    3.74     3.72
#> 20    3.75     3.74
#> 21    3.74     3.73
#> 22    3.80     3.80

ResMod$summary
#> NULL
```

After that we can re-run the model of choice, which in this case is
BLABLA. We are free to change some of the parameters now, so we will try
the sdunif prior, and setting `scale.mod=TRUE` and `sp.copy.fixed=TRUE`:

``` r
ResMod_fin <- inla.rerun.SpANOVA(
  SpANOVA_obj = ResMod, 
  n_mod = 20, 
  gr = graph_states, 
  scale.mod = TRUE, 
  sp.prior = "sdunif", 
  sp.copy.fixed = TRUE
)
```

Since it is an INLA object we can run a summary:

``` r
summary(ResMod)
```

We can also check the spatial effects adjusted:

``` r
plot.SpANOVA(
  obj=ResMod,
  obj_type="SpANOVA",
  fill_by="Spatial",
  n_mod=20,
  sp_object=sp_object_sim,
  breaks=c(-3, -2, -1, -0.5, -0.1, 0.1, 0.5, 1, 2, 3),
  fil_scale=c("#133BF2", "#7189F7", "#FFFFFF", "#FF867A", "#FF2F1B"),
  scale_name="Values",
  sp_null=0.125,
  legend.position="right",
  ncol_fig=2
)
```

As well as the adjusted relative risks:

``` r
plot.SpANOVA(
  obj=ResMod,
  obj_type="SpANOVA",
  fill_by="RR",
  n_mod=20,
  sp_object=sp_object_sim,
  breaks=c(0, 0.5, 0.9, 1.1, 2, 3),
  fil_scale=c("#133BF2", "#7189F7", "#FFFFFF", "#FF867A", "#FF2F1B"),
  scale_name="Values",
  sp_null=0.125,
  legend.position="right",
  ncol_fig=2
)
```

## References

------------------------------------------------------------------------

## On this page

Developed by [Pablo Escobar-Hernández](https://github.com/VdaK1NG),
[Francisco Palmí-Perales](https://github.com/FranciscoPalmiPerales),
[Antonio López-Quílez](https://github.com/antoloqui).

Site built with [pkgdown](https://pkgdown.r-lib.org/) 2.2.0.
