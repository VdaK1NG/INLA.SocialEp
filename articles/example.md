# Usage Example

## Introduction

``` r
library("INLA.SocialEp")
```

## Understanding the inputs

User needs to provide a graph for the function that can be prepared from
an sf object like:

``` r
sec.nb <- spdep::poly2nb(sp_object_sim, snap=0.0000005)
spdep::nb2INLA("MapGraph", sec.nb)
graph_states <- INLA::inla.read.graph("MapGraph")
```

Moreover, the user needs to provide the data with a data.frame that has
the following structure:

- **obs**: contains the observations for the response variable per area
  and group.
- **exp**: contains the expected values per area and group.
- **lev.fac1**: contains the levels of *factor 1* per area and group.
- **lev.fac2**: contains the levels of *factor 2* per area and group.

``` r
data_mod <- data.frame(
  "obs" = c(sp_object_sim$OBS_M6_SIM_G1_v3, sp_object_sim$OBS_M6_SIM_G2_v3, 
            sp_object_sim$OBS_M6_SIM_G3_v3, sp_object_sim$OBS_M6_SIM_G4_v3),
  "exp" = c(sp_object_sim$EXP, sp_object_sim$EXP, sp_object_sim$EXP, sp_object_sim$EXP),
  "lev.fac1" = c(rep("F1L1", graph_states$n*2), rep("F1L2", graph_states$n*2)), 
  "lev.fac2" = c(rep("F2L1", graph_states$n), rep("F2L2", graph_states$n), rep("F2L1", graph_states$n), rep("F2L2", graph_states$n))
)
```

## Customizing parameters

There are three mainly aspects regarding the modelization and
post-processing that can be modified by the user:

- **Prior Specification**: the function allows the user to choose
  between uniform prior (`sp.prior = "pc.prec"`) or Penalized-Complexity
  priors (`sp.prior = "pc.prec"`), with the user being able to modify
  the paremeters of the PC prior using option
  `pc.prec.val = c(1, 0.01)`.
- **Shared Copies in INLA**: you can control if you want the shared
  copies to be fixed with option `sp.copy.fixed=TRUE`.
- **Scaled Spatial Effects**: this option controls whether the random
  spatial effects have their variance scaled (`TRUE`) or not (`FALSE`).

We can now run the model:

``` r
ResMod <- INLA.SocialEp::inla.SpANOVA.2x2(
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
```

And we can check the summary of the results obtained:

``` r
ResMod$Summary
#>    NUMBER                  MODEL      DIC     WAIC    CPU LOOCV LOGCV.3 LOGCV.5
#> 1       1                     M0 92182.37 89865.26   4.66  4.56    4.56    4.56
#> 2       2                     M1 87736.14 85749.56  62.39  3.93    3.91    3.91
#> 3       3      M2-ind(F1L1-F2L1) 91115.03 88700.52  27.85  4.38    4.39    4.39
#> 4       4      M2-ind(F1L2-F2L1) 91116.86 88701.68  21.42  4.38    4.39    4.39
#> 5       5      M2-ind(F2L1-F1L2) 91117.13 88701.91  26.87  4.38    4.39    4.39
#> 6       6      M2-ind(F2L2-F1L2) 91117.41 88702.71  20.85  4.38    4.39    4.39
#> 7       7           M3-F1.(F1L1) 90446.04 88097.93  68.15  4.27    4.27    4.27
#> 8       8           M3-F1.(F1L2) 89684.16 87308.31  68.24  4.16    4.17    4.17
#> 9       9           M4-F2.(F2L1) 90521.51 88118.30  80.47  4.28    4.28    4.28
#> 10     10           M4-F2.(F2L2) 89761.28 87407.14  68.48  4.16    4.18    4.17
#> 11     11 M5-F1.(F1L1)+F2.(F2L1) 88393.32 86095.64 197.44  4.03    4.02    4.01
#> 12     12 M5-F1.(F1L2)+F2.(F2L1) 88102.19 85848.25 170.69  3.99    3.99    3.99
#> 13     13 M5-F1.(F1L1)+F2.(F2L2) 88127.66 85861.21 207.53  3.99    4.00    3.99
#> 14     14 M5-F1.(F1L2)+F2.(F2L2) 88975.35 86618.53 247.20  4.07    4.07    4.06
#> 15     15 M6-F1.(F1L1)*F2.(F2L1) 86042.32 83972.56 406.53  3.75    3.76    3.73
#> 16     16 M6-F1.(F1L2)*F2.(F2L1) 85780.89 83780.32 349.86  3.72    3.77    3.75
#> 17     17 M6-F1.(F1L1)*F2.(F2L2) 85886.89 83945.06 608.46  3.73    3.78    3.75
#> 18     18 M6-F1.(F1L2)*F2.(F2L2) 86825.70 84817.70 402.74  3.84    3.80    3.79
#> 19     19 M6-F2.(F2L1)*F1.(F1L1) 86043.43 83974.37 512.13  3.75    3.76    3.74
#> 20     20 M6-F2.(F2L2)*F1.(F1L1) 85789.72 83790.85 485.41  3.72    3.77    3.75
#> 21     21 M6-F2.(F2L1)*F1.(F1L2) 85823.41 83844.93 435.44  3.73    3.77    3.74
#> 22     22 M6-F2.(F2L2)*F1.(F1L2) 86761.66 84714.51 493.84  3.83    3.81    3.80
#>    LOGCV.10 sp.null
#> 1      4.56       0
#> 2      3.91       0
#> 3      4.39       0
#> 4      4.39       0
#> 5      4.39       0
#> 6      4.39       0
#> 7      4.27       0
#> 8      4.17       0
#> 9      4.28       0
#> 10     4.17       0
#> 11     4.01       0
#> 12     3.99       0
#> 13     3.99       0
#> 14     4.06       0
#> 15     3.72       0
#> 16     3.73       0
#> 17     3.74       0
#> 18     3.79       0
#> 19     3.72       0
#> 20     3.74       0
#> 21     3.73       0
#> 22     3.80       0
```

We can also plot the spatial effects adjusted for each of the models
using the function
[`plotSpANOVA()`](https://vdak1ng.github.io/INLA.SocialEp/reference/plotSpANOVA.md):

``` r
INLA.SocialEp::plotSpANOVA(
  obj=ResMod,
  obj_type="SpANOVA",
  fill_by="Spatial",
  n_mod=20,
  sp_obj=sp_object_sim,
  fil_scale=c("#133BF2", "#7189F7", "#FFFFFF", "#FF867A", "#FF2F1B"),
  scale_name="Values",
  sp_null=0.125,
  legend.position="right",
  ncol_fig=2
)
```

As well as the adjusted relative risks:

``` r
INLA.SocialEp::plotSpANOVA(
  obj=ResMod,
  obj_type="SpANOVA",
  fill_by="RR",
  n_mod=20,
  sp_obj=sp_object_sim,
  fil_scale=c("#133BF2", "#7189F7", "#FFFFFF", "#FF867A", "#FF2F1B"),
  scale_name="Values",
  sp_null=0.125,
  legend.position="right",
  ncol_fig=2
)
```

------------------------------------------------------------------------

## On this page

Developed by [Pablo Escobar-Hernández](https://github.com/VdaK1NG),
[Francisco Palmí-Perales](https://github.com/FranciscoPalmiPerales),
[Antonio López-Quílez](https://github.com/antoloqui).

Site built with [pkgdown](https://pkgdown.r-lib.org/) 2.2.0.
