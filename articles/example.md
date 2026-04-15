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
#>    NUMBER                  MODEL      DIC     WAIC        CPU LOOCV LOGCV.3
#> 1       1                     M0 92182.37 89865.26   4.412914  4.56    4.56
#> 2       2                     M1 87736.16 85749.54  68.748523  3.93    3.91
#> 3       3      M2-ind(F1L1-F2L1) 91115.09 88700.49  30.726008  4.38    4.39
#> 4       4      M2-ind(F1L2-F2L1) 91116.86 88701.68  22.844925  4.38    4.39
#> 5       5      M2-ind(F2L1-F1L2) 91117.28 88701.80  29.020304  4.38    4.39
#> 6       6      M2-ind(F2L2-F1L2) 91117.44 88702.66  22.667647  4.38    4.39
#> 7       7           M3-F1.(F1L1) 90446.04 88097.93  71.984813  4.27    4.27
#> 8       8           M3-F1.(F1L2) 89684.16 87308.31  71.501817  4.16    4.17
#> 9       9           M4-F2.(F2L1) 90521.50 88118.28  86.686704  4.28    4.28
#> 10     10           M4-F2.(F2L2) 89761.28 87407.14  72.742198  4.16    4.18
#> 11     11 M5-F1.(F1L1)+F2.(F2L1) 88393.44 86095.40 208.292671  4.03    4.02
#> 12     12 M5-F1.(F1L2)+F2.(F2L1) 88102.19 85848.25 170.653224  3.99    3.99
#> 13     13 M5-F1.(F1L1)+F2.(F2L2) 88127.10 85861.30 211.727139  3.99    4.00
#> 14     14 M5-F1.(F1L2)+F2.(F2L2) 88975.18 86618.79 246.437363  4.07    4.07
#> 15     15 M6-F1.(F1L1)*F2.(F2L1) 86041.99 83972.57 402.481428  3.75    3.76
#> 16     16 M6-F1.(F1L2)*F2.(F2L1) 85779.43 83779.75 339.454437  3.72    3.77
#> 17     17 M6-F1.(F1L1)*F2.(F2L2) 85888.95 83943.00 533.309143  3.73    3.78
#> 18     18 M6-F1.(F1L2)*F2.(F2L2) 86826.41 84816.98 667.051420  3.84    3.80
#> 19     19 M6-F2.(F2L1)*F1.(F1L1) 86043.02 83967.15 293.247859  3.75    3.76
#> 20     20 M6-F2.(F2L2)*F1.(F1L1) 85788.90 83792.87 502.294037  3.72    3.77
#> 21     21 M6-F2.(F2L1)*F1.(F1L2) 85822.39 83845.67 503.562201  3.73    3.77
#> 22     22 M6-F2.(F2L2)*F1.(F1L2) 86761.81 84713.97 524.410403  3.83    3.81
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
#> 16    3.75     3.73
#> 17    3.75     3.74
#> 18    3.79     3.79
#> 19    3.74     3.72
#> 20    3.75     3.74
#> 21    3.74     3.73
#> 22    3.80     3.80
```

After that we can run the function
[`inla.null.sp()`](https://vdak1ng.github.io/INLA.SocialEp/reference/inla.null.sp.md).
In this case, we will set the threshold for considering a spatial effect
as spurious if the standard deviation is below $0.125$:

``` r
ResMod <- INLA.SocialEp::inla.null.sp(ResMod, thres = 0.125)
```

And we can finally check the summary of the results obtained:

``` r
ResMod$Summary
#>    NUMBER                  MODEL      DIC     WAIC        CPU LOOCV LOGCV.3
#> 1       1                     M0 92182.37 89865.26   4.412914  4.56    4.56
#> 2       2                     M1 87736.16 85749.54  68.748523  3.93    3.91
#> 3       3      M2-ind(F1L1-F2L1) 91115.09 88700.49  30.726008  4.38    4.39
#> 4       4      M2-ind(F1L2-F2L1) 91116.86 88701.68  22.844925  4.38    4.39
#> 5       5      M2-ind(F2L1-F1L2) 91117.28 88701.80  29.020304  4.38    4.39
#> 6       6      M2-ind(F2L2-F1L2) 91117.44 88702.66  22.667647  4.38    4.39
#> 7       7           M3-F1.(F1L1) 90446.04 88097.93  71.984813  4.27    4.27
#> 8       8           M3-F1.(F1L2) 89684.16 87308.31  71.501817  4.16    4.17
#> 9       9           M4-F2.(F2L1) 90521.50 88118.28  86.686704  4.28    4.28
#> 10     10           M4-F2.(F2L2) 89761.28 87407.14  72.742198  4.16    4.18
#> 11     11 M5-F1.(F1L1)+F2.(F2L1) 88393.44 86095.40 208.292671  4.03    4.02
#> 12     12 M5-F1.(F1L2)+F2.(F2L1) 88102.19 85848.25 170.653224  3.99    3.99
#> 13     13 M5-F1.(F1L1)+F2.(F2L2) 88127.10 85861.30 211.727139  3.99    4.00
#> 14     14 M5-F1.(F1L2)+F2.(F2L2) 88975.18 86618.79 246.437363  4.07    4.07
#> 15     15 M6-F1.(F1L1)*F2.(F2L1) 86041.99 83972.57 402.481428  3.75    3.76
#> 16     16 M6-F1.(F1L2)*F2.(F2L1) 85779.43 83779.75 339.454437  3.72    3.77
#> 17     17 M6-F1.(F1L1)*F2.(F2L2) 85888.95 83943.00 533.309143  3.73    3.78
#> 18     18 M6-F1.(F1L2)*F2.(F2L2) 86826.41 84816.98 667.051420  3.84    3.80
#> 19     19 M6-F2.(F2L1)*F1.(F1L1) 86043.02 83967.15 293.247859  3.75    3.76
#> 20     20 M6-F2.(F2L2)*F1.(F1L1) 85788.90 83792.87 502.294037  3.72    3.77
#> 21     21 M6-F2.(F2L1)*F1.(F1L2) 85822.39 83845.67 503.562201  3.73    3.77
#> 22     22 M6-F2.(F2L2)*F1.(F1L2) 86761.81 84713.97 524.410403  3.83    3.81
#>    LOGCV.5 LOGCV.10 sp.null
#> 1     4.56     4.56       0
#> 2     3.91     3.91       0
#> 3     4.39     4.39       0
#> 4     4.39     4.39       0
#> 5     4.39     4.39       0
#> 6     4.39     4.39       0
#> 7     4.27     4.27       0
#> 8     4.17     4.17       0
#> 9     4.28     4.28       0
#> 10    4.17     4.17       0
#> 11    4.01     4.01       0
#> 12    3.99     3.99       0
#> 13    3.99     3.99       0
#> 14    4.06     4.06       0
#> 15    3.73     3.72       0
#> 16    3.75     3.73       0
#> 17    3.75     3.74       0
#> 18    3.79     3.79       0
#> 19    3.74     3.72       0
#> 20    3.75     3.74       0
#> 21    3.74     3.73       0
#> 22    3.80     3.80       0
```

We can also plot the spatial effects adjusted for each of the models
using the function \`plot.SpANOVA()\`\`:

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
INLA.SocialEp::plot.SpANOVA(
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

------------------------------------------------------------------------

## On this page

Developed by [Pablo Escobar-Hernández](https://github.com/VdaK1NG),
[Francisco Palmí-Perales](https://github.com/FranciscoPalmiPerales),
[Antonio López-Quílez](https://github.com/antoloqui).

Site built with [pkgdown](https://pkgdown.r-lib.org/) 2.2.0.
