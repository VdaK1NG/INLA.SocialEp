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

Moreover, the over user needs to provide the data with a data.frame that
has the following structure:

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
#> 1       1                     M0 92182.37 89865.26   4.815679  4.56    4.56
#> 2       2                     M1 87736.21 85749.37  69.937367  3.93    3.91
#> 3       3      M2-ind(F1L1-F2L1) 91115.10 88700.47  29.684855  4.38    4.39
#> 4       4      M2-ind(F1L2-F2L1) 91116.86 88701.68  22.368629  4.38    4.39
#> 5       5      M2-ind(F2L1-F1L2) 91117.07 88701.94  27.600629  4.38    4.39
#> 6       6      M2-ind(F2L2-F1L2) 91117.41 88702.72  22.144730  4.38    4.39
#> 7       7           M3-F1.(F1L1) 90446.04 88097.93  69.985486  4.27    4.27
#> 8       8           M3-F1.(F1L2) 89684.16 87308.31  68.871715  4.16    4.17
#> 9       9           M4-F2.(F2L1) 90521.68 88118.20  84.434294  4.28    4.28
#> 10     10           M4-F2.(F2L2) 89761.28 87407.14  69.950024  4.16    4.18
#> 11     11 M5-F1.(F1L1)+F2.(F2L1) 88392.49 86096.81 198.738812  4.03    4.02
#> 12     12 M5-F1.(F1L2)+F2.(F2L1) 88102.19 85848.25 165.246119  3.99    3.99
#> 13     13 M5-F1.(F1L1)+F2.(F2L2) 88127.04 85861.26 208.429858  3.99    4.00
#> 14     14 M5-F1.(F1L2)+F2.(F2L2) 88975.27 86618.26 233.116748  4.07    4.07
#> 15     15 M6-F1.(F1L1)*F2.(F2L1) 86041.62 83972.66 424.414256  3.75    3.76
#> 16     16 M6-F1.(F1L2)*F2.(F2L1) 85782.34 83779.01 345.212833  3.72    3.77
#> 17     17 M6-F1.(F1L1)*F2.(F2L2) 85900.47 83931.71 377.189164  3.73    3.78
#> 18     18 M6-F1.(F1L2)*F2.(F2L2) 86826.77 84816.29 490.866777  3.84    3.80
#> 19     19 M6-F2.(F2L1)*F1.(F1L1) 86043.52 83971.35 314.585944  3.75    3.76
#> 20     20 M6-F2.(F2L2)*F1.(F1L1) 85791.01 83791.23 464.118779  3.72    3.77
#> 21     21 M6-F2.(F2L1)*F1.(F1L2) 85821.93 83846.91 673.083046  3.73    3.77
#> 22     22 M6-F2.(F2L2)*F1.(F1L2) 86762.86 84714.13 455.947865  3.83    3.81
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
ResMod <- inla.null.sp(ResMod, thres = 0.125)
```

And we can finally check the summary of the results obtained:

``` r
ResMod$Summary
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

------------------------------------------------------------------------

## On this page

Developed by [Pablo Escobar-Hernández](https://github.com/VdaK1NG),
[Francisco Palmí-Perales](https://github.com/FranciscoPalmiPerales),
[Antonio López-Quílez](https://github.com/antoloqui).

Site built with [pkgdown](https://pkgdown.r-lib.org/) 2.2.0.
