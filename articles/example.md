# Usage Example

## Introduction

``` r

library("INLA.SocialEp")
```

## Understanding the inputs

First of all, the user needs to prepare a graph from an sf object using
the following functions from packages *spdep* and *INLA*:

``` r

sec.nb <- spdep::poly2nb(sp_object_sim, snap=0.0000005)
spdep::nb2INLA("MapGraph", sec.nb)
graph_states <- INLA::inla.read.graph("MapGraph")
```

Next, the user needs to make sure the data has the proper format.
Specifically, data for the model needs to be a data.frame that has the
following structure:

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

Moreover, the user needs to specify the following arguments:

- **gr**: graph for the underlying spatial structure generated using the
  function `inla.read.graph()` from the INLA package.
- **fac.names**: name of the factors, with the format `c("F1", "F2")`.
- **lev.fac1**: name of the levels of factor 1, with the format
  `c("F1L1", "F1L2")`.
- **lev.fac2**: name of the levels of factor 2, with the format
  `c("F2L1", "F2L2")`.
- **thres**: threshold used to consider a spatial effect as spurious,
  default is $`0.125`$.
- **verbose.INLA**: option to run INLA with `verbose=TRUE` if the user
  wants to do debugging on the models, default is *FALSE*.

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
  spatial effects have their variance scaled (`scale.mod=TRUE`) or not
  (`scale.mod=FALSE`).

Lastly, the user can specify several options that determine which parts
of the INLA models are saved in the final object. Particularly, the user
has the following options:

- **save.res**: save residuals from the different models, default is
  *FALSE*.
- **save.random**: save random spatial effects from the different
  models, default is *TRUE*.
- **save.fixed**: save fixed effects from the different models, default
  is *TRUE*.
- **save.hyper**: save hyperparameters data from the different models,
  default is *TRUE*.
- **save.mod.data**: save data to rerun the model later using INLA
  function directly, default is *FALSE*.

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

The object generated by the function contains a **Summary** with the
results of the different models, which can be easily accessed by doing
`object$Summary`:

``` r

ResMod$Summary
#>    NUMBER                          MODEL      DIC     WAIC    CPU LOOCV LOGCV.3
#> 1       1                             M0 92182.37 89865.26   4.17  4.56    4.56
#> 2       2                             M1 87736.15 85749.55  65.88  3.93    3.91
#> 3       3                 M2.(F1L1-F2L1) 91115.09 88700.48  29.20  4.38    4.39
#> 4       4                 M2.(F1L2-F2L1) 91116.86 88701.68  21.94  4.38    4.39
#> 5       5                 M2.(F1L1-F2L2) 91117.27 88701.90  28.60  4.38    4.39
#> 6       6                 M2.(F1L2-F2L2) 91117.42 88702.70  22.91  4.38    4.39
#> 7       7               M3.F1(ref: F1L1) 90446.04 88097.93  67.42  4.27    4.27
#> 8       8               M3.F1(ref: F1L2) 89684.16 87308.31  68.68  4.16    4.17
#> 9       9               M4.F2(ref: F2L1) 90521.49 88118.37  82.84  4.28    4.28
#> 10     10               M4.F2(ref: F2L2) 89761.28 87407.14  69.65  4.16    4.18
#> 11     11 M5.F1(ref: F1L1)+F2(ref: F2L1) 88392.88 86095.64 201.43  4.03    4.02
#> 12     12 M5.F1(ref: F1L2)+F2(ref: F2L1) 88102.19 85848.25 165.70  3.99    3.99
#> 13     13 M5.F1(ref: F1L1)+F2(ref: F2L2) 88127.09 85861.29 203.82  3.99    4.00
#> 14     14 M5.F1(ref: F1L2)+F2(ref: F2L2) 88975.64 86618.34 234.53  4.07    4.07
#> 15     15 M6.F1(ref: F1L1)*F2(ref: F2L1) 86042.16 83972.88 458.90  3.75    3.76
#> 16     16 M6.F1(ref: F1L2)*F2(ref: F2L1) 85780.50 83779.56 368.90  3.72    3.77
#> 17     17 M6.F1(ref: F1L1)*F2(ref: F2L2) 85891.09 83944.68 492.61  3.73    3.78
#> 18     18 M6.F1(ref: F1L2)*F2(ref: F2L2) 86826.87 84817.63 475.04  3.84    3.80
#> 19     19 M6.F2(ref: F2L1)*F1(ref: F1L1) 86044.38 83969.84 302.31  3.75    3.76
#> 20     20 M6.F2(ref: F2L2)*F1(ref: F1L1) 85788.94 83789.82 575.16  3.72    3.77
#> 21     21 M6.F2(ref: F2L1)*F1(ref: F1L2) 85822.97 83846.20 451.67  3.73    3.77
#> 22     22 M6.F2(ref: F2L2)*F1(ref: F1L2) 86762.06 84713.91 427.50  3.83    3.81
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
#> 16    3.74     3.73       0
#> 17    3.75     3.74       0
#> 18    3.79     3.79       0
#> 19    3.74     3.72       0
#> 20    3.75     3.74       0
#> 21    3.74     3.73       0
#> 22    3.80     3.80       0
```

If the user chose the option `save.random=TRUE`, then the different
spatial effects adjusted can be accessed and plotted. Package includes
the
[`plotSpANOVA()`](https://vdak1ng.github.io/INLA.SocialEp/reference/plotSpANOVA.md)
function, allowing the user to easily plot both the spatial effects and
the relative risks adjusted. For example, we could plot the spatial
effect for model number 20 using the following code:

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

Moreover, if the user wanted to plot the relative risks adjusted for
that same model, we could use the following code:

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
