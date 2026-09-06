# Apply SpANOVA modelization using a wrapper function in INLA

This function is a wrapper for multidimensional spatial factor models in
INLA, using a sequential shared spatial effects with nested structure as
discussed in AÑADIR REFERENCIA.  

## Usage

``` r
inla.SpANOVA.2x2(
  data,
  gr,
  fac.names = NULL,
  lev.fac1 = NULL,
  lev.fac2 = NULL,
  scale.mod = TRUE,
  sp.prior = "sdunif",
  pc.prec.val = c(1, 0.01),
  sp.copy.fixed = TRUE,
  save.res = FALSE,
  save.random = TRUE,
  save.hyper = TRUE,
  save.fixed = TRUE,
  save.mod.data = FALSE,
  verbose.INLA = FALSE,
  save.cov.mat = TRUE,
  thres = 0.125
)
```

## Arguments

- data:

  Dataframe containing the number of events observed on column obs, the
  expected values on column exp, the level for factor 1 on column
  lev.fac1 and the level for factor 1 on column lev.fac2

- gr:

  Graph for the underlying spatial structure

- fac.names:

  Names of the factors included

- lev.fac1:

  Levels of the first factor included

- lev.fac2:

  Levels of the second factor included

- scale.mod:

  Scale copies of random spatial effects or not, default is TRUE

- sp.prior:

  Select prior for the random spatial effect, options are sdunif and
  pc.prec

- pc.prec.val:

  Define values por the pc prior in case it was chosen. Default values
  are c(1, 0.01)

- sp.copy.fixed:

  Fix copied values for the random spatial effects, default is TRUE

- save.res:

  Save fitted values or not from the different models, default is TRUE

- save.random:

  Save values adjusted from the random spatial effects or not from the
  different models, default is TRUE

- save.hyper:

  Save hyperparameter values from each individual model, default is TRUE

- save.fixed:

  Save values adjusted from the fixed effects or not from the different
  models, default is TRUE

- save.mod.data:

  Save modelling data to run the model afterwards

- verbose.INLA:

  Verbose option for INLA, default is FALSE

- save.cov.mat:

  Save posterior variance-covariance matrix of the random effects

- thres:

  thres value for inla.null.sp function

## Value

List with all the models analyzed and a summary table with the most
common performance metrics.
