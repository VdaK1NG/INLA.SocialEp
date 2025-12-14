# Apply SpANOVA modelization using a wrapper function in INLA

This function is a wrapper for multidimensional spatial factor models in
INLA, using a sequential shared spatial effects with nested structure as
discussed in AÑADIR REFERENCIA.  
  
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
values are O3 and the last n values are O4

## Usage

``` r
inla.SpANOVA.2x2(
  obs,
  exp,
  gr,
  fac.names = NULL,
  lev.fac1 = NULL,
  lev.fac2 = NULL,
  scale.mod = TRUE,
  sp.prior = "sdunif",
  pc.prec.val = c(1, 0.01),
  sp.copy.fixed = TRUE,
  save.res = TRUE,
  save.random = TRUE,
  save.hyper = TRUE
)
```

## Arguments

- obs:

  Vector of observed values

- exp:

  Vector of expected values

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

## Value

List with all the models analyzed and a summary table with the most
common performance metrics.
