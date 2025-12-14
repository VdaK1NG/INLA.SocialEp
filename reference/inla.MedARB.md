# Calculate MedARB between fitted values in INLA and simulated values

This function estimates Median Absolute Relative Bias between fitted
values in INLA and simulated values

## Usage

``` r
inla.MedARB(fit_values, sim_values, n.sim = 1)
```

## Arguments

- fit_values:

  Mean fitted values from INLA models

- sim_values:

  Vector of simulated values

- n.sim:

  Number of simulated datasets

## Value

Value estimated for MedARB
