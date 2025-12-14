# Calculate MRRMSE between fitted values in INLA and simulated values

This function estimates Mean Relative Root Mean Prediction Error between
fitted values in INLA and simulated values

## Usage

``` r
inla.MRRMSE(fit_values, sim_values, n.sim = 1)
```

## Arguments

- fit_values:

  Mean fitted values from INLA models

- sim_values:

  Vector of simulated values

- n.sim:

  Number of simulated datasets

## Value

Value estimated for MRRMSE
