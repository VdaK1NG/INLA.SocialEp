# Extract null effects from SR-ANOVA models

This function extracts the total number of null spatial effects inside a
SR-ANOVA object, for each one of the different specifications included.

## Usage

``` r
inla.null.sp(mod, thres = 0.125)
```

## Arguments

- mod:

  Mod file coming from a SR-ANOVA object

- thres:

  Threshold of sd use for the spatial effect to consider it null

## Value

SR-ANOVA object with the summary table updated to include null spatial
effects
