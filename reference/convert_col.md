# Generate stable color labels

This function generates stable color labels to ease comparison when
representing several figures that we would want to have the same
palette.

## Usage

``` r
convert_col(
  data,
  breaks,
  pal_fun,
  include.lowest = TRUE,
  right = TRUE,
  na.col = NA
)
```

## Arguments

- data:

  Vector of numeric values to be convert into colors

- breaks:

  Vector of breaks to use for cutting

- pal_fun:

  Vector of colors for the palette

- include.lowest:

  Include.lowest for the cut function

- right:

  Right option for the cut function

- na.col:

  Colour to use in case of NAs

## Value

Returns a list including the colors generated (fill_by), the tags for
the intervals (tags) and the colors corresponding for each interval
(colors)
