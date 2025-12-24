# Represent different values from SpANOVA modelling

This function represents the desired value from a SpANOVA object or an
INLA object after running inla.rerun.SpANOVA

## Usage

``` r
# S3 method for class 'SpANOVA'
plot(
  obj,
  obj_type = c("SpANOVA", "INLA"),
  fill_by = c("Spatial", "Heterogeneity", "RR"),
  n_mod,
  sp_object,
  breaks = NA,
  fil_scale = c("#133BF2", "#7189F7", "#FFFFFF", "#FF867A", "#FF2F1B"),
  col_frontiers = "black",
  scale_name = "Values",
  sp_null = 0.125,
  legend.position = "right",
  ncol_fig = 2
)
```

## Arguments

- obj:

  Object to generate the plot from

- obj_type:

  Type of object provided, options are SpANOVA or INLA after running
  inla.rerun.SpANOVA

- fill_by:

  Values to represent, choosing between Spatial, Heterogeneity and RR

- n_mod:

  Number of the model specification that the user wants to represent

- sp_object:

  Spatial object to be plotted

- breaks:

  Breaks for the color palette in case the user wants to modify the
  default ones

- fil_scale:

  Vector of colors in case the user wants to modify the default ones

- col_frontiers:

  Colour for the lines that define each polygon, default is black

- scale_name:

  Name for the filling scale of the plot

- sp_null:

  Threshold desired to considered a spatial effect represented null,
  default is 0.125

- legend.position:

  Position for the legend of the plot

- ncol_fig:

  Number of columns desired for the grid of figures

## Value

Representation of the values desired using ggplot2 as baseline
