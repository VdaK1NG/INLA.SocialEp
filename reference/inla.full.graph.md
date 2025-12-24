# Make a fully connected graph for INLA

This function converts a graph with islands into a fully connected graph
by searching for the nearest area between the different nodes and
connecting them, with the possibility of setting a threshold for
connecting several areas if they are within the limits.

## Usage

``` r
inla.full.graph(
  sp_obj,
  graph,
  thresh = 0.1,
  snap = 5e-07,
  distance_by = c("perimeter", "centroid")
)
```

## Arguments

- sp_obj:

  Spatial object in sp form

- graph:

  Graph object from nb2INLA function

- snap:

  Snap argument for the function poly2nb

- distance_by:

  Choices are perimeter and centroid

## Value

List including three different objects: the adjacency object, the edited
graph and a dataframe with the changes applied
