# Low-Level Polygon Builder Used by `get_range()`

Internal helper that creates buffered polygons from clustered occurrence
points inside a single ecoregion.

## Usage

``` r
conv_function(sp_coord, bwp, bipl, bwpo, temp_dir, g = NULL)
```

## Arguments

- sp_coord:

  Spatial coordinates for one cluster of occurrences.

- bwp:

  Numeric buffer width, in degrees, used for isolated points.

- bipl:

  Numeric increment, in degrees, added to the point buffer for linear
  clusters.

- bwpo:

  Numeric buffer width, in degrees, applied to convex hull polygons.

- temp_dir:

  Character string giving the temporary directory used for intermediate
  convex-hull files.

- g:

  Optional ecoregion identifier used in status messages.
