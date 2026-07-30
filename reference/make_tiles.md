# Create Tiled GBIF Geometry Queries

Divide a study extent into a set of smaller `POLYGON()` query strings
that can be used with the GBIF `geometry` parameter.

## Usage

``` r
make_tiles(geo, ntiles, sext = TRUE)
```

## Arguments

- geo:

  Spatial extent or geometry. Used to define the study area. Accepted
  classes are `Extent`, `SpatExtent`, `SpatialPolygon`,
  `SpatialPolygonDataFrame`, `SpatVector`, and `sf`. If `NULL`, the
  whole globe is used.

- ntiles:

  Numeric. Approximate number of tiles to create.

- sext:

  Logical. Should the corresponding `SpatExtent` objects also be
  returned?

## Value

If `sext = TRUE`, a list containing WKT `POLYGON()` strings and matching
`SpatExtent` objects. Otherwise, a vector of WKT `POLYGON()` strings.

## Details

The input extent is converted into a regular grid of smaller query
polygons. This is mainly intended to support tiled GBIF downloads when
large extents would otherwise return too many records in a single query.

## References

Chauvier, Y., Thuiller, W., Brun, P., Lavergne, S., Descombes, P.,
Karger, D. N., ... & Zimmermann, N. E. (2021). Influence of climate,
soil, and land cover on plant species distribution in the European Alps.
Ecological Monographs, 91(2), e01433.
[doi:10.1002/ecm.1433](https://doi.org/10.1002/ecm.1433)

## See also

[`get_gbif`](https://8ginette8.github.io/gbif.range/reference/get_gbif.md)()
which uses this tiling internally for large extents.

## Examples

``` r
# Load the European Alps Extent
shp.path <- paste0(
  system.file(package = "gbif.range"), 
  "/extdata/shp_lonlat.shp"
)
shp.lonlat <- terra::vect(shp.path)
 
# Apply the function to divide the extent in ~20 fragments
mt <- make_tiles(geo = shp.lonlat, ntiles = 20, sext = TRUE)
```
