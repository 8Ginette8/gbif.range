# Build Custom Ecoregions from Environmental Layers

Cluster multi-layer environmental data to create a custom ecoregion map
that can be used directly in
[`get_range()`](https://8ginette8.github.io/gbif.range/reference/get_range.md).

## Usage

``` r
make_ecoreg(
  env = NULL,
  nclass = NULL,
  path = "",
  name = "",
  format = c("SpatVector", "sf", "SpatRaster"),
  verbose = TRUE,
  ...
)
```

## Arguments

- env:

  Raster stack. Can be any abiotic or biotic factors thought of defining
  ecoregions boundaries. Accepted classes are `SpatRaster`,
  `RasterBrick`, and `RasterStack`.

- nclass:

  Numeric. Number of environmental classes to create.

- path:

  Optional character. Directory where the output should be written.
  Leave empty to return the result directly.

- name:

  Character. Output file name without extension when `path` is used.

- format:

  Output format. One of `"SpatVector"` (default), `"sf"`, or
  `"SpatRaster"`. `"SpatRaster"` returns the raw cluster raster instead
  of converting to polygons.

- verbose:

  Logical. Should progress messages be printed? Default is `TRUE`.

- ...:

  Additional arguments passed to
  [`cluster::clara()`](https://rdrr.io/pkg/cluster/man/clara.html).

## Value

If `path == ""`, returns the generated raster or polygon object.
Otherwise, writes the output to disk as a GeoTIFF or Shapefile.

## Details

This function is useful when the packaged ecoregion layers are too
coarse for a study area or when a custom environmental regionalization
is needed. Clusters are created with the CLARA algorithm on the
multivariate environmental space represented by `env`.

## References

Chauvier, Y., Zimmermann, N. E., Poggiato, G., Bystrova, D., Brun, P., &
Thuiller, W. (2021). Novel methods to correct for observer and sampling
bias in presence-only species distribution models. Global Ecology and
Biogeography, 30(11), 2312-2325.
[doi:10.1111/geb.13383](https://doi.org/10.1111/geb.13383)

Maechler, M., Rousseeuw, P., Struyf, A., Hubert, M., & Hornik, K.
(2021). cluster: Cluster Analysis Basics and Extensions. R package
version 2.1.2. <https://CRAN.R-project.org/package=cluster/>

Reynolds, A. P., Richards, G., de la Iglesia, B., & Rayward-Smith, V. J.
(2006). Clustering rules: A comparison of partitioning and hierarchical
clustering algorithms. Journal of Mathematical Modelling and Algorithms,
5(4), 475-504.
[doi:10.1007/s10852-005-9022-1](https://doi.org/10.1007/s10852-005-9022-1)

Schubert, E., & Rousseeuw, P. J. (2019). Faster k-Medoids clustering:
Improving the PAM, CLARA, and CLARANS algorithms. In G. Amato, C.
Gennaro, V. Oria, & M. Radovanović (Eds.), Similarity search and
applications. SISAP 2019. Lecture Notes in Computer Science (Vol. 11807,
pp. 171-187). Springer.

## See also

[`get_range`](https://8ginette8.github.io/gbif.range/reference/get_range.md)()
to build a range map using the ecoregion layer produced here.

## Examples

``` r
# Open data
rst_path <- paste0(
  system.file(package = "gbif.range"),
  "/extdata/rst_enl.tif"
)
rst <- terra::rast(rst_path)
shp_path <- paste0(
  system.file(package = "gbif.range"),
  "/extdata/shp_lonlat.shp"
)
shp_lonlat <- terra::vect(shp_path)
rst <- terra::crop(rst, shp_lonlat)

# Apply the function by inferring 50 environmental classes
my_eco <- make_ecoreg(env = rst,
  nclass = 50,
  format = "sf"
)
#> CLARA algorithm processing... 
#> Generating polygons... 

# \donttest{
# Downloading in the European Alps the observations of one plant species
obs_paed <- get_gbif(
  sp_name = "Paederota bonarota",
  geo = shp_lonlat,
  grain = 1
)
#> |--------------------------------------------|
#> | Total number (all records)    :       1044 |
#> | Kept records                  :        618 |
#> |--------------------------------------------|
#> | Kept records according to parameters:
#> | spatial_issue = FALSE, has_xy = TRUE by default ('geo' was set)
#> 
#> ...GBIF records of Paederota bonarota: download starting...
#> ------------- #1 (100%..)               
#> 
#> ...Records (XY) filtering summary:
#> ----------------------------------------------
#>                     step removed remaining
#>          Grain filtering      79       539
#>       Duplicated records      12       527
#>          Absence records       0       527
#>          Basis selection      49       478
#>  Establishment selection       0       478
#>               Time frame       0       478
#>        Identical records       0       478
#>         Raster centroids       0       478
#> 
#> Initial records         : 618
#> Total removed           : 140
#> Final records (XY)      : 478
#> ----------------------------------------------
#> Final records (no XY)   : 0

# Guard against an unavailable remote source
if (gbif_have(obs_paed)) {

# Create the range map based on:
# - custom ecoregion at 5 x 5 km resolution
# - smaller buffer because of regional extent
range_paed <- get_range(
  occ_coord = obs_paed,
  ecoreg = my_eco,
  ecoreg_name = "EcoRegion",
  res = 0.05,
  degrees_outlier = 0.5,
  buff_width_point = 0.5,
  buff_incrmt_pts_line = 0.5,
  buff_width_polygon = 0.5
)

# Plot
countries <- terra::vect(
  system.file("extdata", "world_countries.shp", package = "gbif.range")
)
terra::plot(terra::crop(countries,terra::ext(rst)), col = "#bcbddc")
terra::plot(
  merge_range(range_paed),
  add = TRUE,
  col = "darkgreen",
  axes = FALSE,
  legend = FALSE
)
graphics::points(
  obs_paed[, c("decimalLongitude","decimalLatitude")],
  pch = 20,
  col = "#99340470",
  cex = 1
)

}
#> ## Start of computation for species: Paederota bonarota ###
#> 4 outlier's from 466 | proportion from total points: 1%
#> ecoregion 1 of 27: 10
#> ecoregion 2 of 27: 11
#> ecoregion 3 of 27: 13
#> ecoregion 4 of 27: 16
#> ecoregion 5 of 27: 17
#> ecoregion 6 of 27: 18
#> ecoregion 7 of 27: 2
#> ecoregion 8 of 27: 20
#> ecoregion 9 of 27: 27
#> ecoregion 10 of 27: 28
#> ecoregion 11 of 27: 29
#> ecoregion 12 of 27: 3
#> ecoregion 13 of 27: 30
#> ecoregion 14 of 27: 33
#> ecoregion 15 of 27: 34
#> ecoregion 16 of 27: 35
#> ecoregion 17 of 27: 36
#> ecoregion 18 of 27: 37
#> ecoregion 19 of 27: 38
#> ecoregion 20 of 27: 39
#> ecoregion 21 of 27: 4
#> ecoregion 22 of 27: 42
#> ecoregion 23 of 27: 5
#> ecoregion 24 of 27: 6
#> ecoregion 25 of 27: 7
#> ecoregion 26 of 27: 8
#> ecoregion 27 of 27: 9
#> ## End of computation for species: Paederota bonarota ###

# }
```
