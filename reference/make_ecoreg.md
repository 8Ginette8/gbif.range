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
version 2.1.2. <https://CRAN.R-project.org/package=cluster>

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
rst.path <- paste0(
  system.file(package = "gbif.range"),
  "/extdata/rst_enl.tif"
)
rst <- terra::rast(rst.path)
shp.path <- paste0(
  system.file(package = "gbif.range"),
  "/extdata/shp_lonlat.shp"
)
shp.lonlat <- terra::vect(shp.path)
rst <- terra::crop(rst, shp.lonlat)

# Apply the function by inferring 200 environmental classes
my.eco <- make_ecoreg(env = rst,
  nclass = 200,
  format = "sf"
)
#> CLARA algorithm processing... 
#> Generating polygons... 

# \donttest{
# Downloading in the European Alps the observations of one plant species
obs.paed <- get_gbif(
  sp_name = "Paederota bonarota",
  geo = shp.lonlat,
  grain = 1
)
#> |--------------------------------------------|
#> | Total number (all records)    :       1011 |
#> | Kept records                  :        585 |
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
#>          Grain filtering      78       507
#>       Duplicated records      12       495
#>          Absence records       0       495
#>          Basis selection      49       446
#>  Establishment selection       0       446
#>               Time frame       0       446
#>        Identical records       0       446
#>         Raster centroids       0       446
#> 
#> Initial records         : 585
#> Total removed           : 139
#> Final records (XY)      : 446
#> ----------------------------------------------
#> Final records (no XY)   : 0

# Create the range map based on:
# - custom ecoregion at 5 x 5 km resolution
# - smaller buffer because of regional extent
range.arcto <- get_range(
  occ_coord = obs.paed,
  ecoreg = my.eco,
  ecoreg_name = "EcoRegion",
  res = 0.05,
  degrees_outlier = 0.5,
  buff_width_point = 0.5,
  buff_incrmt_pts_line = 0.5,
  buff_width_polygon = 0.5
)
#> ## Start of computation for species:  Paederota bonarota  ### 
#> 4 outlier's from 437 | proportion from total points: 1%
#> ecoregion 1  of  71 :  100 
#> ecoregion 2  of  71 :  107 
#> ecoregion 3  of  71 :  108 
#> ecoregion 4  of  71 :  109 
#> ecoregion 5  of  71 :  110 
#> ecoregion 6  of  71 :  113 
#> ecoregion 7  of  71 :  114 
#> ecoregion 8  of  71 :  119 
#> ecoregion 9  of  71 :  12 
#> ecoregion 10  of  71 :  121 
#> ecoregion 11  of  71 :  122 
#> ecoregion 12  of  71 :  127 
#> ecoregion 13  of  71 :  128 
#> ecoregion 14  of  71 :  129 
#> ecoregion 15  of  71 :  130 
#> ecoregion 16  of  71 :  131 
#> ecoregion 17  of  71 :  132 
#> ecoregion 18  of  71 :  133 
#> ecoregion 19  of  71 :  138 
#> ecoregion 20  of  71 :  139 
#> ecoregion 21  of  71 :  14 
#> ecoregion 22  of  71 :  140 
#> ecoregion 23  of  71 :  142 
#> ecoregion 24  of  71 :  145 
#> ecoregion 25  of  71 :  147 
#> ecoregion 26  of  71 :  148 
#> ecoregion 27  of  71 :  15 
#> ecoregion 28  of  71 :  151 
#> ecoregion 29  of  71 :  153 
#> ecoregion 30  of  71 :  155 
#> ecoregion 31  of  71 :  156 
#> ecoregion 32  of  71 :  157 
#> ecoregion 33  of  71 :  16 
#> ecoregion 34  of  71 :  17 
#> ecoregion 35  of  71 :  19 
#> ecoregion 36  of  71 :  2 
#> ecoregion 37  of  71 :  22 
#> ecoregion 38  of  71 :  28 
#> ecoregion 39  of  71 :  29 
#> ecoregion 40  of  71 :  30 
#> ecoregion 41  of  71 :  31 
#> ecoregion 42  of  71 :  32 
#> ecoregion 43  of  71 :  33 
#> ecoregion 44  of  71 :  35 
#> ecoregion 45  of  71 :  40 
#> ecoregion 46  of  71 :  41 
#> ecoregion 47  of  71 :  44 
#> ecoregion 48  of  71 :  5 
#> ecoregion 49  of  71 :  57 
#> ecoregion 50  of  71 :  59 
#> ecoregion 51  of  71 :  60 
#> ecoregion 52  of  71 :  62 
#> ecoregion 53  of  71 :  64 
#> ecoregion 54  of  71 :  68 
#> ecoregion 55  of  71 :  69 
#> ecoregion 56  of  71 :  70 
#> ecoregion 57  of  71 :  71 
#> ecoregion 58  of  71 :  72 
#> ecoregion 59  of  71 :  77 
#> ecoregion 60  of  71 :  79 
#> ecoregion 61  of  71 :  8 
#> ecoregion 62  of  71 :  81 
#> ecoregion 63  of  71 :  84 
#> ecoregion 64  of  71 :  85 
#> ecoregion 65  of  71 :  86 
#> ecoregion 66  of  71 :  9 
#> ecoregion 67  of  71 :  92 
#> ecoregion 68  of  71 :  94 
#> ecoregion 69  of  71 :  95 
#> ecoregion 70  of  71 :  98 
#> ecoregion 71  of  71 :  99 
#> ## End of computation for species:  Paederota bonarota  ### 

# Plot
countries <- terra::vect(
  system.file("extdata", "world_countries.shp", package = "gbif.range")
)
alps.shp <- terra::aggregate(terra::crop(countries,terra::ext(rst)))
r.arcto <- terra::mask(range.arcto$rangeOutput,alps.shp)
terra::plot(terra::crop(countries,terra::ext(rst)), col = "#bcbddc")
terra::plot(
  r.arcto,
  add = TRUE,
  col = "darkgreen",
  axes = FALSE,
  legend = FALSE
)
graphics::points(
  obs.paed[, c("decimalLongitude","decimalLatitude")],
  pch = 20,
  col = "#99340470",
  cex = 1
)


# }
```
