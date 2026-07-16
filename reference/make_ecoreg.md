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
###########################################
### Example plot
###########################################

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
obs.arcto <- get_gbif(
  sp_name = "Arctostaphylos alpinus",
  geo = shp.lonlat,
  grain = 1
)
#> |--------------------------------------------|
#> | Total number (all records)    :      42972 |
#> | Kept records                  :       6340 |
#> |--------------------------------------------|
#> | Kept records according to parameters:
#> | spatial_issue = FALSE, has_xy = TRUE by default ('geo' was set)
#> 
#> ...GBIF records of Arctostaphylos alpinus: download starting...
#> ------------- #1 (100%..)               
#> 
#> ...Records (XY) filtering summary:
#> -----------------------------------------------
#>                     step removed remaining
#>          Grain filtering    4968      1372
#>       Duplicated records     448       924
#>          Absence records       0       924
#>          Basis selection      87       837
#>  Establishment selection       0       837
#>               Time frame       0       837
#>        Identical records       0       837
#>         Raster centroids       0       837
#> 
#> Initial records         : 6340
#> Total removed           : 5503
#> Final records (XY)      : 837
#> -----------------------------------------------
#> Final records (no XY)   : 0

# Create the range map based on our custom ecoregion at 5 x 5 km resolution
range.arcto <- get_range(
  occ_coord = obs.arcto,
  ecoreg = my.eco,
  ecoreg_name = "EcoRegion",
  res = 0.05
)
#> ## Start of computation for species:  Arctostaphylos alpinus  ### 
#> 0 outlier's from 807 | proportion from total points: 0%
#> ecoregion 1  of  83 :  101 
#> ecoregion 2  of  83 :  102 
#> ecoregion 3  of  83 :  103 
#> ecoregion 4  of  83 :  104 
#> ecoregion 5  of  83 :  105 
#> ecoregion 6  of  83 :  106 
#> ecoregion 7  of  83 :  107 
#> ecoregion 8  of  83 :  108 
#> ecoregion 9  of  83 :  109 
#> ecoregion 10  of  83 :  110 
#> ecoregion 11  of  83 :  111 
#> ecoregion 12  of  83 :  112 
#> ecoregion 13  of  83 :  113 
#> ecoregion 14  of  83 :  114 
#> ecoregion 15  of  83 :  115 
#> ecoregion 16  of  83 :  116 
#> ecoregion 17  of  83 :  117 
#> ecoregion 18  of  83 :  118 
#> ecoregion 19  of  83 :  119 
#> ecoregion 20  of  83 :  120 
#> ecoregion 21  of  83 :  121 
#> ecoregion 22  of  83 :  122 
#> ecoregion 23  of  83 :  125 
#> [ecoreg = 23] 4 points lying on one line. Using buffer width of 611.1105 km
#> ecoregion 24  of  83 :  126 
#> ecoregion 25  of  83 :  127 
#> ecoregion 26  of  83 :  128 
#> ecoregion 27  of  83 :  129 
#> ecoregion 28  of  83 :  130 
#> ecoregion 29  of  83 :  131 
#> ecoregion 30  of  83 :  132 
#> ecoregion 31  of  83 :  133 
#> ecoregion 32  of  83 :  136 
#> ecoregion 33  of  83 :  138 
#> ecoregion 34  of  83 :  139 
#> ecoregion 35  of  83 :  140 
#> ecoregion 36  of  83 :  142 
#> ecoregion 37  of  83 :  143 
#> ecoregion 38  of  83 :  145 
#> ecoregion 39  of  83 :  146 
#> ecoregion 40  of  83 :  147 
#> ecoregion 41  of  83 :  148 
#> ecoregion 42  of  83 :  151 
#> ecoregion 43  of  83 :  152 
#> ecoregion 44  of  83 :  153 
#> ecoregion 45  of  83 :  155 
#> ecoregion 46  of  83 :  156 
#> ecoregion 47  of  83 :  16 
#> ecoregion 48  of  83 :  191 
#> ecoregion 49  of  83 :  27 
#> ecoregion 50  of  83 :  28 
#> ecoregion 51  of  83 :  29 
#> ecoregion 52  of  83 :  30 
#> ecoregion 53  of  83 :  36 
#> ecoregion 54  of  83 :  40 
#> ecoregion 55  of  83 :  41 
#> ecoregion 56  of  83 :  44 
#> ecoregion 57  of  83 :  56 
#> ecoregion 58  of  83 :  57 
#> ecoregion 59  of  83 :  58 
#> ecoregion 60  of  83 :  59 
#> ecoregion 61  of  83 :  62 
#> ecoregion 62  of  83 :  67 
#> ecoregion 63  of  83 :  68 
#> ecoregion 64  of  83 :  69 
#> ecoregion 65  of  83 :  71 
#> ecoregion 66  of  83 :  72 
#> ecoregion 67  of  83 :  74 
#> ecoregion 68  of  83 :  75 
#> ecoregion 69  of  83 :  76 
#> ecoregion 70  of  83 :  77 
#> ecoregion 71  of  83 :  84 
#> ecoregion 72  of  83 :  85 
#> ecoregion 73  of  83 :  86 
#> ecoregion 74  of  83 :  88 
#> ecoregion 75  of  83 :  89 
#> ecoregion 76  of  83 :  92 
#> ecoregion 77  of  83 :  93 
#> ecoregion 78  of  83 :  94 
#> ecoregion 79  of  83 :  95 
#> ecoregion 80  of  83 :  96 
#> ecoregion 81  of  83 :  97 
#> ecoregion 82  of  83 :  98 
#> ecoregion 83  of  83 :  99 
#> ## End of computation for species:  Arctostaphylos alpinus  ### 

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
  obs.arcto[, c("decimalLongitude","decimalLatitude")],
  pch = 20,
  col = "#99340470",
  cex = 1
)


###########################################
### Package manuscript plot (Fig 1a, left)
###########################################

# Root and package
root_dir <- tempdir()
if (!dir.exists(file.path(root_dir, "fig_plots"))) {
    dir.create(file.path(root_dir, "fig_plots"))
}

# Assign colors to ecoregions
my.eco = terra::vect(my.eco)
eco.alps <- terra::crop(my.eco, terra::ext(shp.lonlat))
eco.alps2 <- terra::intersect(eco.alps, alps.shp)
col.palette <- grDevices::colorRampPalette(c("#a6cee3", "#1f78b4",
  "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c",
  "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a",
  "#ffff99", "#b15928"))
colcol <- col.palette(length(my.eco))
set.seed(7)
eco.alps2$color <- sample(
  paste0(colcol, ""),
  length(eco.alps2),
  replace = FALSE
)

# Extract ecoregions values for points
pt.col <- terra::extract(
  x = eco.alps2,
  y = as.data.frame(obs.arcto[, c("decimalLongitude","decimalLatitude")])
)
pt.plot <- obs.arcto[!is.na(pt.col$color),
                      c("decimalLongitude","decimalLatitude")]
pt.col2 <- pt.col[!is.na(pt.col$color), "color"]
pt.col3 <- grDevices::adjustcolor(pt.col2, red.f = 0.6, green.f = 0.6, blue.f = 0.6)

# Plot
grDevices::png(
  paste0(
    root_dir,
    "/fig_plots/fig1_arcto.png"
  ),
  width = 100,
  height = 70,
  unit = "cm",
  res = 100,
  pointsize = 110
)
oldpar <- graphics::par(mfrow = c(1,1), mar = c(5,5,5,20), lwd = 10, cex = 1)
terra::plot(
  eco.alps2,
  col = paste0(eco.alps2$color, "99"),
  border = NA,
  axes = FALSE
)
terra::plot(
  terra::as.polygons(r.arcto),
  border = "black",
  lwd = 7,
  col = "#00000099",
  add = TRUE
)
graphics::points(pt.plot, col = pt.col2, pch = 16, cex = 0.7)
graphics::points(pt.plot, col = pt.col3, pch = 16, cex = 0.4)
grDevices::dev.off()
#> agg_record_1ded7b32c317 
#>                       2 
graphics::par(oldpar)

# }
```
