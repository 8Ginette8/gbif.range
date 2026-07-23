# Filter GBIF Records by Grid Cell

Reduce the spatial density of a `getGBIF` object by retaining a single
record per grid cell and, optionally, removing cells with too few
records.

## Usage

``` r
obs_filter(gbifs, grid, threshold = NULL)
```

## Arguments

- gbifs:

  `getGBIF` object containing one or more species.

- grid:

  Raster. Defining the target spatial resolution and extent. Accepted
  classes are `SpatRaster`, `RasterLayer`, `RasterBrick`, and
  `RasterStack`.

- threshold:

  Optional integer. Specifying the minimum number of records a cell must
  contain to be retained.

## Value

A data frame with the columns `Species`, `x`, and `y`, representing the
filtered coordinates.

## Details

The function first collapses each species to one occurrence per grid
cell. If `threshold` is supplied, cells with fewer than that many
original records are then discarded.

## See also

[`get_gbif`](https://8ginette8.github.io/gbif.range/reference/get_gbif.md)()
to produce the `getGBIF` object filtered by this function.

## Examples

``` r
# \donttest{
# Load data
shp.path <- paste0(
  system.file(package = "gbif.range"),
  "/extdata/shp_lonlat.shp"
)
shp.lonlat <- terra::vect(shp.path)
rst.path <- paste0(
  system.file(package = "gbif.range"),
  "/extdata/rst_enl.tif"
)
rst <- terra::rast(rst.path)

# Download observations for two plant species in the European Alps
obs.paed <- get_gbif(
  sp_name = "Paederota bonarota",
  geo = shp.lonlat
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
#>          Grain filtering       0       585
#>       Duplicated records      31       554
#>          Absence records       0       554
#>          Basis selection      74       480
#>  Establishment selection       0       480
#>               Time frame       0       480
#>        Identical records       0       480
#>         Raster centroids       0       480
#> 
#> Initial records         : 585
#> Total removed           : 105
#> Final records (XY)      : 480
#> ----------------------------------------------
#> Final records (no XY)   : 0
obs.saxi <- get_gbif(
  sp_name = "Saxifraga cernua",
  geo = shp.lonlat
)
#> |--------------------------------------------|
#> | Total number (all records)    :      20266 |
#> | Kept records                  :        407 |
#> |--------------------------------------------|
#> | Kept records according to parameters:
#> | spatial_issue = FALSE, has_xy = TRUE by default ('geo' was set)
#> 
#> ...GBIF records of Saxifraga cernua: download starting...
#> ------------- #1 (100%..)               
#> 
#> ...Records (XY) filtering summary:
#> ----------------------------------------------
#>                     step removed remaining
#>          Grain filtering       5       402
#>       Duplicated records     286       116
#>          Absence records       0       116
#>          Basis selection      62        54
#>  Establishment selection       0        54
#>               Time frame       0        54
#>        Identical records       0        54
#>         Raster centroids       0        54
#> 
#> Initial records         : 407
#> Total removed           : 353
#> Final records (XY)      : 54
#> ----------------------------------------------
#> Final records (no XY)   : 0

# Test plot
terra::plot(shp.lonlat)
graphics::points(
  obs.paed[, c("decimalLongitude","decimalLatitude")],
  pch = 20,
  col = "#238b4550",
  cex = 1
)
graphics::points(
  obs.saxi[, c("decimalLongitude","decimalLatitude")],
  pch = 20,
  col = "#99000d50",
  cex = 1
)


# Combine both datasets
both.sp <- rbind(obs.paed,obs.saxi)

# Run function
obs.filt <- obs_filter(gbifs = both.sp, grid = rst, threshold = 4)

# Check new points
terra::plot(shp.lonlat)
graphics::points(
  obs.filt[obs.filt$Species%in%"Arctostaphylos alpinus",c("x","y")],
  pch = 20,
  col = "#238b4550",
  cex = 1
)
graphics::points(
  obs.filt[obs.filt$Species%in%"Saxifraga cernua",c("x","y")],
  pch = 20,
  col = "#99000d50",
  cex = 1
)


# }
```
