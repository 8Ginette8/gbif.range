# Split Data into Approximately Balanced Folds

Create a fold-assignment vector for random or spatially structured
cross-validation.

## Usage

``` r
make_blocks(
  nfolds = 5,
  df = data.frame(),
  nblocks = nfolds * 2,
  npoints = NA,
  pres = numeric()
)
```

## Arguments

- nfolds:

  Numeric. Number of folds to create.

- df:

  Optional `data.frame`. Columns define the structure used to build
  folds. When omitted, random folds are created from `npoints`.

- nblocks:

  Numeric. Number of initial clusters used to build the folds. Must be
  at least `nfolds`.

- npoints:

  Optional numeric. Number of observations to split when `df` is not
  supplied.

- pres:

  Optional binary vector. Used to restrict CLARA clustering to a subset
  of rows and assign the remainder by nearest neighbours.

## Value

An integer vector of length `nrow(df)` or `npoints`, giving the fold
assignment for each observation.

## Details

If `df` has one column, folds are based on quantile bins. If `df` has
two or more columns, folds are based on CLARA clustering. Remaining
clusters are assigned to folds with an optimization step that tries to
balance fold sizes as evenly as possible.

## References

Brun, P., Thuiller, W., Chauvier, Y., Pellissier, L., Wüest, R. O.,
Wang, Z., & Zimmermann, N. E. (2020). Model complexity affects species
distribution projections under climate change. Journal of Biogeography,
47(1), 130-142.
[doi:10.1111/jbi.13700](https://doi.org/10.1111/jbi.13700)

## See also

[`cv_range`](https://8ginette8.github.io/gbif.range/reference/cv_range.md)()
which uses this function internally to create cross-validation folds.

## Examples

``` r
# \donttest{
# Downloading worldwide all the observations of great panda 
obs.pd <- get_gbif(sp_name = "Ailuropoda melanoleuca")
#> |--------------------------------------------|
#> | Total number (all records)    :        300 |
#> | Kept records                  :         66 |
#> |--------------------------------------------|
#> | Kept records according to parameters:
#> | spatial_issue = FALSE, has_xy = TRUE
#> 
#> ...GBIF records of Ailuropoda melanoleuca: download starting...
#> ------------- #1 (100%..)               
#> 
#> ...Records (XY) filtering summary:
#> ---------------------------------------------
#>                     step removed remaining
#>          Grain filtering       6        60
#>       Duplicated records      13        47
#>          Absence records       0        47
#>          Basis selection      10        37
#>  Establishment selection       0        37
#>               Time frame       0        37
#>        Identical records       0        37
#>         Raster centroids       0        37
#> 
#> Initial records         : 66
#> Total removed           : 29
#> Final records (XY)      : 37
#> ---------------------------------------------
#> Final records (no XY)   : 0

# Create a vector of folds (n = 5) spatially blocked (n = 10)
block.pd <- make_blocks(
    nfolds = 5,
    df = obs.pd[, c("decimalLatitude","decimalLongitude")],
    nblocks = 5
)

# Plot one colour per fold
countries <- terra::vect(
  system.file("extdata", "world_countries.shp", package = "gbif.range")
)
countries.focus <- terra::crop(
    countries,
    terra::ext(73.0, 135.0, 18.0, 54.0)
)
terra::plot(countries.focus, col = "#bcbddc")
graphics::points(
    obs.pd[, c("decimalLongitude","decimalLatitude")],
    pch = 20,
    col = block.pd,
    cex = 1
)


# }
```
